#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include "zfpcompress.h"
#include "sz.h"
#include "mpi.h"

int main(int argc, char* argv[])
{
	// default settings
	bool dp = false;
	uint nx = 0;
	uint ny = 0;
	uint nz = 0;
	uint minbits = 0;
	uint maxbits = 0;
	uint maxprec = 0;
	int minexp = INT_MIN;
	char* inpath = 0;
	char zippath[640];
	char* outpath = 0;
	int status;

	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


	// parse arguments
	switch (argc) {
		case 11:
			outpath = argv[10];
			// FALLTHROUGH
		case 10:
			inpath = argv[9];
			// FALLTHROUGH
		case 9:
			if (sscanf(argv[8], "%d", &minexp) != 1)
				goto usage;
			// FALLTHROUGH
		case 8:
			if (sscanf(argv[7], "%u", &maxprec) != 1)
				goto usage;
			// FALLTHROUGH
		case 7:
			if (sscanf(argv[6], "%u", &maxbits) != 1)
				goto usage;
			// FALLTHROUGH
		case 6:
			if (sscanf(argv[5], "%u", &minbits) != 1)
				goto usage;
			// FALLTHROUGH
		case 5:
			if (sscanf(argv[4], "%u", &nz) != 1)
				goto usage;
			// FALLTHROUGH
		case 4:
			if (sscanf(argv[3], "%u", &ny) != 1)
				goto usage;
			// FALLTHROUGH
		case 3:
			if (sscanf(argv[2], "%u", &nx) != 1)
				goto usage;
			if (argv[1] == std::string("-f"))
				dp = false;
			else if (argv[1] == std::string("-d"))
				dp = true;
			else
				goto usage;
			break;
		default:
usage:
			std::cerr << "Usage: zfp <-f|-d> <nx> [ny nz minbits maxbits maxprec minexp infile outfile]" << std::endl;
			std::cerr << "-f : single precision (float type)" << std::endl;
			std::cerr << "-d : double precision (double type)" << std::endl;
			std::cerr << "nx, ny, nz : grid dimensions (set nz = 0 for 2D, ny = nz = 0 for 1D)" << std::endl;
			std::cerr << "minbits : min # bits per 4^d values in d dimensions (= maxbits for fixed rate)" << std::endl;
			std::cerr << "maxbits : max # bits per 4^d values in d dimensions" << std::endl;
			std::cerr << "maxprec : max # bits of precision per value (0 for fixed rate)" << std::endl;
			std::cerr << "minexp : min bit plane coded (error tolerance = 2^minexp; -1024 for fixed rate)" << std:: endl;
			std::cerr << "infile : optional floating-point input file to compress" << std::endl;
			std::cerr << "outfile : optional output file for reconstructed data" << std::endl;
			std::cerr << "Examples:" << std::endl;
			std::cerr << "zfp -f 100 100 100 1024 1024 : 2x fixed-rate compression of 100x100x100 floats" << std::endl;
			std::cerr << "zfp -d 1000000 0 0 128 128 : 2x fixed-rate compression of stream of 1M doubles" << std::endl;
			std::cerr << "zfp -d 1000 1000 0 0 0 32 : 32-bit precision compression of 1000x1000 doubles" << std::endl;
			std::cerr << "zfp -d 1000000 0 0 0 0 0 -16 : compression of 1M doubles with < 2^-16 error" << std::endl;
			return EXIT_FAILURE;
	}

	// effective array dimensions
	uint mx = std::max(nx, 1u);
	uint my = std::max(ny, 1u);
	uint mz = std::max(nz, 1u);

	// array dimensionality
	uint dims = 3;
	if (nz == 0) {
		if (ny == 0) {
			if (nx == 0) {
				std::cerr << "cannot compress zero-size array" << std::endl;
				return EXIT_FAILURE;
			}
			else
				dims = 1;
		}
		else
			dims = 2;
	}
	else
		dims = 3;

	// number of floating-point values per block
	uint blocksize = 1u << (2 * dims);

	// number of blocks
	uint blocks = ((mx + 3) / 4) * ((my + 3) / 4) * ((mz + 3) / 4);

	// size of floating-point type in bytes
	size_t typesize = dp ? sizeof(double) : sizeof(float);

	// correct compression parameters if zero initialized
	if (maxbits == 0)
		maxbits = blocksize * CHAR_BIT * typesize;
	if (maxprec == 0)
		maxprec = CHAR_BIT * typesize;

	// allocate space for uncompressed and compressed fields
	size_t insize = mx * my * mz * typesize;
	void* f = new unsigned char[insize];
	void* g = new unsigned char[insize];
	size_t outsize = (blocks * std::min(maxbits, blocksize * maxprec) + CHAR_BIT - 1) / CHAR_BIT;
	outsize = std::max(outsize, 2 * insize);
	unsigned char* zip = new unsigned char[outsize];

	double start, end;
	double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

	if (world_rank == 0) printf ("Start ... \n");

	// initialize uncompressed field
	// read from file
	start = MPI_Wtime();
	FILE* file = fopen(inpath, "rb");
	if (!file) {
		std::cerr << "cannot open file" << std::endl;
		return EXIT_FAILURE;
	}
	if (fread(f, typesize, mx * my * mz, file) != mx * my * mz) {
		std::cerr << "cannot read file" << std::endl;
		return EXIT_FAILURE;
	}
	fclose(file);
	end = MPI_Wtime();
	costReadOri += end - start;
	MPI_Barrier(MPI_COMM_WORLD);


	/* calculate value-range-based relative error bound */
	unsigned int i = 0;
	double relEB = pow(2, minexp);
	double absEB = 0;
	double range;

	if (!dp)
	{
		float *data = (float*)malloc(mx*my*mz*sizeof(float));
		memcpy(data, f, insize);
		float maxval, minval;
		maxval = data[0];
		minval = data[0];
		for (i = 0; i < mx*my*mz; i++)
		{
			if (data[i] > maxval) maxval = data[i];
			if (data[i] < minval) minval = data[i];
		}
		range = maxval - minval;
		absEB = relEB*range;
		free(data);
	}
	else
	{
		double *data = (double*)malloc(mx*my*mz*sizeof(double));
		memcpy(data, f, insize);
		double maxval, minval;
		maxval = data[0];
		minval = data[0];
		for (i = 0; i < mx*my*mz; i++)
		{
			if (data[i] > maxval) maxval = data[i];
			if (data[i] < minval) minval = data[i];
		}
		range = maxval - minval;
		absEB = relEB*range;
		free(data);
	}

	frexp(absEB, &minexp);

	// ZFP compress data
	start = MPI_Wtime();
	outsize = ZFP::compress(f, zip, dp, nx, ny, nz, minbits, maxbits, maxprec, minexp);
	end = MPI_Wtime();
	costComp += end - start;
	MPI_Barrier(MPI_COMM_WORLD);

	// SZ compress data
	//  SZ_Init("sz.config");
	//  size_t outsize_sz;
	//  if (!dp)
	//  {
	//	  unsigned char *bytes = SZ_compress_args(SZ_FLOAT, (float*)f, &outsize_sz, REL, 0.0001, 0.0001, 0.0001, 1, 0, 0, nz, ny, nx);
	//	  free (bytes);
	//  }

	// ZFP decompress data
	start = MPI_Wtime();
	sprintf (zippath, "%s.zfp", inpath);
	writeByteData(zip, outsize, zippath, &status);
	end = MPI_Wtime();
	costWriteZip += end - start;
	MPI_Barrier(MPI_COMM_WORLD);

#if IT_SEEMS_TOO_GOOD_TO_BE_TRUE
	// for skeptics: relocate compressed data
	unsigned char* copy = new unsigned char[outsize];
	std::copy(zip, zip + outsize, copy);
	delete[] zip;
	zip = copy;
#endif

	// read compressed data
	start = MPI_Wtime();
	size_t insize2;	
	unsigned char *zip2 = readByteData(zippath, &insize2, &status);	
	end = MPI_Wtime();
	costReadZip += end - start;
	MPI_Barrier(MPI_COMM_WORLD);

	// decompress data
	start = MPI_Wtime();  
	ZFP::decompress(zip2, g, dp, nx, ny, nz, minbits, maxbits, maxprec, minexp);
	end = MPI_Wtime();
	costDecomp += end - start;
	MPI_Barrier(MPI_COMM_WORLD);

	// write reconstructed data
	start = MPI_Wtime();	  
	file = fopen(outpath, "wb");
	fwrite(g, typesize, mx * my * mz, file);
	fclose(file);
	end = MPI_Wtime();
	costWriteOut += end - start;
	MPI_Barrier(MPI_COMM_WORLD);

	double globalcostReadOri, globalcostReadZip, globalcostWriteZip, globalcostWriteOut, globalcostComp, globalcostDecomp;

	MPI_Reduce(&costReadOri, &globalcostReadOri, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costReadZip, &globalcostReadZip, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costWriteZip, &globalcostWriteZip, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costWriteOut, &globalcostWriteOut, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costComp, &globalcostComp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costDecomp, &globalcostDecomp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (world_rank == 0)
	{
		printf ("Finish.\n");
		printf ("Timecost to read original files = %.2f seconds\n", globalcostReadOri/world_size);
		printf ("Timecost to read compressed files = %.2f seconds\n", globalcostReadZip/world_size);
		printf ("Timecost to write compressed files = %.2f seconds\n", globalcostWriteZip/world_size);
		printf ("Timecost to write decompressed files = %.2f seconds\n", globalcostWriteOut/world_size);
		printf ("Timecost of compression using %d processes = %.2f seconds\n", world_size, globalcostComp/world_size);
		printf ("Timecost of decompression using %d processes = %.2f seconds\n\n", world_size, globalcostDecomp/world_size);
	}

	// clean up
	delete[] static_cast<unsigned char*>(f);
	delete[] static_cast<unsigned char*>(g);
	delete[] zip;
	delete[] zip2;

	MPI_Finalize();

	return EXIT_SUCCESS;
}
