#ifndef INTCODEC64_H
#define INTCODEC64_H

#include <climits>
#include "intcodec.h"

#include "fixedpoint32.h"

// embedded codec for general integer types using blocks of 64 integers
template <class BitStream, typename Int, typename UInt>
class IntCodec64 : public IntCodec<BitStream, Int, UInt> {
public:
  // encode a sequence of 64 values and write to stream
  static inline void encode(BitStream& bitstream, const Int* data, uint minbits = 0, uint maxbits = UINT_MAX, uint maxprec = UINT_MAX);

  // decode a sequence of 64 values read from stream
  static inline void decode(BitStream& bitstream, Int* data, uint minbits = 0, uint maxbits = UINT_MAX, uint maxprec = UINT_MAX);

  static inline struct RD estimate(const Int* data, uint minbits = 0, uint maxbits = UINT_MAX, uint maxprec = UINT_MAX, uint emax = UINT_MAX);

protected:
  using IntCodec<BitStream, Int, UInt>::encode;
  using IntCodec<BitStream, Int, UInt>::decode;
  using IntCodec<BitStream, Int, UInt>::width;
};

// encode a set of 64 integers
template <class BitStream, typename Int, typename UInt>
void
IntCodec64<BitStream, Int, UInt>::encode(BitStream& bitstream, const Int* data, uint minbits, uint maxbits, uint maxprec)
{
  // compute bit width of each of the 9 groups
  UInt m = 0;
  uint64 w = 0;
  w = (w << 6) + width(data + 60,  4, m);
  w = (w << 6) + width(data + 54,  6, m);
  w = (w << 6) + width(data + 44, 10, m);
  w = (w << 6) + width(data + 32, 12, m);
  w = (w << 6) + width(data + 20, 12, m);
  w = (w << 6) + width(data + 10, 10, m);
  w = (w << 6) + width(data +  4,  6, m);
  w = (w << 6) + width(data +  1,  3, m);
  w = (w << 6) + width(data +  0,  1, m);
  encode(bitstream, data, minbits, maxbits, maxprec, 0x46acca631ull, w);
}

// decode a set of 64 integers
template <class BitStream, typename Int, typename UInt>
void
IntCodec64<BitStream, Int, UInt>::decode(BitStream& bitstream, Int* data, uint minbits, uint maxbits, uint maxprec)
{
  decode(bitstream, data, minbits, maxbits, maxprec, 0x46acca631ull, 64);
}


// estimate bit-rate and psnr of a set of 64 integers
template <class BitStream, typename Int, typename UInt>
struct RD
IntCodec64<BitStream, Int, UInt>::estimate(const Int* data, uint minbits, uint maxbits, uint maxprec, uint emax)
{
	struct RD rd;
	UInt m = 0;
	uint nb = 0;
	int nbSample = 9;
	int i, j;
	double mse = 0;

	uint *e = (uint*) malloc((nbSample+1)*sizeof(uint));
	uint *l = (uint*) malloc((nbSample+1)*sizeof(uint));

	//initialize group length
	l[0] = 0;
	l[1] = 1;
	l[2] = 4;
	l[3] = 10;
	l[4] = 20;
	l[5] = 32;
	l[6] = 44;
	l[7] = 54;
	l[8] = 60;
	l[9] = 64;

	//estimate bit-rate
	// count number of significant bits for sampling coefficients
	for (i = nbSample-1; i >= 0; i--)
		e[i] = width(data+l[i], l[i+1]-l[i], m);

	// estimate number of significant bits to be encoded
	const uint intprec = CHAR_BIT * sizeof(Int);
	uint minbitplane = intprec > maxprec ? intprec - maxprec : 0;

	for (i = 0; i < nbSample; i++)
		nb += (l[i+1]-l[i]) * (e[i] < minbitplane ? 0 : (e[i] - minbitplane) );

	// count bits for group tests
	for (i = intprec, j = 0; i-- > minbitplane;)
	{
		while(j < nbSample)
		{
			if (e[j] > i)
			{
				j++;
				nb++;
			}
			else
			{
				nb++;
				break;
			}
		}
	}

	// count bits for sign
	Int data_trun[64];
	Int data_abs[64];
	for (i = 0; i < 64; i++)
	{
		Int dat = data[i];
		dat = (dat < 0 ? -dat : +dat);
		data_abs[i] = dat;
		dat = dat >> minbitplane;
		if (dat >= 1) nb++;
		data_trun[i] = dat;
	}

	//estimate MSE
	uint coeffSampleDistance = 1;
	for (i = 0; i < 64; i += coeffSampleDistance)
	{
		Int dat = data_trun[i];
		dat = dat << minbitplane;
		FixedPoint::FixedPoint32<4> diff = FixedPoint::FixedPoint32<4>::reinterpret(data_abs[i] - dat);
		float fpdiff = diff.ldexp(emax);
		mse += fpdiff*fpdiff;
	}

	rd.bit_rate = nb;
	rd.mse = mse*coeffSampleDistance;

	return rd;
}

#endif
