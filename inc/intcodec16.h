#ifndef INTCODEC16_H
#define INTCODEC16_H

#include <climits>
#include "intcodec.h"
#include "intrinsics.h"

#include "fixedpoint32.h"

// embedded codec for general integer types using blocks of 16 integers
template <class BitStream, typename Int, typename UInt>
class IntCodec16 : public IntCodec<BitStream, Int, UInt> {
public:
  // encode a sequence of 16 values and write to stream
  static inline void encode(BitStream& bitstream, const Int* data, uint minbits = 0, uint maxbits = UINT_MAX, uint maxprec = UINT_MAX);

  // decode a sequence of 16 values read from stream
  static inline void decode(BitStream& bitstream, Int* data, uint minbits = 0, uint maxbits = UINT_MAX, uint maxprec = UINT_MAX);

  static inline struct RD estimate(const Int* data, uint minbits = 0, uint maxbits = UINT_MAX, uint maxprec = UINT_MAX, uint emax = UINT_MAX);

protected:
  using IntCodec<BitStream, Int, UInt>::encode;
  using IntCodec<BitStream, Int, UInt>::decode;
  using IntCodec<BitStream, Int, UInt>::width;
};

// encode a set of 16 integers
template <class BitStream, typename Int, typename UInt>
void
IntCodec16<BitStream, Int, UInt>::encode(BitStream& bitstream, const Int* data, uint minbits, uint maxbits, uint maxprec)
{
  // compute bit width of each of the 6 groups
  UInt m = 0;
  uint64 w = 0;
  w = (w << 6) + width(data + 13,  3, m);
//  printf ("%llu ", w & 0x3fu);
  w = (w << 6) + width(data + 10,  3, m);
//  printf ("%llu ", w & 0x3fu);
  w = (w << 6) + width(data +  6,  4, m);
//  printf ("%llu ", w & 0x3fu);
  w = (w << 6) + width(data +  3,  3, m);
//  printf ("%llu ", w & 0x3fu);
  w = (w << 6) + width(data +  1,  2, m);
//  printf ("%llu ", w & 0x3fu);
  w = (w << 6) + width(data +  0,  1, m);
//  printf ("%llu\n", w & 0x3fu);
  encode(bitstream, data, minbits, maxbits, maxprec, 0x334321u, w);
}

// decode a set of 16 integers
template <class BitStream, typename Int, typename UInt>
void
IntCodec16<BitStream, Int, UInt>::decode(BitStream& bitstream, Int* data, uint minbits, uint maxbits, uint maxprec)
{
  decode(bitstream, data, minbits, maxbits, maxprec, 0x334321u, 16);
}

// estimate bit-rate and psnr of a set of 16 integers
template <class BitStream, typename Int, typename UInt>
struct RD
IntCodec16<BitStream, Int, UInt>::estimate(const Int* data, uint minbits, uint maxbits, uint maxprec, uint emax)
{
	struct RD rd;
	UInt m = 0;
	uint nb = 0;	// estimated number of bits
	int nbSample = 6;
	int i, j;
	double mse = 0;

	uint *e = (uint*) malloc((nbSample+1)*sizeof(uint));
	uint *l = (uint*) malloc((nbSample+1)*sizeof(uint));

	//initialize group length
	l[0] = 0;
	l[1] = 1;
	l[2] = 3;
	l[3] = 6;
	l[4] = 10;
	l[5] = 13;
	l[6] = 16;

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
	Int data_trun[16];
	Int data_abs[16];
	for (i = 0; i < 16; i++)
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
	for (i = 0; i < 16; i += coeffSampleDistance)
	{
		Int dat = data_trun[i];
		dat = dat << minbitplane;
		FixedPoint::FixedPoint32<3> diff = FixedPoint::FixedPoint32<3>::reinterpret(data_abs[i] - dat);
		float fpdiff = diff.ldexp(emax);
		mse += fpdiff*fpdiff;
	}

	rd.bit_rate = nb;
	rd.mse = mse*coeffSampleDistance;

	return rd;
}


#endif
