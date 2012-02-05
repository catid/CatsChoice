#include "Platform.hpp"
using namespace cat;

#include <iostream>
using namespace std;

/*
	So I found a problem with my favorite PRNG CatsChoice.
	One of the two seed inputs produces predictable outputs
	for the first few outputs, which means I need to search
	for new PRNG parameters for the MWC generators.  These
	generators have large period if:

		A = 32-bit generator parameter.
		(A << 32) - 1 and (A << 31) - 1 are 64-bit primes.

	And the period will be (A << 31) - 1.

	I want to find the largest few values of A that satisfy
	this requirements.  Furthermore, I want A to have a minimal
	number of prime factors; if this is the case, then the
	output will look random in addition to having a large period.

	So I set out to devise a fast 64-bit primality tester.
	The best approach I know of for these small sizes is the
	probabilistic Rabin-Miller primality test.  After a few
	iterations I can conclude the number is probably prime
	and then look at the prime factorization of A to select
	candidates that may go well together.  Wolfram Alpha was
	used to check that the values are really prime.

	I came up with these results:

	-- Candidate 0xffffbe17.  Factors = 3, 1431650141
	-- Candidate 0xffff4b9f.  Factors = 3, 1431640373
	-- Candidate 0xffff0207.  Factors = 3, 1431634093
	-- Candidate 0xfffe1495.  Factors = 3, 1431613831
	-- Candidate 0xfffd8b79.  Factors = 3, 1431602131
	-- Candidate 0xfffd6389.  Factors = 3, 1431598723
	-- Candidate 0xfffd21a7.  Factors = 3, 1431593101
	-- Candidate 0xfffd1361.  Factors = 3, 1431591883
	...

	So I guess that 3 is always a factor of A.

	I then produced CatsChoice-like generators using pairs
	of these A values and tested them with BigCrush.
	These are the pairs that passed the test:

		A1 = 0xffffbe17, A2 = 0xffff4b9f
		A1 = 0xffff0207, A2 = 0xfffe1495
		A1 = 0xfffd8b79, A2 = 0xfffd6389
		A1 = 0xfffd21a7, A2 = 0xfffd1361

	In practice the generator is seeded using two numbers,
	one is fixed and the other increments by 1 each time.
	So the actual generator usage is very different from
	the steady state after seeding.  I then retested the
	above generators with BigCrush by fixing one of the
	seeds = 0 and producing output by incrementing a
	counter for the other seed.  These are the pairs that
	passed this stringent test:

		A1 = 0xffffbe17, A2 = 0xffff4b9f
		A1 = 0xffff0207, A2 = 0xfffe1495
		A1 = 0xfffd8b79, A2 = 0xfffd6389
		A1 = 0xfffd21a7, A2 = 0xfffd1361

	Then I swapped which seed was held fixed and did it
	again.  These are the pairs that passed both tests:

		A1 = 0xffffbe17, A2 = 0xffff4b9f
		A1 = 0xffff0207, A2 = 0xfffe1495
		A1 = 0xfffd8b79, A2 = 0xfffd6389
		A1 = 0xfffd21a7, A2 = 0xfffd1361

	It was a lot of work but I now have a very good and
	fast generator that approximates a hash function and
	can take 32-bit (x,y) coordinate input.

	Say hello to the new CatsChoice PRNG:
*/


/*
	This is a unified implementation of my favorite generator
	that is designed to generate up to 2^^32 numbers per seed.

	Its period is about 2^^126 and passes all BigCrush tests.
	It is the fastest generator I could find that passes all tests.

	Furthermore, the input seeds are hashed to avoid linear
	relationships between the input seeds and the low bits of
	the first few outputs.
*/
class CAT_EXPORT CatsChoice
{
	u64 _x, _y;

public:
	CAT_INLINE void Initialize(u32 seed_x, u32 seed_y)
	{
		// Based on the final mixing function of MurmurHash3
		static const u64 MURMUR_FMIX_K1 = 0xff51afd7ed558ccdULL;
		static const u64 MURMUR_FMIX_K2 = 0xc4ceb9fe1a85ec53ULL;

		// Mix bits of seed x
		u64 x = 0x9368e53c2f6af274ULL ^ seed_x;
		x ^= x >> 33;
		x *= MURMUR_FMIX_K1;
		x ^= x >> 33;
		x *= MURMUR_FMIX_K2;
		x ^= x >> 33;
		_x = x;

		// Mix bits of seed y
		u64 y = 0x586dcd208f7cd3fdULL ^ seed_y;
		y ^= y >> 33;
		y *= MURMUR_FMIX_K1;
		y ^= y >> 33;
		y *= MURMUR_FMIX_K2;
		y ^= y >> 33;
		_y = y;
	}

	CAT_INLINE void Initialize(u32 seed)
	{
		Initialize(seed, seed);
	}

	CAT_INLINE u32 Next()
	{
		_x = (u64)0xffff0207 * (u32)_x + (u32)(_x >> 32);
		_y = (u64)0xffff4b9f * (u32)_y + (u32)(_y >> 32);
		return (u32)_x + (u32)_y;
	}
};



static const u16 PRIME_LIST[] = {
	2,      3,      5,      7,     11,     13,     17,     19,     23,     29, 
	31,     37,     41,     43,     47,     53,     59,     61,     67,     71, 
	73,     79,     83,     89,     97,    101,    103,    107,    109,    113, 
	127,    131,    137,    139,    149,    151,    157,    163,    167,    173, 
	179,    181,    191,    193,    197,    199,    211,    223,    227,    229,
	233,    239,    241,    251,    257,    263,    269,    271,    277,    281, 
	283,    293,    307,    311,    313,    317,    331,    337,    347,    349, 
	353,    359,    367,    373,    379,    383,    389,    397,    401,    409, 
	419,    421,    431,    433,    439,    443,    449,    457,    461,    463, 
	467,    479,    487,    491,    499,    503,    509,    521,    523,    541, 
	547,    557,    563,    569,    571,    577,    587,    593,    599,    601, 
	607,    613,    617,    619,    631,    641,    643,    647,    653,    659, 
	661,    673,    677,    683,    691,    701,    709,    719,    727,    733, 
	739,    743,    751,    757,    761,    769,    773,    787,    797,    809, 
	811,    821,    823,    827,    829,    839,    853,    857,    859,    863, 
	877,    881,    883,    887,    907,    911,    919,    929,    937,    941, 
	947,    953,    967,    971,    977,    983,    991,    997,   1009,   0
};

bool NoSmallPrimeFactors(u32 n)
{
	// Trial divide small primes
	const u16 *lp = PRIME_LIST;
	u16 p;

	while (*lp)
	{
		p = *lp++;

		if (n % p == 0)
			return false;
	}

	return true;
}

// Fairly slow 32-bit prime test
bool IsPrime32(u32 n)
{
	// Trial divide small primes
	const u16 *lp = PRIME_LIST;
	u16 p;

	while (*lp)
	{
		p = *lp++;

		if (n % p == 0)
			return false;
	}

	// Trial divide remaining odd numbers up to square root
	u16 lim = (u16)sqrt((double)n);

	while (p <= lim)
	{
		if (n % p == 0)
			return false;

		p += 2;
	}

	return true;
}

// Very slow 64-bit deterministic prime test
bool IsPrime64(u64 n)
{
	// Trial divide small primes
	const u16 *lp = PRIME_LIST;
	u32 p;

	while (*lp)
	{
		p = *lp++;

		if (n % p == 0)
			return false;
	}

	// Trial divide remaining odd numbers up to 2^32
	while (p != 1)
	{
		if (n % p == 0)
			return false;

		p += 2;
	}

	return true;
}

/*
	Link to big_x64 library
*/
#pragma comment(lib, "big_x64")

extern "C" u64 _mul_mod64(u64 a, u64 b, u64 m);

CAT_INLINE u64 MulMod64(u64 a, u64 b, u64 m)
{
	return _mul_mod64(a, b, m);
}

// returns b ^ e (Mod m)
u64 ExpMod64(u64 b, u64 e, u64 m)
{
	// validate arguments
	if (b == 0 || m <= 1) return 0;
	if (e == 0) return 1;

	// find high bit of exponent
	u64 mask = 0x8000000000000000ULL;
	while ((e & mask) == 0) mask >>= 1;

	// seen 1 set bit, so result = base so far
	u64 r = b;

	// For each bit
	while (mask >>= 1)
	{
		// r = r^2 (mod m)
		r = MulMod64(r, r, m);

		// if exponent bit is set, r = r*b (mod m)
		if (e & mask)
			r = MulMod64(r, b, m);
	}

	return r;
}

// Rabin-Miller method for finding a strong pseudo-prime
bool RabinMillerPrimeTest64(
	CatsChoice &prng,
	const u64 n,	// Number to check for primality
	u32 k)			// Confidence level (40 is pretty good)
{
	// n1 = n - 1
	u64 n1 = n - 1;

	// d = n1
	u64 d = n1;

	// Remove factors of two from d
	while (!(d & 1))
		d >>= 1;

	// Repeat k times
	while (k--)
	{
		u64 a = prng.Next() | ((u64)prng.Next() << 32);
		a %= n;

		// a = a ^ d (Mod n)
		a = ExpMod64(a, d, n);

		u64 t = d;
		while (t != n1 && a != 1 && a != n1)
		{
			// a = a^2 (Mod n), non-critical path
			a = MulMod64(a, a, n);

			// t <<= 1
			t <<= 1;
		}

		if (a != n1 && (t & 1) == 0)
			return false;
	}

	return true;
}

int main()
{
	CatsChoice prng;
	prng.Initialize(0);

	for (int ii = 1; ii; ++ii)
	{
		prng.Initialize(ii, 0);

		cout << (prng.Next() & 1) << endl;
	}

	// Find values that work
	for (u32 a = 0xffffffff; a; a -= 2)
	{
		//if (NoSmallPrimeFactors(a))
		{
			//cout << "Prime A = 0x" << hex << a << endl;

			u64 A32 = ((u64)a << 32) - 1;

			if (RabinMillerPrimeTest64(prng, A32, 10))
			{
				//cout << "- Prime A<<32-1 = 0x" << hex << a << endl;

				u64 A31 = ((u64)a << 31) - 1;

				if (RabinMillerPrimeTest64(prng, A31, 10))
				{
					int factors[2];
					int count = 0;

					u32 n = a;
					for (u32 p = 3; p <= n; p += 2)
					{
						while (n % p == 0)
						{
							if (count == 2)
							{
								count = -1;
								p = n;
								break;
							}

							factors[count++] = p;

							n /= p;
						}
					}

					if (count > 0)
					{
						if (count != 2)
						{
							cout << "WOW!" << endl;
						}

						cout << "-- Candidate 0x" << hex << a << dec << ".  Factors = " << factors[0] << ", " << factors[1] << endl;
					}
				}
			}
		}
	}

	cin.get();

	return 0;
}
