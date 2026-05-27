//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#include "vectorhash_priv.h"
#include "vectorhash_core.h"
#include "vectorhash_finalize.h"
#include "vectorhash_sse2.h"

//-----------------------------------------------------------------------------
// Platform-specific functions and macros

static const size_t nreg128 = vh_nint/4;

#ifdef VH_INTEL
void EXT(VectorHashBody128)(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[])
{
	v4si s[nreg128], x1[nreg128], x2[nreg128];

	VEC( nreg128, s[j] = _mm_xor_si128(h1[j], h2[j]) );
	VEC( nreg128, s[j] = _mm_xor_si128(s[j], *data++) );
	VEC( nreg128, x1[j] = _mm_slli_epi32(h1[j], 11) );
	VEC( nreg128, x2[j] = _mm_srli_epi32(h1[j], 21) );
	VEC( nreg128, x1[j] = _mm_or_si128(x1[j], x2[j]) );
	VEC( nreg128, x1[j] = _mm_xor_si128(x1[j], s[j]) );
	VEC( nreg128, x2[j] = _mm_slli_epi32(s[j], 14) );
	VEC( nreg128, h1[j] = _mm_xor_si128(x1[j], x2[j]) );
	VEC( nreg128, x1[j] = _mm_slli_epi32(s[j], 19) );
	VEC( nreg128, x2[j] = _mm_srli_epi32(s[j], 13) );
	VEC( nreg128, h2[j] = _mm_or_si128(x1[j], x2[j]) );

	VEC( nreg128, s[j] = _mm_xor_si128(h2[j], h3[j]) );
	VEC( nreg128, s[j] = _mm_xor_si128(s[j], *data++) );
	VEC( nreg128, x1[j] = _mm_slli_epi32(h2[j], 11) );
	VEC( nreg128, x2[j] = _mm_srli_epi32(h2[j], 21) );
	VEC( nreg128, x1[j] = _mm_or_si128(x1[j], x2[j]) );
	VEC( nreg128, x1[j] = _mm_xor_si128(x1[j], s[j]) );
	VEC( nreg128, x2[j] = _mm_slli_epi32(s[j], 14) );
	VEC( nreg128, h2[j] = _mm_xor_si128(x1[j], x2[j]) );
	VEC( nreg128, x1[j] = _mm_slli_epi32(s[j], 19) );
	VEC( nreg128, x2[j] = _mm_srli_epi32(s[j], 13) );
	VEC( nreg128, h3[j] = _mm_or_si128(x1[j], x2[j]) );

	VEC( nreg128, s[j] = _mm_xor_si128(h3[j], h4[j]) );
	VEC( nreg128, s[j] = _mm_xor_si128(s[j], *data++) );
	VEC( nreg128, x1[j] = _mm_slli_epi32(h3[j], 11) );
	VEC( nreg128, x2[j] = _mm_srli_epi32(h3[j], 21) );
	VEC( nreg128, x1[j] = _mm_or_si128(x1[j], x2[j]) );
	VEC( nreg128, x1[j] = _mm_xor_si128(x1[j], s[j]) );
	VEC( nreg128, x2[j] = _mm_slli_epi32(s[j], 14) );
	VEC( nreg128, h3[j] = _mm_xor_si128(x1[j], x2[j]) );
	VEC( nreg128, x1[j] = _mm_slli_epi32(s[j], 19) );
	VEC( nreg128, x2[j] = _mm_srli_epi32(s[j], 13) );
	VEC( nreg128, h4[j] = _mm_or_si128(x1[j], x2[j]) );

	VEC( nreg128, s[j] = _mm_xor_si128(h4[j], h1[j]) );
	VEC( nreg128, s[j] = _mm_xor_si128(s[j], *data++) );
	VEC( nreg128, x1[j] = _mm_slli_epi32(h4[j], 11) );
	VEC( nreg128, x2[j] = _mm_srli_epi32(h4[j], 21) );
	VEC( nreg128, x1[j] = _mm_or_si128(x1[j], x2[j]) );
	VEC( nreg128, x1[j] = _mm_xor_si128(x1[j], s[j]) );
	VEC( nreg128, x2[j] = _mm_slli_epi32(s[j], 14) );
	VEC( nreg128, h4[j] = _mm_xor_si128(x1[j], x2[j]) );
	VEC( nreg128, x1[j] = _mm_slli_epi32(s[j], 19) );
	VEC( nreg128, x2[j] = _mm_srli_epi32(s[j], 13) );
	VEC( nreg128, h1[j] = _mm_or_si128(x1[j], x2[j]) );
}
#else
void EXT(VectorHashBody128)(const v4si*, v4si[], v4si[], v4si[], v4si[])
{
	(void)0;
}
#endif

void EXT(VectorHash128)(const void* buffer, size_t len, uint32_t seed, void* out, size_t hw)
{
	v4si h1[nreg128], h2[nreg128], h3[nreg128], h4[nreg128];
	stateinit( (uint32_t*)h1, seed, vh_nint );
	stateinit( (uint32_t*)h2, seed, vh_nint );
	stateinit( (uint32_t*)h3, seed, vh_nint );
	stateinit( (uint32_t*)h4, seed, vh_nint );

	size_t nblocks = len/blocksize;
	const v4si* data = (const v4si*)buffer;
	for( size_t i=0; i < nblocks; i++ )
	{
		EXT(VectorHashBody128)(data, h1, h2, h3, h4);
		data += 4*nreg128;
	}

	// pad the remaining characters and process...
	v4si buf[blocksize/sizeof(v4si)];
	pad_buffer( (const uint8_t*)data, (uint8_t*)buf, len-nblocks*blocksize, blocksize );
	EXT(VectorHashBody128)(buf, h1, h2, h3, h4);

	uint32_t* z1 = (uint32_t*)h1;
	uint32_t* z2 = (uint32_t*)h2;
	uint32_t* z3 = (uint32_t*)h3;
	uint32_t* z4 = (uint32_t*)h4;
	EXT(VectorHashFinalize)(len, z1, z2, z3, z4, out, hw);
}
