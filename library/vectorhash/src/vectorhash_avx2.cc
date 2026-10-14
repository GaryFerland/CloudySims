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
#include "vectorhash_avx2.h"

//-----------------------------------------------------------------------------
// Platform-specific functions and macros

static const size_t nreg256 = vh_nint/8;

#ifdef VH_INTEL
void EXT(VectorHashBody256)(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[])
{
	v8si s[nreg256], x1[nreg256], x2[nreg256];

	VEC( nreg256, s[j] = _mm256_xor_si256(h1[j], h2[j]) );
	VEC( nreg256, s[j] = _mm256_xor_si256(s[j], *data++) );
	VEC( nreg256, x1[j] = _mm256_slli_epi32(h1[j], 11) );
	VEC( nreg256, x2[j] = _mm256_srli_epi32(h1[j], 21) );
	VEC( nreg256, x1[j] = _mm256_or_si256(x1[j], x2[j]) );
	VEC( nreg256, x1[j] = _mm256_xor_si256(x1[j], s[j]) );
	VEC( nreg256, x2[j] = _mm256_slli_epi32(s[j], 14) );
	VEC( nreg256, h1[j] = _mm256_xor_si256(x1[j], x2[j]) );
	VEC( nreg256, x1[j] = _mm256_slli_epi32(s[j], 19) );
	VEC( nreg256, x2[j] = _mm256_srli_epi32(s[j], 13) );
	VEC( nreg256, h2[j] = _mm256_or_si256(x1[j], x2[j]) );

	VEC( nreg256, s[j] = _mm256_xor_si256(h2[j], h3[j]) );
	VEC( nreg256, s[j] = _mm256_xor_si256(s[j], *data++) );
	VEC( nreg256, x1[j] = _mm256_slli_epi32(h2[j], 11) );
	VEC( nreg256, x2[j] = _mm256_srli_epi32(h2[j], 21) );
	VEC( nreg256, x1[j] = _mm256_or_si256(x1[j], x2[j]) );
	VEC( nreg256, x1[j] = _mm256_xor_si256(x1[j], s[j]) );
	VEC( nreg256, x2[j] = _mm256_slli_epi32(s[j], 14) );
	VEC( nreg256, h2[j] = _mm256_xor_si256(x1[j], x2[j]) );
	VEC( nreg256, x1[j] = _mm256_slli_epi32(s[j], 19) );
	VEC( nreg256, x2[j] = _mm256_srli_epi32(s[j], 13) );
	VEC( nreg256, h3[j] = _mm256_or_si256(x1[j], x2[j]) );

	VEC( nreg256, s[j] = _mm256_xor_si256(h3[j], h4[j]) );
	VEC( nreg256, s[j] = _mm256_xor_si256(s[j], *data++) );
	VEC( nreg256, x1[j] = _mm256_slli_epi32(h3[j], 11) );
	VEC( nreg256, x2[j] = _mm256_srli_epi32(h3[j], 21) );
	VEC( nreg256, x1[j] = _mm256_or_si256(x1[j], x2[j]) );
	VEC( nreg256, x1[j] = _mm256_xor_si256(x1[j], s[j]) );
	VEC( nreg256, x2[j] = _mm256_slli_epi32(s[j], 14) );
	VEC( nreg256, h3[j] = _mm256_xor_si256(x1[j], x2[j]) );
	VEC( nreg256, x1[j] = _mm256_slli_epi32(s[j], 19) );
	VEC( nreg256, x2[j] = _mm256_srli_epi32(s[j], 13) );
	VEC( nreg256, h4[j] = _mm256_or_si256(x1[j], x2[j]) );

	VEC( nreg256, s[j] = _mm256_xor_si256(h4[j], h1[j]) );
	VEC( nreg256, s[j] = _mm256_xor_si256(s[j], *data++) );
	VEC( nreg256, x1[j] = _mm256_slli_epi32(h4[j], 11) );
	VEC( nreg256, x2[j] = _mm256_srli_epi32(h4[j], 21) );
	VEC( nreg256, x1[j] = _mm256_or_si256(x1[j], x2[j]) );
	VEC( nreg256, x1[j] = _mm256_xor_si256(x1[j], s[j]) );
	VEC( nreg256, x2[j] = _mm256_slli_epi32(s[j], 14) );
	VEC( nreg256, h4[j] = _mm256_xor_si256(x1[j], x2[j]) );
	VEC( nreg256, x1[j] = _mm256_slli_epi32(s[j], 19) );
	VEC( nreg256, x2[j] = _mm256_srli_epi32(s[j], 13) );
	VEC( nreg256, h1[j] = _mm256_or_si256(x1[j], x2[j]) );
}
#else
void EXT(VectorHashBody256)(const v8si*, v8si[], v8si[], v8si[], v8si[])
{
	(void)0;
}
#endif

void EXT(VectorHash256)(const void* buffer, size_t len, uint32_t seed, void* out, size_t hw)
{
	v8si h1[nreg256], h2[nreg256], h3[nreg256], h4[nreg256];
	stateinit( (uint32_t*)h1, seed, vh_nint );
	stateinit( (uint32_t*)h2, seed, vh_nint );
	stateinit( (uint32_t*)h3, seed, vh_nint );
	stateinit( (uint32_t*)h4, seed, vh_nint );

	size_t nblocks = len/blocksize;
	const v8si* data = (const v8si*)buffer;		
	for( size_t i=0; i < nblocks; i++ )
	{
		EXT(VectorHashBody256)(data, h1, h2, h3, h4);
		data += 4*nreg256;
	}

	// pad the remaining characters and process...
	v8si buf[blocksize/sizeof(v8si)];
	pad_buffer( (const uint8_t*)data, (uint8_t*)buf, len-nblocks*blocksize, blocksize );
	EXT(VectorHashBody256)(buf, h1, h2, h3, h4);

	uint32_t* z1 = (uint32_t*)h1;
	uint32_t* z2 = (uint32_t*)h2;
	uint32_t* z3 = (uint32_t*)h3;
	uint32_t* z4 = (uint32_t*)h4;
	EXT(VectorHashFinalize)(len, z1, z2, z3, z4, out, hw);
}
