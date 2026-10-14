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
#include "vectorhash_avx512.h"

//-----------------------------------------------------------------------------
// Platform-specific functions and macros

static const size_t nreg512 = vh_nint/16;

#ifdef VH_INTEL
void EXT(VectorHashBody512)(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[])
{
	v16si s[nreg512], x1[nreg512], x2[nreg512];

	VEC( nreg512, s[j] = _mm512_xor_si512(h1[j], h2[j]) );
	VEC( nreg512, s[j] = _mm512_xor_si512(s[j], *data++) );
	VEC( nreg512, x1[j] = _mm512_rol_epi32(h1[j], 11) );
	VEC( nreg512, x1[j] = _mm512_xor_si512(x1[j], s[j]) );
	VEC( nreg512, x2[j] = _mm512_slli_epi32(s[j], 14) );
	VEC( nreg512, h1[j] = _mm512_xor_si512(x1[j], x2[j]) );
	VEC( nreg512, h2[j] = _mm512_rol_epi32(s[j], 19) );

	VEC( nreg512, s[j] = _mm512_xor_si512(h2[j], h3[j]) );
	VEC( nreg512, s[j] = _mm512_xor_si512(s[j], *data++) );
	VEC( nreg512, x1[j] = _mm512_rol_epi32(h2[j], 11) );
	VEC( nreg512, x1[j] = _mm512_xor_si512(x1[j], s[j]) );
	VEC( nreg512, x2[j] = _mm512_slli_epi32(s[j], 14) );
	VEC( nreg512, h2[j] = _mm512_xor_si512(x1[j], x2[j]) );
	VEC( nreg512, h3[j] = _mm512_rol_epi32(s[j], 19) );

	VEC( nreg512, s[j] = _mm512_xor_si512(h3[j], h4[j]) );
	VEC( nreg512, s[j] = _mm512_xor_si512(s[j], *data++) );
	VEC( nreg512, x1[j] = _mm512_rol_epi32(h3[j], 11) );
	VEC( nreg512, x1[j] = _mm512_xor_si512(x1[j], s[j]) );
	VEC( nreg512, x2[j] = _mm512_slli_epi32(s[j], 14) );
	VEC( nreg512, h3[j] = _mm512_xor_si512(x1[j], x2[j]) );
	VEC( nreg512, h4[j] = _mm512_rol_epi32(s[j], 19) );

	VEC( nreg512, s[j] = _mm512_xor_si512(h4[j], h1[j]) );
	VEC( nreg512, s[j] = _mm512_xor_si512(s[j], *data++) );
	VEC( nreg512, x1[j] = _mm512_rol_epi32(h4[j], 11) );
	VEC( nreg512, x1[j] = _mm512_xor_si512(x1[j], s[j]) );
	VEC( nreg512, x2[j] = _mm512_slli_epi32(s[j], 14) );
	VEC( nreg512, h4[j] = _mm512_xor_si512(x1[j], x2[j]) );
	VEC( nreg512, h1[j] = _mm512_rol_epi32(s[j], 19) );
}
#else
void EXT(VectorHashBody512)(const v16si*, v16si[], v16si[], v16si[], v16si[])
{
	(void)0;
}
#endif

void EXT(VectorHash512)(const void* buffer, size_t len, uint32_t seed, void* out, size_t hw)
{
	v16si h1[nreg512], h2[nreg512], h3[nreg512], h4[nreg512];
	stateinit( (uint32_t*)h1, seed, vh_nint );
	stateinit( (uint32_t*)h2, seed, vh_nint );
	stateinit( (uint32_t*)h3, seed, vh_nint );
	stateinit( (uint32_t*)h4, seed, vh_nint );

	size_t nblocks = len/blocksize;
	const v16si* data = (const v16si*)buffer;
	for( size_t i=0; i < nblocks; i++ )
	{
		EXT(VectorHashBody512)(data, h1, h2, h3, h4);
		data += 4*nreg512;
	}

	// pad the remaining characters and process...
	v16si buf[blocksize/sizeof(v16si)];
	pad_buffer( (const uint8_t*)data, (uint8_t*)buf, len-nblocks*blocksize, blocksize );
	EXT(VectorHashBody512)(buf, h1, h2, h3, h4);

	uint32_t* z1 = (uint32_t*)h1;
	uint32_t* z2 = (uint32_t*)h2;
	uint32_t* z3 = (uint32_t*)h3;
	uint32_t* z4 = (uint32_t*)h4;
	EXT(VectorHashFinalize)(len, z1, z2, z3, z4, out, hw);
}
