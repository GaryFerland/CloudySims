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
#include "vectorhash_scalar.h"

void EXT(VectorHashBody32)(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[])
{
	uint32_t s[vh_nint];

	VEC( vh_nint, s[j] = h1[j] ^ h2[j] ^ getblock32( *data++ ) );
	VEC( vh_nint, h1[j] = ROTL32(h1[j], 11) ^ s[j] ^ (s[j] << 14) );
	VEC( vh_nint, h2[j] = ROTL32(s[j], 19) );

	VEC( vh_nint, s[j] = h2[j] ^ h3[j] ^ getblock32( *data++ ) );
	VEC( vh_nint, h2[j] = ROTL32(h2[j], 11) ^ s[j] ^ (s[j] << 14) );
	VEC( vh_nint, h3[j] = ROTL32(s[j], 19) );

	VEC( vh_nint, s[j] = h3[j] ^ h4[j] ^ getblock32( *data++ ) );
	VEC( vh_nint, h3[j] = ROTL32(h3[j], 11) ^ s[j] ^ (s[j] << 14) );
	VEC( vh_nint, h4[j] = ROTL32(s[j], 19) );

	VEC( vh_nint, s[j] = h4[j] ^ h1[j] ^ getblock32( *data++ ) );
	VEC( vh_nint, h4[j] = ROTL32(h4[j], 11) ^ s[j] ^ (s[j] << 14) );
	VEC( vh_nint, h1[j] = ROTL32(s[j], 19) );
}

void EXT(VectorHash32)(const void* buffer, size_t len, uint32_t seed, void* out, size_t hw)
{
	uint32_t h1[vh_nint], h2[vh_nint], h3[vh_nint], h4[vh_nint];
	stateinit( h1, seed, vh_nint );
	stateinit( h2, seed, vh_nint );
	stateinit( h3, seed, vh_nint );
	stateinit( h4, seed, vh_nint );

	size_t nblocks = len/blocksize;
	const uint32_t* data = (const uint32_t*)buffer;
	for( size_t i=0; i < nblocks; i++ )
	{
		EXT(VectorHashBody32)(data, h1, h2, h3, h4);
		data += 4*vh_nint;
	}

	// pad the remaining characters and process...
	uint32_t buf[blocksize/sizeof(uint32_t)];
	pad_buffer( (const uint8_t*)data, (uint8_t*)buf, len-nblocks*blocksize, blocksize );
	EXT(VectorHashBody32)(buf, h1, h2, h3, h4);

	EXT(VectorHashFinalize)(len, h1, h2, h3, h4, out, hw);
}
