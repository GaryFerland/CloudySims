//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_CORE_H
#define VECTORHASH_CORE_H

#include <cstdint>
#include "vectorhash.h"
#include "vectorhash_priv.h"
#include "vectorhash_avx512.h"
#include "vectorhash_avx2.h"
#include "vectorhash_sse2.h"
#include "vectorhash_scalar.h"

void stateinit(uint32_t st[], uint32_t& seed, size_t lvh_nint);
is_type GetSIMDVersion();
void VectorHashBody32(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[],
					  size_t hash_width);
void VectorHashBody128(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[], size_t hash_width);
void VectorHashBody256(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[], size_t hash_width);
void VectorHashBody512(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[], size_t hash_width);
void VectorHash(const void* buf, size_t len, uint32_t seed, void* out, is_type SIMDversion, size_t hash_width);

// This routine is needed because the standard says that integer overflow results in undefined behavior.
// This routine looks like a lot of overhead, but a good compiler will optimize this into a single
// hardware instruction provided the hardware does the right thing (which is usually the case).
inline uint32_t vh_add(uint32_t a, uint32_t b)
{
	uint64_t c = a;
	c += uint64_t(b);
	return uint32_t(c & UINT64_C(0xffffffff));
}

inline uint32_t fmix32(uint32_t h0, uint32_t h1 = 0xd86b048b)
{
	h1 ^= h0;
	h0 = ROTL32(h0, 11) ^ h1 ^ ROTL32(h1, 13);
	h1 = ROTL32(h1, 19);
	return vh_add(h0, h1);
}

inline uint32_t pow2roundup(uint32_t x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}

#endif
