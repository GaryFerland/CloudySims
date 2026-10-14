//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_SSE2_H
#define VECTORHASH_SSE2_H

#include "vectorhash_priv.h"

#ifdef VH_INTEL
#include <immintrin.h>
typedef __m128i v4si;
#else
typedef uint32_t v4si;
#endif

void VectorHashBody128_32(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[]);
void VectorHashBody128_64(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[]);
void VectorHashBody128_128(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[]);
void VectorHashBody128_256(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[]);
void VectorHashBody128_512(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[]);
void VectorHashBody128_1024(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[]);

void VectorHash128_32(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash128_64(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash128_128(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash128_256(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash128_512(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash128_1024(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);

#endif
