//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_AVX2_H
#define VECTORHASH_AVX2_H

#include "vectorhash_priv.h"

#ifdef VH_INTEL
#include <immintrin.h>
typedef __m256i v8si;
#else
typedef uint32_t v8si;
#endif

void VectorHashBody256_32(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[]);
void VectorHashBody256_64(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[]);
void VectorHashBody256_128(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[]);
void VectorHashBody256_256(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[]);
void VectorHashBody256_512(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[]);
void VectorHashBody256_1024(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[]);

void VectorHash256_32(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash256_64(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash256_128(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash256_256(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash256_512(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash256_1024(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);

#endif
