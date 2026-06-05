//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_AVX512_H
#define VECTORHASH_AVX512_H

#include "vectorhash_priv.h"

#ifdef VH_INTEL
#include <immintrin.h>
typedef __m512i v16si;
#else
typedef uint32_t v16si;
#endif

void VectorHashBody512_32(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[]);
void VectorHashBody512_64(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[]);
void VectorHashBody512_128(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[]);
void VectorHashBody512_256(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[]);
void VectorHashBody512_512(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[]);
void VectorHashBody512_1024(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[]);

void VectorHash512_32(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash512_64(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash512_128(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash512_256(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash512_512(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash512_1024(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);

#endif
