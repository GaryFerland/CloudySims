//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_SCALAR_H
#define VECTORHASH_SCALAR_H

#include "vectorhash_priv.h"

void VectorHashBody32_32(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[]);
void VectorHashBody32_64(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[]);
void VectorHashBody32_128(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[]);
void VectorHashBody32_256(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[]);
void VectorHashBody32_512(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[]);
void VectorHashBody32_1024(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[]);

void VectorHash32_32(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash32_64(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash32_128(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash32_256(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash32_512(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);
void VectorHash32_1024(const void* buf, size_t len, uint32_t seed, void* out, size_t hw);

#endif
