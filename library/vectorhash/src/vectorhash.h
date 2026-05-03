//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_H
#define VECTORHASH_H

#include <stdint.h>

// what SIMD instruction set should be used?
typedef enum { IS_INVALID=-1, IS_SCALAR=0, IS_SSE2, IS_AVX2, IS_AVX512 } is_type;

#ifdef __cplusplus
extern "C" {
#endif

// version where the SIMD support is automatically detected
void VectorHash(const void* buf, size_t len, uint32_t seed, void* out, size_t hash_width);

// version where you set the SIMD instructions by hand -- usually this is not the routine you want!
void VectorHashSIMD(const void* buf, size_t len, uint32_t seed, void* out, is_type SIMDversion, size_t hash_width);

#ifdef __cplusplus
}
#endif

#endif
