//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_FINALIZE_H
#define VECTORHASH_FINALIZE_H

void VectorHashFinalize_32(size_t len, uint32_t* h1, uint32_t* h2, uint32_t* h3, uint32_t* h4, void* out, size_t hw);
void VectorHashFinalize_64(size_t len, uint32_t* h1, uint32_t* h2, uint32_t* h3, uint32_t* h4, void* out, size_t hw);
void VectorHashFinalize_128(size_t len, uint32_t* h1, uint32_t* h2, uint32_t* h3, uint32_t* h4, void* out, size_t hw);
void VectorHashFinalize_256(size_t len, uint32_t* h1, uint32_t* h2, uint32_t* h3, uint32_t* h4, void* out, size_t hw);
void VectorHashFinalize_512(size_t len, uint32_t* h1, uint32_t* h2, uint32_t* h3, uint32_t* h4, void* out, size_t hw);
void VectorHashFinalize_1024(size_t len, uint32_t* h1, uint32_t* h2, uint32_t* h3, uint32_t* h4, void* out, size_t hw);

#endif
