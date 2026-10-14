//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstdlib>

#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#else
#define _POSIX_MAPPED_FILES 0
#endif

#if _POSIX_MAPPED_FILES > 0
#include <sys/mman.h>
#endif

#include <fcntl.h>

#include "../cpuid/cpuinfo.hpp"
#include "vectorhash.h"
#include "vectorhash_priv.h"
#include "vectorhash_binary.h"
#include "vectorhash_core.h"
#include "vectorhash_finalize.h"
#include "vectorhash_avx512.h"
#include "vectorhash_avx2.h"
#include "vectorhash_sse2.h"
#include "vectorhash_scalar.h"

void stateinit(uint32_t st[], uint32_t& seed, size_t lvh_nint)
{
	for( size_t i=0; i < lvh_nint; i++ )
		st[i] = seed;
	seed = fmix32(seed);
}

API_EXPORT is_type GetSIMDVersion()
{
#ifdef VH_INTEL
	// determine which SIMD instruction set should be used to calculate the checksum
	//
	// -1: not initialized yet
	//  0: use scalar version of the algorithm
	//  1: use SSE2 instructions
	//  2: use AVX2 instructions
	//  3: use AVX512F instructions

	cpuid::cpuinfo m_cpuid;
	if( m_cpuid.has_avx512_f() )
		return IS_AVX512;
	else if( m_cpuid.has_avx2() )
		return IS_AVX2;
	else if( m_cpuid.has_sse2() )
		return IS_SSE2;
	else
		return IS_SCALAR;
#else
	return IS_SCALAR;
#endif
}

void VectorHashBody32(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[], size_t hw)
{
	uint32_t rhw =  pow2roundup(hw);
	if( rhw == 32 )
		VectorHashBody32_32(data, h1, h2, h3, h4);
	else if( rhw == 64 )
		VectorHashBody32_64(data, h1, h2, h3, h4);
	else if( rhw == 128 )
		VectorHashBody32_128(data, h1, h2, h3, h4);
	else if( rhw == 256 )
		VectorHashBody32_256(data, h1, h2, h3, h4);
	else if( rhw == 512 )
		VectorHashBody32_512(data, h1, h2, h3, h4);
	else if( rhw == 1024 )
		VectorHashBody32_1024(data, h1, h2, h3, h4);
	else
	{
		cout << "Internal error: impossible value for rounded hash width: " << rhw << "." << endl;
		exit(1);
	}
}

void VectorHashBody128(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[], size_t hw)
{
	uint32_t rhw =  pow2roundup(hw);
	if( rhw == 32 )
		VectorHashBody128_32(data, h1, h2, h3, h4);
	else if( rhw == 64 )
		VectorHashBody128_64(data, h1, h2, h3, h4);
	else if( rhw == 128 )
		VectorHashBody128_128(data, h1, h2, h3, h4);
	else if( rhw == 256 )
		VectorHashBody128_256(data, h1, h2, h3, h4);
	else if( rhw == 512 )
		VectorHashBody128_512(data, h1, h2, h3, h4);
	else if( rhw == 1024 )
		VectorHashBody128_1024(data, h1, h2, h3, h4);
	else
	{
		cout << "Internal error: impossible value for rounded hash width: " << rhw << "." << endl;
		exit(1);
	}
}

void VectorHashBody256(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[], size_t hw)
{
	uint32_t rhw =  pow2roundup(hw);
	if( rhw == 32 )
		VectorHashBody256_32(data, h1, h2, h3, h4);
	else if( rhw == 64 )
		VectorHashBody256_64(data, h1, h2, h3, h4);
	else if( rhw == 128 )
		VectorHashBody256_128(data, h1, h2, h3, h4);
	else if( rhw == 256 )
		VectorHashBody256_256(data, h1, h2, h3, h4);
	else if( rhw == 512 )
		VectorHashBody256_512(data, h1, h2, h3, h4);
	else if( rhw == 1024 )
		VectorHashBody256_1024(data, h1, h2, h3, h4);
	else
	{
		cout << "Internal error: impossible value for rounded hash width: " << rhw << "." << endl;
		exit(1);
	}
}

void VectorHashBody512(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[], size_t hw)
{
	uint32_t rhw =  pow2roundup(hw);
	if( rhw == 32 )
		VectorHashBody512_32(data, h1, h2, h3, h4);
	else if( rhw == 64 )
		VectorHashBody512_64(data, h1, h2, h3, h4);
	else if( rhw == 128 )
		VectorHashBody512_128(data, h1, h2, h3, h4);
	else if( rhw == 256 )
		VectorHashBody512_256(data, h1, h2, h3, h4);
	else if( rhw == 512 )
		VectorHashBody512_512(data, h1, h2, h3, h4);
	else if( rhw == 1024 )
		VectorHashBody512_1024(data, h1, h2, h3, h4);
	else
	{
		cout << "Internal error: impossible value for rounded hash width: " << rhw << "." << endl;
		exit(1);
	}
}

static void VectorHash32(const void* buf, size_t len, uint32_t seed, void* out, size_t hw)
{
	uint32_t rhw =  pow2roundup(hw);
	if( rhw == 32 )
		VectorHash32_32(buf, len, seed, out, hw);
	else if( rhw == 64 )
		VectorHash32_64(buf, len, seed, out, hw);
	else if( rhw == 128 )
		VectorHash32_128(buf, len, seed, out, hw);
	else if( rhw == 256 )
		VectorHash32_256(buf, len, seed, out, hw);
	else if( rhw == 512 )
		VectorHash32_512(buf, len, seed, out, hw);
	else if( rhw == 1024 )
		VectorHash32_1024(buf, len, seed, out, hw);
	else
	{
		cout << "Internal error: impossible value for rounded hash width: " << rhw << "." << endl;
		exit(1);
	}
}

static void VectorHash128(const void* buf, size_t len, uint32_t seed, void* out, size_t hw)
{
	uint32_t rhw =  pow2roundup(hw);
	if( rhw == 32 )
		VectorHash128_32(buf, len, seed, out, hw);
	else if( rhw == 64 )
		VectorHash128_64(buf, len, seed, out, hw);
	else if( rhw == 128 )
		VectorHash128_128(buf, len, seed, out, hw);
	else if( rhw == 256 )
		VectorHash128_256(buf, len, seed, out, hw);
	else if( rhw == 512 )
		VectorHash128_512(buf, len, seed, out, hw);
	else if( rhw == 1024 )
		VectorHash128_1024(buf, len, seed, out, hw);
	else
	{
		cout << "Internal error: impossible value for rounded hash width: " << rhw << "." << endl;
		exit(1);
	}
}

static void VectorHash256(const void* buf, size_t len, uint32_t seed, void* out, size_t hw)
{
	uint32_t rhw =  pow2roundup(hw);
	if( rhw == 32 )
		VectorHash256_32(buf, len, seed, out, hw);
	else if( rhw == 64 )
		VectorHash256_64(buf, len, seed, out, hw);
	else if( rhw == 128 )
		VectorHash256_128(buf, len, seed, out, hw);
	else if( rhw == 256 )
		VectorHash256_256(buf, len, seed, out, hw);
	else if( rhw == 512 )
		VectorHash256_512(buf, len, seed, out, hw);
	else if( rhw == 1024 )
		VectorHash256_1024(buf, len, seed, out, hw);
	else
	{
		cout << "Internal error: impossible value for rounded hash width: " << rhw << "." << endl;
		exit(1);
	}
}

static void VectorHash512(const void* buf, size_t len, uint32_t seed, void* out, size_t hw)
{
	uint32_t rhw =  pow2roundup(hw);
	if( rhw == 32 )
		VectorHash512_32(buf, len, seed, out, hw);
	else if( rhw == 64 )
		VectorHash512_64(buf, len, seed, out, hw);
	else if( rhw == 128 )
		VectorHash512_128(buf, len, seed, out, hw);
	else if( rhw == 256 )
		VectorHash512_256(buf, len, seed, out, hw);
	else if( rhw == 512 )
		VectorHash512_512(buf, len, seed, out, hw);
	else if( rhw == 1024 )
		VectorHash512_1024(buf, len, seed, out, hw);
	else
	{
		cout << "Internal error: impossible value for rounded hash width: " << rhw << "." << endl;
		exit(1);
	}
}

API_EXPORT void VectorHashSIMD(const void* buf, size_t len, uint32_t seed, void* out, is_type SIMDversion, size_t hw)
{
	auto ibuf = reinterpret_cast<uintptr>(buf);
	// check if the alignment of the pointer in buf is OK
	// AVX512 requires 64-byte alignment
	// AVX2   requires 32-byte alignment
	// SSE2   requires 16-byte alignment
	// Scalar requires no special alignment
	if( SIMDversion >= IS_AVX512 && (ibuf&0x3f) == 0 )
	{
		VectorHash512(buf, len, seed, out, hw);
	}
	else if( SIMDversion >= IS_AVX2 && (ibuf&0x1f) == 0 )
	{
		VectorHash256(buf, len, seed, out, hw);
	}
	else if( SIMDversion >= IS_SSE2 && (ibuf&0x0f) == 0 )
	{
		VectorHash128(buf, len, seed, out, hw);
	}
	else if( SIMDversion >= IS_SCALAR )
	{
		VectorHash32(buf, len, seed, out, hw);
	}
	else
	{
		cout << "Internal error: impossible value for SIMD version: " << SIMDversion << "." << endl;
		exit(1);
	}
}

API_EXPORT void VectorHash(const void* buf, size_t len, uint32_t seed, void* out, size_t hw)
{
	VectorHashSIMD(buf, len, seed, out, GetSIMDVersion(), hw);
}

API_EXPORT int VectorHashStream(const vh_params& vhp, FILE* io, vector<uint32_t>& state)
{
	if( fseek( io, 0, SEEK_END ) != 0 )
		return -1;
	long fsize = ftell(io);
	if( fsize < 0 )
		return -2;
	state.resize(vhp.vh_nstate);
#if _POSIX_MAPPED_FILES > 0
	int fd = fileno(io);
	char* map = ( fsize > 0 ) ? (char*)mmap( NULL, fsize, PROT_READ, MAP_SHARED, fd, 0 ) : nullptr;
	if( fsize > 0 && map == MAP_FAILED )
		return -3;
	VectorHashSIMD( map, fsize, vhp.seed, state.data(), vhp.SIMDversion, vhp.vh_hash_width );
	munmap(map, fsize);
#else
	if( fseek( io, 0, SEEK_SET ) != 0 )
		return -1;
	void* map = NULL;
	if( fsize > 0 )
	{
		if( posix_memalign( &map, vh_hwreg_width/8, fsize ) != 0 )
			return -10;
		if( fread( map, fsize, 1, io ) != 1 )
			return -4;
	}
	VectorHashSIMD( map, fsize, vhp.seed, state.data(), vhp.SIMDversion, vhp.vh_hash_width );
	if( map != NULL )
		posix_memalign_free( map );
#endif

	return 0;
}

API_EXPORT int VectorHashStdin(const vh_params& vhp, vector<uint32_t>& state)
{
	void* h1;
	if( posix_memalign( &h1, vh_hwreg_width/8, vhp.vh_nint*sizeof(uint32_t) ) != 0 )
		return -10;
	void* h2;
	if( posix_memalign( &h2, vh_hwreg_width/8, vhp.vh_nint*sizeof(uint32_t) ) != 0 )
		return -10;
	void* h3;
	if( posix_memalign( &h3, vh_hwreg_width/8, vhp.vh_nint*sizeof(uint32_t) ) != 0 )
		return -10;
	void* h4;
	if( posix_memalign( &h4, vh_hwreg_width/8, vhp.vh_nint*sizeof(uint32_t) ) != 0 )
		return -10;

	uint32_t* z1 = (uint32_t*)h1;
	uint32_t* z2 = (uint32_t*)h2;
	uint32_t* z3 = (uint32_t*)h3;
	uint32_t* z4 = (uint32_t*)h4;

	uint32_t seed = vhp.seed;
	stateinit( z1, seed, vhp.vh_nint );
	stateinit( z2, seed, vhp.vh_nint );
	stateinit( z3, seed, vhp.vh_nint );
	stateinit( z4, seed, vhp.vh_nint );

	void* map = NULL;
	if( posix_memalign( &map, vh_hwreg_width/8, vhp.blocksize ) != 0 )
		return -10;
	state.resize(vhp.vh_nstate);
	size_t len = 0, bsize;
	bool lgContinue = true;
	while( lgContinue )
	{
		bsize = fread( map, 1, vhp.blocksize, stdin );
		if( bsize < vhp.blocksize ) {
			pad_buffer( (const uint8_t*)map, (uint8_t*)map, bsize, vhp.blocksize );
			lgContinue = false;
		}
		if( vhp.SIMDversion == IS_AVX512 )
			VectorHashBody512((const v16si*)map, (v16si*)h1, (v16si*)h2, (v16si*)h3, (v16si*)h4, vhp.vh_hash_width);
		else if( vhp.SIMDversion == IS_AVX2 )
			VectorHashBody256((const v8si*)map, (v8si*)h1, (v8si*)h2, (v8si*)h3, (v8si*)h4, vhp.vh_hash_width);
		else if( vhp.SIMDversion == IS_SSE2 )
			VectorHashBody128((const v4si*)map, (v4si*)h1, (v4si*)h2, (v4si*)h3, (v4si*)h4, vhp.vh_hash_width);
		else if( vhp.SIMDversion == IS_SCALAR )
			VectorHashBody32((const uint32_t*)map, (uint32_t*)h1, (uint32_t*)h2, (uint32_t*)h3, (uint32_t*)h4,
							 vhp.vh_hash_width);
		else {
			cerr << "Internal error: impossible value for SIMD version: " << vhp.SIMDversion << "." << endl;
			exit(1);
		}
		len += bsize;
	}
	posix_memalign_free( map );

	if( vhp.vh_nstate == 1 )
		VectorHashFinalize_32(len, z1, z2, z3, z4, state.data(), vhp.vh_hash_width);
	else if( vhp.vh_nstate == 2 )
		VectorHashFinalize_64(len, z1, z2, z3, z4, state.data(), vhp.vh_hash_width);
	else if( vhp.vh_nstate == 4 )
		VectorHashFinalize_128(len, z1, z2, z3, z4, state.data(), vhp.vh_hash_width);
	else if( vhp.vh_nstate == 8 )
		VectorHashFinalize_256(len, z1, z2, z3, z4, state.data(), vhp.vh_hash_width);
	else if( vhp.vh_nstate == 16 )
		VectorHashFinalize_512(len, z1, z2, z3, z4, state.data(), vhp.vh_hash_width);
	else if( vhp.vh_nstate == 32 )
		VectorHashFinalize_1024(len, z1, z2, z3, z4, state.data(), vhp.vh_hash_width);
	else {
		cerr << "Internal error: impossible value for vh_nstate: " << vhp.vh_nstate << "." << endl;
		exit(1);
	}

	posix_memalign_free( h4 );
	posix_memalign_free( h3 );
	posix_memalign_free( h2 );
	posix_memalign_free( h1 );

	return 0;
}
