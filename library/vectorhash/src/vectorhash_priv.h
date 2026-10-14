//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_PRIV_H
#define VECTORHASH_PRIV_H

#include <string>
#include <type_traits>
#include <cstdlib>
#include <cstdint>

using namespace std;

static const string vh_version( "1.2" );

#define API_EXPORT __attribute__((visibility("default")))

#undef VH_INTEL
#ifdef MSVC

#if defined(_M_IX86) || defined(_M_X64)
#define VH_INTEL 1
#endif

#else

#if defined(i386) || defined(__i386__) || defined(__i386) || defined(__x86_64__)
#define VH_INTEL 1
#endif

#endif

#if defined(__CYGWIN__)
// _aligned_malloc / _aligned_free defined in Windows, but not Cygwin, reported by Richard Rudy
inline int posix_memalign(void **p, size_t a, size_t s)
{
	*p = aligned_alloc(s, a);
	return ( *p == NULL ) ? errno : 0;
}

inline void posix_memalign_free(void *p)
{
	free(p);
}
#elif defined(_MSC_VER)
#include <malloc.h>

// posix_memalign not defined on windows
inline int posix_memalign(void **p, size_t a, size_t s)
{
	*p = _aligned_malloc(s, a);
	return ( *p == NULL ) ? errno : 0;
}

inline void posix_memalign_free(void *p)
{
	_aligned_free(p);
}
#else
inline void posix_memalign_free(void *p)
{
	free(p);
}
#endif

//-----------------------------------------------------------------------------
// tunable parameters

// width of the hash (in bits)
// this value MUST be a power of 2 with a minimum of 32
#if defined(VH32)
static const size_t vh_hash_width = 32;
#define EXT(X) X##_32
#elif defined(VH64)
static const size_t vh_hash_width = 64;
#define EXT(X) X##_64
#elif defined(VH128)
static const size_t vh_hash_width = 128;
#define EXT(X) X##_128
#elif defined(VH256)
static const size_t vh_hash_width = 256;
#define EXT(X) X##_256
#elif defined(VH512)
static const size_t vh_hash_width = 512;
#define EXT(X) X##_512
#elif defined(VH1024)
static const size_t vh_hash_width = 1024;
#define EXT(X) X##_1024
#else
static const size_t vh_hash_width = 32;
#define EXT(X) X##_32
#endif

// width of largest hardware register that is supported (in bits)
static const size_t vh_hwreg_width = 512;

//-----------------------------------------------------------------------------
// derived quantities

// width of the virtual SIMD register supported in VectorHash (in bits)
// must be at least as wide as the largest hardware register that is used
static const size_t vh_virtreg_width = ( 2*vh_hash_width > vh_hwreg_width ) ? 2*vh_hash_width : vh_hwreg_width;
// width of the hash in uint32_t elements
static const size_t vh_nstate = vh_hash_width/32;
// number of uint32_t's that fit in the virtual register
static const size_t vh_nint = vh_virtreg_width/32;

// the file is read with this blocksize (in bytes)
static const size_t blocksize = 4*vh_nint*sizeof(uint32_t);

// unfortunately uintptr_t is optional, so we cannot rely on that...
typedef conditional<sizeof(int*) == 8, uint64_t, uint32_t>::type uintptr;

static_assert( sizeof(int*) == sizeof(uintptr), "definition of uintptr failed" );

void stateinit(uint32_t st[], uint32_t& seed);

//-----------------------------------------------------------------------------

#if defined(_MSC_VER)

#define ROTL32(x,y)	_rotl(x,y)

#else

inline uint32_t rotl32 ( uint32_t x, int r )
{
	// this expression is missing various safety checks for better speed
	return (x << r) | (x >> (32 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)

#endif

inline void pad_buffer(const uint8_t* src, uint8_t* buf, size_t len, size_t bufsz)
{
	for( size_t i=0; i < len; i++ )
		buf[i] = src[i];
	uint32_t c = 0;
	const uint32_t mask = 0xff;
	while( len < bufsz )
		buf[len++] = uint8_t(++c & mask);
}

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

inline uint32_t getblock32 ( uint32_t p )
{
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
	return (((p & 0xff000000) >> 24) | ((p & 0x00ff0000) >>  8) |
		((p & 0x0000ff00) <<  8) | ((p & 0x000000ff) << 24));
#else
	return p;
#endif
}

#define VEC(N, X) for( size_t j=0; j < (N); j++ ) { X; }

#endif
