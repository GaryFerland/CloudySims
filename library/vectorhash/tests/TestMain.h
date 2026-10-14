//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef TESTMAIN_H
#define TESTMAIN_H

#include <cstdint>
#include <cctype>
#include <string>
#include <fstream>
#include <UnitTest++/UnitTest++.h>
#include "vectorhash.h"
#include "vectorhash_core.h"

using namespace std;

extern void* buffer;
extern is_type SIMDversion;
extern uint32_t cksum[1024/32];

inline uint32_t hexval(unsigned char c)
{
	if( isdigit(c) )
		return uint32_t(c - '0');
	else
		return uint32_t(c - 'a' + 10);
}

inline uint32_t ReadChunk(uint32_t pos, const string& hashstr)
{
	if( pos*8 >= hashstr.size() )
		return 0;
	uint32_t res = 0;
	for( uint32_t i=0; i < 8; ++i )
	{
		res = res << 4;
		res += hexval(hashstr[pos*8+i]);
	}
	return res;
}

inline bool CheckHash(const uint32_t hash[], const string& hashstr)
{
	if( hashstr.size() < 8 || hashstr.size() > 256 || (hashstr.size()&0x7) != 0 )
		return false;
	uint32_t numchunk = hashstr.size()/8;
	for( uint32_t i=0; i < numchunk; ++i )
		if( ReadChunk(i, hashstr) != hash[i] )
			return false;
	return true;
}

inline bool ReadBuffer(const string& fname, uint32_t size, void* buf)
{
	ifstream io(fname, ifstream::binary);
	if( !io.is_open() )
		return false;
	io.read((char*)buf, size);
	return io.good();
}

#endif
