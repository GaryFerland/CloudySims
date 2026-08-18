//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#include "TestMain.h"
#include "vectorhash_core.h"

namespace {

	uint32_t countbits(uint32_t n)
	{
		uint32_t res = 0;
		for( uint32_t k=0; k < 32; ++k )
		{
			if( (n&1) == 1 )
				++res;
			n >>= 1;
		}
		return res;
	}

	bool validate_result(uint32_t i, uint32_t n)
	{
		if( n < i )
			return false;
		if( n/2 >= i )
			return false;
		if( countbits(n) != 1 )
			return false;
		return true;
	}

	TEST(TestVHAdd)
	{
		for( uint32_t i=1; i <= 1024; ++i )
			for( uint32_t j=1; j <= 1024; ++j )
				CHECK( vh_add(i, j) == i+j );
		// check that overflow is handled correctly
		CHECK( vh_add(0xffffffffu, 1) == 0 );
		CHECK( vh_add(0xffffffffu, 2) == 1 );
		CHECK( vh_add(0xffffffffu, 3) == 2 );
	}

	TEST(TestPow2Roundup)
	{
		for( uint32_t i=1; i <= 1024; ++i )
			CHECK( validate_result(i, pow2roundup(i)) );
	}

}
