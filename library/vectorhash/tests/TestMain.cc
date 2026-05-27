//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#include <cstdlib>
#include <UnitTest++/TestReporterStdout.h>
#include <cpuinfo.hpp>
#include "TestMain.h"

void* buffer;
is_type SIMDversion;
uint32_t cksum[1024/32];

int main ()
{
	void *buffer_raw = malloc((1<<20) + 63);

	// make sure buffer is aligned on a 64-byte boundary
	uintptr ibuf = reinterpret_cast<uintptr>(buffer_raw);
	uintptr mask = 0x3f;
	ibuf = (ibuf+63)&(~mask);
	buffer = reinterpret_cast<void*>(ibuf);

	// determine SIMD capability, this determines which tests can be run
	SIMDversion = GetSIMDVersion();

	int retval = UnitTest::RunAllTests();

	free(buffer_raw);
	return retval;
}
