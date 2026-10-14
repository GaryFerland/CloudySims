//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include "vectorhash_priv.h"
#include "vectorhash_core.h"
#include "vectorhash_binary.h"

static const string SIMDname[] = { "Scalar", "SSE2", "AVX2", "AVX512" };

API_EXPORT bool vh_params::set_hash_width(size_t hw)
{
	// width of the hash (in bits)
	// this value MUST be a multiple of 32 between 32 and 1024
	if( hw < 32 || hw > 1024 || (hw & size_t{0x1f}) != 0 )
		return false;
	// initially we will compute a hash that is a power of two wide
	// as a last step that width will be reduced to the requested width
	vh_hash_width = pow2roundup(hw);
	// width of the virtual SIMD register supported in VectorHash (in bits)
	// must be at least as wide as the largest hardware register that is used
	vh_virtreg_width = ( 2*vh_hash_width > vh_hwreg_width ) ? 2*vh_hash_width : vh_hwreg_width;
	// width of the rounded hash in uint32_t elements
	vh_nstate = vh_hash_width/32;
	// number of uint32_t's that fit in the virtual register
	vh_nint = vh_virtreg_width/32;
	// the file is read with this blocksize (in bytes)
	blocksize = 4*vh_nint*sizeof(uint32_t);
	// the actual requested hash width
	vh_hash_width = hw;
	// width of the actual hash in uint32_t elements
	vh_nhash = hw/32;
	// update the name as well...
	ostringstream oss;
	oss << "VH" << vh_hash_width;
	name = oss.str();
	return true;
}

API_EXPORT void vh_params::SetSIMDVersion()
{
	// check if SIMD version was already forced by a command line option
	if( SIMDversion == IS_INVALID )
	{
		SIMDversion = GetSIMDVersion();
		if( lgVerbose )
		{
			cout << "found SIMD capability: " << SIMDname[SIMDversion] << "." << endl;
		}
	}
	else
	{
		if( lgVerbose )
		{
			cout << "SIMD capability was set on the command line to: " << SIMDname[SIMDversion] << "." << endl;
		}
	}			
}
