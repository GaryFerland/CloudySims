//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#ifndef VECTORHASH_BINARY_H
#define VECTORHASH_BINARY_H

#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>

#include "vectorhash.h"
// #include "vectorhash_priv.h"

struct vh_params {
	string cmd;
	string name;
	bool lgBSDstyle;
	bool lgCheckMode;
	bool lgIgnoreMissing;
	bool lgBinarySet;
	bool lgTextSet;
	bool lgBinary;
	bool lgQuiet;
	bool lgStatusOnly;
	bool lgStrict;
	bool lgWarnSyntax;
	bool lgVerbose;
	bool lgZero;
	is_type SIMDversion;
	int returncode;
	uint32_t seed;
	size_t vh_hash_width;
	size_t vh_virtreg_width;
	size_t vh_nstate;
	size_t vh_nhash;
	size_t vh_nint;
	size_t blocksize;
	bool set_hash_width(size_t hw);
	char sentinel() const { return lgBinary ? '*' : ' '; }
	string option() const { return lgBinary ? "rb" : "r"; }
	void SetSIMDVersion();
	vh_params() : lgBSDstyle(false), lgCheckMode(false), lgIgnoreMissing(false), lgBinarySet(false),
				  lgTextSet(false), lgBinary(false), lgQuiet(false), lgStatusOnly(false), lgStrict(false),
				  lgWarnSyntax(false), lgVerbose(false), lgZero(false), SIMDversion(IS_INVALID),
				  returncode(0), seed(0xfd4c799d)
	{
		(void)set_hash_width(32);
	}
};

API_EXPORT int VectorHashStream(const vh_params& vhp, FILE* io, vector<uint32_t>& state);
API_EXPORT int VectorHashStdin(const vh_params& vhp, vector<uint32_t>& state);

#endif
