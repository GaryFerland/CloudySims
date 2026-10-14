//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include "vectorhash_priv.h"
#include "vectorhash_core.h"

inline size_t operator ""_z (unsigned long long n)
{
	return n;
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

void EXT(VectorHashFinalize)(size_t len, uint32_t* h1, uint32_t* h2, uint32_t* h3, uint32_t* h4, void* out, size_t hw)
{
	uint32_t lenc = fmix32(len & 0x7fffffff);

	VEC( vh_nint, h1[j] ^= lenc );
	VEC( vh_nint, h2[j] ^= lenc );
	VEC( vh_nint, h3[j] ^= lenc );
	VEC( vh_nint, h4[j] ^= lenc );

	for( size_t j=1; j < vh_nint; j++ )
	{
		h1[0] = vh_add(h1[0], h2[j]);
		h2[0] = vh_add(h2[0], h3[j]);
		h3[0] = vh_add(h3[0], h4[j]);
		h4[0] = vh_add(h4[0], h1[j]);
		h1[j] = vh_add(h1[j], h2[0]);
		h2[j] = vh_add(h2[j], h3[0]);
		h3[j] = vh_add(h3[j], h4[0]);
		h4[j] = vh_add(h4[j], h1[0]);
	}

	VEC( vh_nint, h1[j] = vh_add(h1[j], h2[j]) );
	VEC( vh_nint, h1[j] = vh_add(h1[j], h3[j]) );
	VEC( vh_nint, h1[j] = vh_add(h1[j], h4[j]) );
	VEC( vh_nint, h2[j] = vh_add(h2[j], h1[j]) );
	VEC( vh_nint, h3[j] = vh_add(h3[j], h1[j]) );
	VEC( vh_nint, h4[j] = vh_add(h4[j], h1[j]) );
	VEC( vh_nint, h1[j] = fmix32(h1[j], h3[j]) );
	VEC( vh_nint, h2[j] = fmix32(h2[j], h4[j]) );
	VEC( vh_nint, h3[j] = fmix32(h3[j], h1[j]) );
	VEC( vh_nint, h4[j] = fmix32(h4[j], h2[j]) );
	VEC( vh_nint, h1[j] = vh_add(h1[j], h2[j]) );
	VEC( vh_nint, h1[j] = vh_add(h1[j], h3[j]) );
	VEC( vh_nint, h1[j] = vh_add(h1[j], h4[j]) );
	VEC( vh_nint, h2[j] = vh_add(h2[j], h1[j]) );
	VEC( vh_nint, h3[j] = vh_add(h3[j], h1[j]) );
	VEC( vh_nint, h4[j] = vh_add(h4[j], h1[j]) );

	for( size_t j=1; j < vh_nint; j++ )
	{
		h1[0] = vh_add(h1[0], h4[j]);
		h2[0] = vh_add(h2[0], h1[j]);
		h3[0] = vh_add(h3[0], h2[j]);
		h4[0] = vh_add(h4[0], h3[j]);
		h1[j] = vh_add(h1[j], h4[0]);
		h2[j] = vh_add(h2[j], h1[0]);
		h3[j] = vh_add(h3[j], h2[0]);
		h4[j] = vh_add(h4[j], h3[0]);
	}

	size_t lc_nstate = max(vh_nstate, 4_z);
	vector<uint32_t> lres(lc_nstate);
	vector<uint32_t*> z(lc_nstate);
	size_t j = 0;
	size_t nn = vh_nint*4/lc_nstate;
	for( size_t i=0; i < lc_nstate/4; i++ )
		z[j++] = &h1[i*nn];
	for( size_t i=0; i < lc_nstate/4; i++ )
		z[j++] = &h2[i*nn];
	for( size_t i=0; i < lc_nstate/4; i++ )
		z[j++] = &h3[i*nn];
	for( size_t i=0; i < lc_nstate/4; i++ )
		z[j++] = &h4[i*nn];

	// fold the result to the desired hash width
	size_t ns2 = lc_nstate/2;
	for( size_t i=0; i < lc_nstate; i++ )
	{
		lres[i] = fmix32(z[i][0] ^ z[(i+ns2)%lc_nstate][nn-1]);
		for( size_t j=1; j < nn; j++ )
			lres[i] ^= fmix32(z[(j+i)%lc_nstate][j] ^ z[(j+i+ns2)%lc_nstate][nn-1-j]);
	}

	while( lc_nstate > vh_nstate )
	{
		for( size_t i=0; i < lc_nstate/2; ++i )
			lres[i] = fmix32(lres[i], lres[lc_nstate-1-i]);
		lc_nstate /= 2;
	}

	size_t nw = hw/32;
	if( nw < vh_nstate )
	{
		if( nw < vh_nstate/2 )
		{
			cerr << "Internal error: impossible value for hash width\n";
			exit(1);
		}
		for( size_t i=0; i < nw; ++i )
			lres[i] = fmix32(lres[i], lres[i+vh_nstate-nw]);
	}

	uint32_t* res = (uint32_t*)out;
	for( size_t i=0; i < nw; i++ )
		res[i] = lres[i];
}

