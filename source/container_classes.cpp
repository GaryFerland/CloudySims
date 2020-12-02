/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "container_classes.h"

#ifdef _MSC_VER
/* disable "'extern' before template explicit instantiation" */
#	pragma warning( disable : 4231 )
#endif

/* Explicit instantiations for debugging purposes */
INSTANTIATE_MULTI_ARR( bool, lgBOUNDSCHECKVAL )
INSTANTIATE_MULTI_ARR( char, lgBOUNDSCHECKVAL )
INSTANTIATE_MULTI_ARR( int, lgBOUNDSCHECKVAL )
INSTANTIATE_MULTI_ARR( long, lgBOUNDSCHECKVAL )
#ifndef FLT_IS_DBL
INSTANTIATE_MULTI_ARR( realnum,lgBOUNDSCHECKVAL )
#endif
INSTANTIATE_MULTI_ARR( double, lgBOUNDSCHECKVAL )

/** dump the array to a file in binary format; the file must
 *  already have been opened prior to calling this method */
void do_dump_state(const void* buf, size_t nelem, size_t size, FILE* out, int32 magic)
{
	DEBUG_ENTRY( "do_dump_state()" );

	bool lgErr = ( fwrite( &magic, sizeof(int32), 1, out ) != 1 );
	int32 help = (int32)sizeof(size_t);
	lgErr = lgErr || ( fwrite( &help, sizeof(int32), 1, out ) != 1 );
	lgErr = lgErr || ( fwrite( &size, sizeof(size_t), 1, out ) != 1 );
	lgErr = lgErr || ( fwrite( buf, size, nelem, out ) != nelem );
	if( lgErr )
	{
		fprintf( ioQQQ, " I/O error while dumping state!\n" );
		cdEXIT(EXIT_FAILURE);
	}
}

/** restore the array from a file in binary format; the file must already
 *  have been opened prior to calling this method and the array must already
 *  have been allocated in exactly the same way as when it was dumped; some
 *  checks are performed, but not every error is excluded */
void do_restore_state(void* buf, size_t nelem, size_t size, FILE *in, int32 magic)
{
	DEBUG_ENTRY( "do_restore_state()" );

	int32 help = 0;
	size_t help2 = 0;
	bool lgErr = ( fread( &help, sizeof(int32), 1, in ) != 1 );
	// this checks for correct version and prevents mixing up old style and new style data
	// it also prevents mixing up data from big-endian and little-endian machines.
	lgErr = lgErr || ( help != magic );
	lgErr = lgErr || ( fread( &help, sizeof(int32), 1, in ) != 1 );
	// this prevents mixing up data from 32-bit and 64-bit systems
	lgErr = lgErr || ( help != (int32)sizeof(size_t) );
	lgErr = lgErr || ( fread( &help2, sizeof(size_t), 1, in ) != 1 );
	// this may guard against reading an older, incompatible version of the array
	lgErr = lgErr || ( help2 != size );
	lgErr = lgErr || ( fread( buf, size, nelem, in ) != nelem );
	if( lgErr )
	{
		fprintf( ioQQQ, " Error while restoring state!\n" );
		cdEXIT(EXIT_FAILURE);
	}
}
