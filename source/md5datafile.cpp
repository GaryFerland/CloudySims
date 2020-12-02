/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "thirdparty.h"

int main(int argc, char** argv)
{
	if( argc != 2 )
	{
		printf( "usage: %s <filename>\n", argv[0] );
		return 1;
	}
	string md5sum = MD5datafile( argv[1] );
	printf( "%s  %s\n", md5sum.c_str(), argv[1] );
	return 0;
}
