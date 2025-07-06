/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <stdint.h>
#include "vectorhash.h"

using namespace std;

int main(int argc, char** argv)
{
	if( argc < 2 )
	{
		cout << "usage: " << argv[0] << " <file>..." << endl;
		return 1;
	}
	for( int i=1; i < argc; ++i )
	{
		string vh128sum = VHstream(argv[i]);
		if( vh128sum == string() )
		{
			cerr << argv[0] << ": " << argv[i] << ": an error occurred while computing checksum" << endl;
			return 1;
		}
		cout << vh128sum << "  " << argv[i] << endl;
	}
	return 0;
}
