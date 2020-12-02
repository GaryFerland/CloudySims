/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*LoadIsotopes	read in the nuclear isotope data and allocate space */
#include "cddefines.h"
#include "abund.h"
#include "atmdat.h"

/*LoadIsotopes	read in the nuclear isotope data and allocate space */
void LoadIsotopes ( )
{
	DEBUG_ENTRY( "SetIsotopeFractions()" );

	FILE *ioDATA;
	string chFile = "isotopes.dat";

	// first try local directory, then data/abundances
	if( (ioDATA = open_data( chFile.c_str(), "r", AS_LOCAL_ONLY_TRY ) ) == NULL )
	{
		char chPath[FILENAME_PATH_LENGTH_2] = { 0 };	/*entire path to file including file name*/
		strcat( chPath, chFile.c_str() );
		ioDATA = open_data( chPath, "r" );	// will abort if not found
	}

	char chLine[INPUT_LINE_LENGTH];			/*to store file lines*/
	bool lgEOL;

	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		if( chLine[0]=='*' )
			break;

		/* skip comment */
		if( chLine[0]=='#' )
			continue;
		if( chLine[0]=='\n' || chLine[0]=='\0' )
		{
			fprintf(ioQQQ,
				"PROBLEM in LoadIsotopes: Encountered unexpected empty line in %s.\n",
				chFile.c_str());
			cdEXIT(EXIT_FAILURE);
		}


		long i = 1;
		int ielem = (int)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL) - 1;
		ASSERT( ielem >= 0 );

		int Aiso  = (int) FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
		ASSERT( Aiso >  0  );

		FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);

		realnum mass = (realnum) FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
		ASSERT( mass >  0. );
	
	
		realnum spin = 0., magm = 0.;
		double tmp = FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
		if( tmp )
		{
			spin = (realnum) tmp;
			ASSERT( spin >= 0. );

			magm = (realnum) FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
			ASSERT( magm != 0. );
		}
		//	printf("ielem = %d\t Aiso = %d\t mass = %9.6f\t spin = %3.1f\t magm = %10.6f\n",
		//		ielem, Aiso, mass, spin, magm);

		abund.IsoAbn[ielem].setData( Aiso, mass, spin, magm);
	}

	fclose( ioDATA );
	return;
}
