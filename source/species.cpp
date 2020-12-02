/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "species.h"
#include "taulines.h"
#include "input.h"
#include "dense.h"
#include "atmdat.h"
#include "elementnames.h"
#include "version.h"
#include "save.h"
#include "mole.h"
#include "service.h"
#include "parse_species.h"
#include "prt.h"

/*File nemala.cpp was developed by Humeshkar B Nemala as a part of his thesis work during the Summer of 2007*/
/* Initially the code has been developed to read in energy levels,radiative and
 * collisional data from the CHIANTI and LEIDEN databases. The idea is to extend it to more databases.
 * In the case of the Leiden database there is a single .dat file which has the energy levels information,
 * radiative and collisional data, with the data corresponding to each collider coming one after the other.
 * In the case of CHIANTI, the energy levels data, radiative data and collision data are present in seperate files.
 * While LEIDEN gives collisional rate coefficients, CHIANTI gives collisional strengths.
 * In the case of CHIANTI only two colliders are used:electrons and protons. They appear as separate files.
 * The electron collision strengths files are always expected to be there. A flag is set and data processed 
 * if the file on proton collision strengths is available.*/

/* There is an initialization file called species.ini which tells Cloudy what kind of data is to be used */
/* Structures are created separately to hold the transition data,radiative and collisional data */
/* The collisional structures are different for different databases depending upon whether */
/* collisional strengths or collisional rate coefficients are used.Finally a superstructure is constructed to hold */
/* the total collisional rate obtained by considering all the colliders */
/* The colliders considered are electron,proton,Atomic Hydrogen,He,He+,He++,Ortho Molecular Hydrogen,Para Molecular Hydrogen and Molecular Hydrogen */
STATIC void states_popfill(void);
STATIC void states_nelemfill(void);
STATIC void database_prep(int);
STATIC void trim_levels(long);
STATIC void set_fractionation( species *sp );
STATIC void states_propprint(void);

#define DEBUGSTATE false
void database_readin( void )
{
	int i,intNoSp;

	FILE *ioMASTERLIST, *ioVERSION;

	char *chToken;

	char chLine[FILENAME_PATH_LENGTH_2],
		chDLine[FILENAME_PATH_LENGTH_2],
		chPath[FILENAME_PATH_LENGTH_2] = "";

	const int MAX_NUM_SPECIES = 1000;

	char chLabels[MAX_NUM_SPECIES][CHARS_SPECIES];
	char chLabelsOrig[MAX_NUM_SPECIES][CHARS_SPECIES];
	char chPaths[MAX_NUM_SPECIES][FILENAME_PATH_LENGTH_2];

	static int nCalled = 0;
	long nSpeciesLAMDA, nSpeciesSTOUT, nSpeciesCHIANTI;

	DEBUG_ENTRY( "database_readin()" );

	/* only do this once. */
	if( nCalled > 0 )
	{
		return;
	}

	/* this is first call, increment the nCalled counterso never do this again */
	++nCalled;

	// read masterlists, count number of species
	nSpecies = 0;

	////////////////////////////////////
	//  
	// Read LAMDA masterlist 
	//
	//////////////////////////////////

	/* count how many lines are in the file, ignoring all lines
	 * starting with '#':This would give the number of molecules */
	nSpeciesLAMDA = 0;

	if( atmdat.lgLamdaOn )
	{
		long numModelsNotUsed = 0;
		strcpy( chPath, "lamda" );
		strcat( chPath, input.chDelimiter );
		strcat( chPath, "masterlist" );
		strcat( chPath, input.chDelimiter );
		strcat( chPath, atmdat.chLamdaFile );

		ioMASTERLIST = open_data( chPath, "r" );

		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of LAMDA masterlist.\n");
			cdEXIT(EXIT_FAILURE);
		}

		do
		{	
			if ((chLine[0]!='#') && (chLine[0]!='\n')&&(chLine[0]!='\t')&&(chLine[0]!='\r'))
			{	
				strcpy(chDLine, chLine);
				chToken = strtok(chDLine," \t\n");
				if( findspecies( chToken ) != null_mole  || 
					( chToken[1]=='-' && findspecies( chToken+2 ) != null_mole ) )
				{
					ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
					ASSERT( nSpeciesLAMDA + 1 <= MAX_NUM_SPECIES );
					ASSERT( strlen(chToken) < CHARS_SPECIES );
					strcpy( chLabels[nSpecies], chToken );
					chLabels[nSpecies][CHARS_SPECIES-1] = '\0';

					// path is, for example, LAMDA/no.dat
					strcpy( chPaths[nSpecies], "lamda" );
					strcat( chPaths[nSpecies], input.chDelimiter );
					chToken = strtok( NULL," \t\n" );
					strcat( chPaths[nSpecies], chToken );
					++nSpecies;
					++nSpeciesLAMDA;
				}
				else
					++numModelsNotUsed;
			}
		}
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) != NULL );

		/* \todo 1 - save this and stuff as note since not really a "PROBLEM" but worth reporting */
		//if( !t_version::Inst().lgRelease && numModelsNotUsed > 0 )
		//	fprintf( ioQQQ, "\n PROBLEM - %li LAMDA models could not be found in chemistry network.\n\n\n", numModelsNotUsed );

		fclose(ioMASTERLIST);
	}

	/* Print LAMDA molecule list if save data sources is on*/
	if( save.lgSDSOn && atmdat.lgLamdaOn)
	{
		fprintf(save.ipSDSFile, "##################################################\n");
		fprintf( save.ipSDSFile,"LAMDA (2005, A&A, 432, 369) molecules in this run.\n");
		for( int i=0; i<nSpeciesLAMDA; i++)
		{
			fprintf( save.ipSDSFile,"%s\t\t",chLabels[i]);
			if( (i+1)%5 == 0)
			{
				fprintf( save.ipSDSFile, "\n");
			}
		}
		fprintf(save.ipSDSFile,"\n\n");
	}

	//////////////////////////////////
	//  
	// Read CDMS/JPL masterlist
	//
	// These data files are in LAMDA format
	//
	//////////////////////////////////

	if( atmdat.lgCalpgmOn )
	{
		strcpy( chPath, "cdms+jpl" );
		strcat( chPath, input.chDelimiter );
		strcat( chPath, "masterlist" );

		ioMASTERLIST = open_data( chPath, "r" );

		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of CDMS/JPL masterlist.\n");
			cdEXIT(EXIT_FAILURE);
		}

		do
		{	
			if ((chLine[0]!='#') && (chLine[0]!='\n')&&(chLine[0]!='\t')&&(chLine[0]!='\r'))
			{	
				strcpy(chDLine, chLine);
				chToken = strtok(chDLine," \t\n");
				// hacks for alternative dialects...
				if( strcmp( chToken, "SH" ) == 0 )
					strcpy( chToken, "HS" );
				if( strcmp( chToken, "SH+" ) == 0 )
					strcpy( chToken, "HS+" );
				if( strcmp( chToken, "SD" ) == 0 )
					strcpy( chToken, "DS" );
				if( strcmp( chToken, "CCH" ) == 0 )
					strcpy( chToken, "C2H" );
				if( strcmp( chToken, "CCD" ) == 0 )
					strcpy( chToken, "C2D" );
				if( strcmp( chToken, "^17OO" ) == 0 )
					strcpy( chToken, "O^17O" );
				if( strcmp( chToken, "H^18O" ) == 0 )
					strcpy( chToken, "^18OH" );
				if( strcmp( chToken, "HCCD" ) == 0 )
					strcpy( chToken, "C2HD" );
				if( strcmp( chToken, "^13CCCH" ) == 0 )
					strcpy( chToken, "^13CC2H" );
				if( strcmp( chToken, "CC^13CH" ) == 0 )
					strcpy( chToken, "C2^13CH" );
				if( strcmp( chToken, "H^13CCCN" ) == 0 )
					strcpy( chToken, "H^13CC2N" );
				if( strcmp( chToken, "HCC^13CN" ) == 0 )
					strcpy( chToken, "HC2^13CN" );
				if( strcmp( chToken, "HCCC^15N" ) == 0 )
					strcpy( chToken, "HC3^15N" );
				// this molecule is cyclic, so these two are identical
				if( strcmp( chToken, "Si^13CC" ) == 0 )
					strcpy( chToken, "SiC^13C" );
				if( findspecies( chToken ) != null_mole  || 
				    ( chToken[1]=='-' && findspecies( chToken+2 ) != null_mole ) )
				{
					ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
					ASSERT( nSpeciesLAMDA + 1 <= MAX_NUM_SPECIES );
					strcpy( chLabels[nSpecies], chToken );
					chLabels[nSpecies][CHARS_SPECIES-1] = '\0';

					strcpy( chPaths[nSpecies], "cdms+jpl" );
					strcat( chPaths[nSpecies], input.chDelimiter );
					chToken = strtok( NULL," \t\n" );
					strcat( chPaths[nSpecies], chToken );
					++nSpecies;
					++nSpeciesLAMDA;
				}
				else
				{
					if( !t_version::Inst().lgRelease )
						fprintf( ioQQQ, "Warning: CDMS/JPL species %s not found\n", chToken );
				}
			}
		}
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) != NULL );

		fclose(ioMASTERLIST);
	}

	////////////////////////////////////
	//
	// Read STOUT masterlist and VERSION
	//
	///////////////////////////////////
	nSpeciesSTOUT = 0;

	//numLevels: index is nSpecies, value is the number of levels
	vector<long> numLevels(MAX_NUM_SPECIES,0L);

	if( atmdat.lgStoutOn )
	{
		// default location of Stout masterlist file
		strcpy( chPath, "stout" );
		strcat( chPath, input.chDelimiter );
		strcat( chPath, "masterlist" );
		strcat( chPath, input.chDelimiter );

		strcat( chPath, atmdat.chStoutFile );

		// first try local directory, then data/SED
		if( (ioMASTERLIST = open_data( atmdat.chStoutFile, "r", AS_LOCAL_ONLY_TRY ) ) == NULL )
		{
			ioMASTERLIST = open_data( chPath, "r" );
		}

		// magic number
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of stout.ini.\n");
			cdEXIT(EXIT_FAILURE);
		}

		bool lgEOL1=true, lgEOL2=true, lgEOL3=true;
		long int nMonRdST=-1, nDayRdST=-1;
		long int ipST = 1;
		long int nYrRdST = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOL1);
		if( !lgEOL1 )
		{
			nMonRdST = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOL2);
			if( !lgEOL2 )
				nDayRdST = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOL3);
		}
		if( lgEOL3 )
		{
			fprintf(ioQQQ,"PROBLEM, there must be three magic numbers on the first line of the stout masterlist file.\n");
			cdEXIT(EXIT_FAILURE);
		}

		static long int nYrST =11 , nMonST = 10, nDayST = 25;
		if( ( nYrRdST != nYrST ) || ( nMonRdST != nMonST ) || ( nDayRdST != nDayST ) )
		{
			fprintf( ioQQQ,
				" I expected to find the number %2.2li %2.2li %2.2li and got %2.2li %2.2li %2.2li instead.\n" ,
				nYrST , nMonST , nDayST , nYrRdST , nMonRdST , nDayRdST );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of CHIANTI masterlist.\n");
			cdEXIT(EXIT_FAILURE);
		}

		do
		{
			strcpy(chDLine, chLine);
			// three delimiters for tokesn, space, tab, newline
			// species name can have any number of columns
			// we will split line into two tokens, the name,
			// and the remainder of the lines
			chToken = strtok(chDLine," \t\n");

			if ((chLine[0]!='#') && (chLine[0]!='\n')&&(chLine[0]!='\t')&&(chLine[0]!='\r'))
			{
				ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
				ASSERT( nSpeciesSTOUT + 1 <= MAX_NUM_SPECIES );

				// first token is the species name
				strcpy( chLabels[nSpecies], chToken );
				strcpy( chLabelsOrig[nSpecies], chLabels[nSpecies] );

				// second optional token is lower limit to number of levels
				// first get full string after arbitrary length species name
				char *chNumLevs = strtok(NULL,"\n");
				if( chNumLevs != NULL )
				{
					long i = 1;
					bool lgEOL;
					long numLevs = (long)FFmtRead(chNumLevs,&i,sizeof(chLine),&lgEOL);
					// was there a lower bound to the number of levels
					if( !lgEOL )
					{
						if( numLevs > 0 )
						{
							numLevels[nSpecies] = numLevs;
						}
						else
						{
							fprintf(ioQQQ,"PROBLEM the limit to the number of levels must be positive, it was %li\n", numLevs);
							fprintf(ioQQQ,"The species was %s\n",chToken );
							cdEXIT( EXIT_FAILURE );
						}
					}
				}

				bool skipSpecies = false;

				//Check for duplicate species within Stout masterlist
				for( int j = nSpeciesLAMDA; j < nSpecies; j++)
				{
					if( strcmp( chLabelsOrig[j], chLabelsOrig[nSpecies] ) == 0)
					{
						fprintf(ioQQQ,"%s appears multiple times in %s.\n",chLabels[nSpecies],atmdat.chStoutFile);
						skipSpecies = true;
						break;
					}
				}

				if( skipSpecies )
					continue;

				char *chElement, chTokenTemp[7];
				strcpy( chTokenTemp, chToken );
				(void) strtok(chTokenTemp," \n");
				chElement = strtok(chTokenTemp,"_");
				uncaps( chElement );

				// path is, for example, CHIANTI/ar/ar_10/ar_10
				// we will append extensions later
				strcpy( chPaths[nSpecies], "stout" );
				strcat( chPaths[nSpecies], input.chDelimiter );
				strcat( chPaths[nSpecies], chElement );
				strcat( chPaths[nSpecies], input.chDelimiter );
				strcat( chPaths[nSpecies], chLabels[nSpecies] );
				strcat( chPaths[nSpecies], input.chDelimiter );
				strcat( chPaths[nSpecies], chLabels[nSpecies] );

				ASSERT( isalpha(chToken[0]) );
				long cursor=0;
				chLabels[nSpecies][0] = chToken[0];
				if( isalpha(chToken[1]) )
				{
					chLabels[nSpecies][1] = chToken[1];
					cursor = 2;
				}
				else
				{
					chLabels[nSpecies][1] = ' ';
					cursor = 1;
				}

				ASSERT( chToken[cursor]=='_' );
				++cursor;
				ASSERT( isdigit(chToken[cursor]) );

				if( isdigit(chToken[cursor+1]) )
				{
					chLabels[nSpecies][2] = chToken[cursor++];
					chLabels[nSpecies][3] =	chToken[cursor++];
				}
				else
				{
					chLabels[nSpecies][2] = ' ';
					chLabels[nSpecies][3] =	chToken[cursor++];
				}
				chLabels[nSpecies][4] = '\0';
				ASSERT( chToken[cursor]=='\0' || chToken[cursor]=='d' );

				// now capitalize the first letter
				chLabels[nSpecies][0] = toupper( chLabels[nSpecies][0] );
				++nSpecies;
				++nSpeciesSTOUT;
			}
		}
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) != NULL );
		fclose(ioMASTERLIST);
	}


	////////////////////////////////////
	//  
	// Read CHIANTI masterlist and VERSION
	//
	///////////////////////////////////

	nSpeciesCHIANTI = 0;

	if( atmdat.lgChiantiOn )
	{
		char chPathSave[FILENAME_PATH_LENGTH_2];
		strcpy( chPath, "chianti" );
		strcat( chPath, input.chDelimiter );
		//Preserve the path /chianti/ with chPathSave
		//Start reading in the chianti version number
		strcpy( chPathSave , chPath );
		strcat(chPath,"VERSION");
		ioVERSION = open_data(chPath,"r");
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioVERSION ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of the Chianti VERSION.\n");
			cdEXIT(EXIT_FAILURE);
		}
		fclose(ioVERSION);
		// chianti version - string since can contain letters
		for( i=0; i < atmdat.iVersionLength-1; i++ )
		{
			if( isprint(chLine[i]) )
				atmdat.chVersion[i] = chLine[i];
			else
			{
				atmdat.chVersion[i] = '\0';
				break;
			}
		}
		// make sure string is null-terminated
		atmdat.chVersion[atmdat.iVersionLength-1] = '\0';
		//Restore the previous chPath
		strcpy(chPath,chPathSave);
		// Read in the masterlist
		strcat( chPath, "masterlist" );
		strcat( chPath, input.chDelimiter );
		// save copy
		strcpy( chPathSave , chPath );

		// our subset of Chianti
		strcat( chPath, atmdat.chCloudyChiantiFile );

		// first try local directory, then data/chianti
		if( (ioMASTERLIST = open_data( atmdat.chCloudyChiantiFile, "r", AS_LOCAL_ONLY_TRY ) ) == NULL )
		{
			ioMASTERLIST = open_data( chPath, "r" );
		}

		// magic number
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of CloudyChianti.ini.\n");
			cdEXIT(EXIT_FAILURE);
		}

		bool lgEOL1=true, lgEOL2=true, lgEOL3=true;
		long int nMonRd=-1, nDayRd=-1;
		long int ipST = 1;
		long int nYrRd = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOL1);
		if( !lgEOL1 )
		{
			nMonRd = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOL2);
			if( !lgEOL2 )
				nDayRd = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOL3);
		}

		/* magic numbers for this version of Chianti masterlist */
		static long int nYr=11 , nMon = 10, nDay = 3;
		if( lgEOL3 )
		{
			fprintf(ioQQQ,"PROBLEM, there must be three magic numbers on the first line of the chianti masterlist file.\n");
			fprintf( ioQQQ,
				" I expected to find the numbers %2.2li %2.2li %2.2li.\n" ,
				nYr , nMon , nDay );
			cdEXIT(EXIT_FAILURE);
		}

		if( ( nYrRd != nYr ) || ( nMonRd != nMon ) || ( nDayRd != nDay ) )
		{
			fprintf( ioQQQ,
				" database_readin: the Chianti masterlist file is not the current version.\n" );
			fprintf( ioQQQ,
				" database_readin obtain the current version from the Cloudy web site.\n" );
			fprintf( ioQQQ,
				" I expected to find the number %2.2li %2.2li %2.2li and got %2.2li %2.2li %2.2li instead.\n" ,
				nYr , nMon , nDay , nYrRd , nMonRd , nDayRd );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}

		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of CHIANTI masterlist.\n");
			cdEXIT(EXIT_FAILURE);
		}

		do
		{

			if ((chLine[0]!='#') && (chLine[0]!='\n')&&(chLine[0]!='\t')&&(chLine[0]!='\r'))
			{
				// break line into two chunks, first with species which can have number number of
				// characters, followed by optional chunk with limit to number of levels
				strcpy(chDLine, chLine);
				chToken = strtok(chDLine," \t\n");

				fixit("insert logic here to exclude some ions");
				// (for example, iso sequences)
				// exclude for now the satellite lines (denoted by a "d" after the label
				if( chToken[strlen(chToken)-1] != 'd' )
				{
					ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
					ASSERT( nSpeciesCHIANTI + 1 <= MAX_NUM_SPECIES );
					strcpy( chLabels[nSpecies], chToken );
					strcpy( chLabelsOrig[nSpecies], chLabels[nSpecies]);

					// second optional token is lower limit to number of levels
					// first get full string after arbitrary length species name
					char *chNumLevs = strtok(NULL,"\n");
					if( chNumLevs != NULL )
					{
						long i = 1;
						bool lgEOL;
						long numLevs = (long)FFmtRead(chNumLevs,&i,sizeof(chLine),&lgEOL);
						// was there a lower bound to the number of levels
						if( !lgEOL )
						{
							if( numLevs > 0 )
							{
								numLevels[nSpecies] = numLevs;
							}
							else
							{
								fprintf(ioQQQ,"PROBLEM the limit to the number of levels must be positive, it was %li\n", numLevs);
								fprintf(ioQQQ,"The species was %s\n",chToken );
								cdEXIT( EXIT_FAILURE );
							}
						}
					}

					bool skipSpecies = false;

					//Check for duplicate species with Stout
					for( int j = nSpeciesLAMDA; j < (nSpecies - nSpeciesCHIANTI); j++)
					{
						if( strcmp( chLabelsOrig[j], chLabelsOrig[nSpecies] ) == 0)
						{
							fprintf(ioQQQ,"Skipping the Chianti version of %s, using Stout version\n",chLabels[nSpecies]);
							skipSpecies = true;
							break;
						}
					}
					//Check for duplicate species within Chianti masterlist
					for( int j = nSpecies - nSpeciesCHIANTI; j < nSpecies; j++)
					{
						if( strcmp( chLabelsOrig[j], chLabelsOrig[nSpecies] ) == 0)
						{
							fprintf(ioQQQ,"%s appears multiple times in %s.\n",chLabels[nSpecies],atmdat.chCloudyChiantiFile);
							skipSpecies = true;
							break;
						}
					}
					if( skipSpecies)
						continue;

					char *chElement, chTokenTemp[7];
					strcpy( chTokenTemp, chToken );
					(void) strtok(chTokenTemp," \n");
					chElement = strtok(chTokenTemp,"_");
					uncaps( chElement );

					// path is, for example, CHIANTI/ar/ar_10/ar_10
					// we will append extensions later
					strcpy( chPaths[nSpecies], "chianti" );
					strcat( chPaths[nSpecies], input.chDelimiter );
					strcat( chPaths[nSpecies], chElement );
					strcat( chPaths[nSpecies], input.chDelimiter );
					strcat( chPaths[nSpecies], chLabels[nSpecies] );
					strcat( chPaths[nSpecies], input.chDelimiter );
					strcat( chPaths[nSpecies], chLabels[nSpecies] );

					ASSERT( isalpha(chToken[0]) );
					long cursor=0;
					chLabels[nSpecies][0] = chToken[0];
					if( isalpha(chToken[1]) )
					{
						chLabels[nSpecies][1] = chToken[1];
						cursor = 2;
					}
					else
					{
						chLabels[nSpecies][1] = ' ';
						cursor = 1;
					}

					ASSERT( chToken[cursor]=='_' );
					++cursor;
					ASSERT( isdigit(chToken[cursor]) );

					if( isdigit(chToken[cursor+1]) )
					{
						chLabels[nSpecies][2] = chToken[cursor++];
						chLabels[nSpecies][3] =	chToken[cursor++];
					}
					else
					{
						chLabels[nSpecies][2] = ' ';
						chLabels[nSpecies][3] =	chToken[cursor++];
					}
					chLabels[nSpecies][4] = '\0';
					ASSERT( chToken[cursor]=='\0' || chToken[cursor]=='d' );

					// now capitalize the first letter
					chLabels[nSpecies][0] = toupper( chLabels[nSpecies][0] );
					++nSpecies;
					++nSpeciesCHIANTI;
				}
			}
		}
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) != NULL );

		fclose(ioMASTERLIST);
	}

	/* no species found, nothing to do */
	if( nSpecies==0 )
		return;

	/*Initialization of the dBaseSpecies Structure*/
	dBaseSpecies.resize(nSpecies);

	/*Initialization of the collisional rates array structure*/
	AtmolCollRateCoeff.reserve( nSpecies );
	AtmolCollSplines.resize(nSpecies);
	StoutCollData.resize(nSpecies);

	/*Mallocing here takes care of the number of colliders*/
	for( i=0; i<nSpecies; i++ )
	{
		AtmolCollRateCoeff.reserve( i, ipNCOLLIDER );
	}
	AtmolCollRateCoeff.alloc();

	// malloc state and transition arrays
	dBaseStates.resize(nSpecies);
	ipdBaseTrans.resize(nSpecies);

	for( i = 0; i < nSpecies; i++ )
	{
		dBaseTrans.push_back(TransitionList("dBaseTrans",&dBaseStates[i]));
		// label should be a minimum of 4 characters long
		size_t los = strlen(chLabels[i]);
		ASSERT( los >= 1 && los <= CHARS_SPECIES );
		dBaseSpecies[i].chLabel = new char[los+1];
		strcpy(dBaseSpecies[i].chLabel,chLabels[i]);
		dBaseSpecies[i].chLabel[los]='\0';
		trimTrailingWhiteSpace( dBaseSpecies[i].chLabel );
		dBaseSpecies[i].lgActive = true;

		// was minimum number of levels specified
		if( numLevels[i] > 0 )
		{
			dBaseSpecies[i].numLevels_masterlist = numLevels[i];
		}

		/* set type and isotopologue fractions */
		set_fractionation( &dBaseSpecies[i] );

		// set_fractionation trims off "p-","o-", etc.  Now have set label.  Check size.
		los = (int)strlen( dBaseSpecies[i].chLabel );
		ASSERT( los < CHARS_SPECIES );

		if( i<nSpeciesLAMDA )
		{
			// Read in LAMDA data files
			atmdat_LAMDA_readin( i, chPaths[i] );
		}
		else if( i < nSpeciesLAMDA + nSpeciesSTOUT )
		{
			// Read in STOUT data files
			atmdat_STOUT_readin( i, chPaths[i] );
		}
		else if( i < nSpeciesLAMDA + nSpeciesSTOUT + nSpeciesCHIANTI )
		{
			// Read in CHIANTI data files
			atmdat_CHIANTI_readin( i, chPaths[i] );
		}
		else
			TotalInsanity();
	}

	speciesCheck();

	states_popfill();
	states_nelemfill();

	for( long i=nSpeciesLAMDA; i<nSpeciesLAMDA+nSpeciesSTOUT; i++ )
	{
		strcpy(atmdat.chdBaseSources[dBaseStates[i][0].nelem()-1][dBaseStates[i][0].IonStg()-1],"Stout");
		atmdat.lgdBaseSourceExists[dBaseStates[i][0].nelem()-1][dBaseStates[i][0].IonStg()-1] = true;
	}
	for( long i=nSpeciesLAMDA+nSpeciesSTOUT; i<nSpeciesLAMDA+nSpeciesSTOUT+nSpeciesCHIANTI; i++ )
	{
		strcpy(atmdat.chdBaseSources[dBaseStates[i][0].nelem()-1][dBaseStates[i][0].IonStg()-1],"Chianti");
		atmdat.lgdBaseSourceExists[dBaseStates[i][0].nelem()-1][dBaseStates[i][0].IonStg()-1] = true;
	}

	if( save.lgSDSOn )
	{
		fprintf(save.ipSDSFile, "##################################################\n");
		fprintf( save.ipSDSFile,"Atomic model for each species used in this run.\n");
		fprintf( save.ipSDSFile,"Chianti (C), Stout(S), Iso-sequences (I), old internal treatment( ).\n\n");

		fprintf( save.ipSDSFile,"Ion");
		for(int i=0; i<LIMELM; i++)
		{
			fprintf( save.ipSDSFile,"%4d",i);
		}
		fprintf( save.ipSDSFile,"\n");

		for( int i=0; i<LIMELM; i++)
		{
			fprintf( save.ipSDSFile,"%s ",elementnames.chElementSym[i]);
			for(int j=0; j<i+1; j++)
			{
				fprintf( save.ipSDSFile,"   %c",atmdat.chdBaseSources[i][j][0]);
			}
			fprintf( save.ipSDSFile,"\n");
		}
	}

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\n\nDEBUG: Below are the contents of chdBaseSources[][]. It should contain a database name for each species.\n");
		fprintf( ioQQQ,"Ion");
		for(int i=0; i<LIMELM; i++)
		{
			fprintf( ioQQQ,"\t%i",i);
		}
		fprintf( ioQQQ,"\n");
		for( int i = 0; i < LIMELM; i++ )
		{
			fprintf( ioQQQ,"%s",elementnames.chElementSym[i]);
			for( int j = 0; j < LIMELM+1; j++ )
			{
				fprintf(ioQQQ,"\t%s",atmdat.chdBaseSources[i][j]);
			}
			fprintf(ioQQQ, "\n");
		}
		fprintf(ioQQQ,"\n\n");
	}

	/*Setting nelem of the states to an arbitrary value*/
	/*Also trim the highest levels if there are no valid Auls */
	/*Loop over species*/
	for( intNoSp=0; intNoSp<nSpecies; intNoSp++ )
	{
		database_prep(intNoSp);
		AllTransitions.push_back(dBaseTrans[intNoSp]);
		trim_levels(intNoSp);
	}

	/*To print the states*/
	if(DEBUGSTATE)
		states_propprint();
	return;
}

/** trim_levels: Trim the highest energy levels until there is a transition out of the highest
 * level with an Aul value that is not the default. The default value of 1e-30 is used when no Aul is
 * given in the database. */
STATIC void trim_levels(long ipSpecies)
{
	DEBUG_ENTRY( "trim_levels()" );

	bool lgLevelsToTrim = true;
	double aul;
	char* spectralLabel = dBaseSpecies[ipSpecies].chLabel;
	string speciesLabel = dBaseStates[ipSpecies].chLabel();

	long totalNumLevels = dBaseSpecies[ipSpecies].numLevels_max;

	while( lgLevelsToTrim )
	{
		long ipHi = dBaseSpecies[ipSpecies].numLevels_max-1;
		lgLevelsToTrim = true;

		if( dBaseSpecies[ipSpecies].numLevels_max == 0)
		{
			fprintf(ioQQQ,"PROBLEM: Spectrum %s (species: %s) has no transition probabilities out of the first %li levels.\n",
				spectralLabel, speciesLabel.c_str(), totalNumLevels);
			fprintf(ioQQQ,"Consider allowing Cloudy to use more levels (see Hazy 1 SPECIES STOUT/CHIANTI LEVELS MAX), add more low-level"
					" transition probabilities, or disable %s in the masterlist.\n\n", spectralLabel);
			cdEXIT(EXIT_FAILURE);
		}

		for( int ipLo = 0; ipLo < ipHi; ipLo++)
		{
			TransitionList::iterator tr = dBaseTrans[ipSpecies].begin()+ipdBaseTrans[ipSpecies][ipHi][ipLo];
			aul = tr->Emis().Aul();
			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"trim_levels():\t%s\t%i\t%li\t%e\n", spectralLabel, ipLo+1, ipHi+1, aul);
			}

			if( aul > atmdat.aulThreshold )
			{
				lgLevelsToTrim = false;
				break;
			}

		}
		if( lgLevelsToTrim )
		{
			--dBaseSpecies[ipSpecies].numLevels_max;
			dBaseSpecies[ipSpecies].numLevels_local = dBaseSpecies[ipSpecies].numLevels_max;
			if( atmdat.lgChiantiPrint || atmdat.lgStoutPrint || atmdat.lgLamdaPrint)
			{
				fprintf(ioQQQ,"Spectrum %s (species: %s) trimmed to %li levels (original %li) to have positive Aul.\n",
						spectralLabel,
						speciesLabel.c_str(),
						dBaseSpecies[ipSpecies].numLevels_local,
						totalNumLevels);
			}
		}
	}
}

STATIC void set_fractionation( species *sp )
{
	DEBUG_ENTRY( "set_fractionation()" );

	char chToken[3];

	sp->fracIsotopologue = 1.f;
	//types include "p-", "o-", "e-", and "a-"
	strncpy( chToken, sp->chLabel, 2 );
	chToken[2] = '\0';
	if( strcmp( "p-", chToken )==0 )
		sp->fracType = 0.25f;
	else if( strcmp( "o-", chToken )==0 )
		sp->fracType = 0.75f;
	else if( strcmp( "e-", chToken )==0 )
		sp->fracType = 0.5f;
	else if( strcmp( "a-", chToken )==0 ) 
		sp->fracType = 0.5f;
	else
		sp->fracType = 1.0f;

	fixit("what fraction should e-type and a-type Methanol have?  Assume 50/50 for now.");

	// Now scrape the type specifier off the label.
	if( sp->chLabel[1]=='-')
		memmove(sp->chLabel,sp->chLabel+2,strlen(sp->chLabel+2)+1);

	return;
}

/*This function zeros the population of all states */
STATIC void states_popfill( void)
{
	DEBUG_ENTRY( "states_popfill()" );

	for( long i=0; i<nSpecies; i++)
	{
		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++)
		{
			dBaseStates[i][j].Pop() = 0.;
		}
	}
	return;
}

void parsespect(char* chLabel, long& nelem, long& IonStg)
{
	DEBUG_ENTRY( "parsespect()" );	
	nelem = -1;
	IonStg = -1;
	if ( strlen(chLabel) != 4 || 
		  ! (isalpha(chLabel[0])) ||
		  ! (chLabel[1] == ' ' || isalpha(chLabel[1])) ||
		  ! (chLabel[2] == ' ' || isdigit(chLabel[2])) ||
		  ! (chLabel[3] == ' ' || isdigit(chLabel[3])))
	{
		// Invalid spectrum -- return error state
		return;
	}
	char chToken[3];
	strncpy( chToken, chLabel, 2 );
	chToken[2] = '\0';
	for( long ipElement=0; ipElement<LIMELM; ipElement++ )
	{
		if( strcmp( elementnames.chElementSym[ipElement], chToken )==0 )
		{
			nelem = ipElement;
			break;
		}
	}
	strncpy( chToken, chLabel + 2, 2 );
	IonStg = atoi(chToken);
}

string makeChemical( long nelem, long ion )
{
	DEBUG_ENTRY("makeChemical()");

	string chLabelChemical = elementnames.chElementSym[nelem];
	if( elementnames.chElementSym[nelem][1]==' ' )
		chLabelChemical = elementnames.chElementSym[nelem][0];

	// size has to be > 20 to prevent warnings about writing into an array that may be too small
	char chStage[21] = {'\0'};
	if( ion==1 )
		chStage[0] = '+';
	else if( ion>1 )
		sprintf( chStage, "+%li", ion );

	return chLabelChemical + chStage;
}

void makeChemical(char* chLabelChemical, long nelem, long ion)
{
	DEBUG_ENTRY("makeChemical()");
	string chemLab = makeChemical( nelem, ion );
	strncpy( chLabelChemical, chemLab.c_str(), size_t(CHARS_SPECIES) );
}

STATIC void spectral_to_chemical( char *chLabelChemical, char* chLabel, long &nelem, long &IonStg )
{
	DEBUG_ENTRY( "spectral_to_chemical()" );

	parsespect( chLabel, nelem, IonStg );
	ASSERT( nelem >= 0 && nelem < LIMELM );
	ASSERT( IonStg >= 1 && IonStg <= nelem+2 );

	//Prevent importing of iso-sequences from Chianti
	if( nelem - (IonStg-1) < NISO )
	{
		fprintf(ioQQQ, " PROBLEM: Cannot use Chianti model for %s%li\n",elementnames.chElementSym[nelem],IonStg);
		fprintf(ioQQQ, "  Iso-sequences are handled by our own model.\n");
		cdEXIT(EXIT_FAILURE);
	}

	makeChemical(chLabelChemical, nelem, IonStg-1);

	return;
}

void spectral_to_chemical( char *chLabelChemical, char* chLabel )
{
	DEBUG_ENTRY( "spectral_to_chemical()" );

	long nelem, IonStg;
	return spectral_to_chemical( chLabelChemical, chLabel, nelem, IonStg );
}

void chemical_to_spectral( const string chLabelChem, string &chLabelSpec )
{
	DEBUG_ENTRY( "chemical_to_spectral()" );

	size_t plus_sign_pos = chLabelChem.find_first_of( '+' );

	if( plus_sign_pos == string::npos )
	{
		/* Both 'H' and 'HCl' go through this branch */
		chLabelSpec = chLabelChem;

		if( chLabelSpec.length() == 1 )
			chLabelSpec += " ";

		if( chLabelSpec.length() == 2 &&
			isElementSym( chLabelSpec.c_str() ) )
		{
			chLabelSpec += " 1";
		}
	}
	else
	{
		/* Both 'C+2' and 'LiH+' go through this branch */
		string elm = chLabelChem.substr( 0, plus_sign_pos );

		if( elm.length() == 1 )
		{
			elm += " ";
		}

		if( ! isElementSym( elm.c_str() ) )
		{
			chLabelSpec = chLabelChem;
		}
		else
		{
			chLabelSpec = elm;
			int ionstg = atoi( chLabelChem.substr( plus_sign_pos+1 ).c_str() );
			if( ionstg == 0 )
				ionstg = 1;
			ionstg++;
			if( ionstg < 10 )
				chLabelSpec += " ";
			stringstream tmp;
			tmp << ionstg;
			chLabelSpec += tmp.str();
		}
	}
}
 
/*This function fills the nelem and IonStg fields */
STATIC void states_nelemfill(void)
{
	DEBUG_ENTRY( "states_nelemfill()" );

	for( long i=0; i<nSpecies; i++ )
	{
		long nelem = 0, IonStg;
		char chLabelChemical[CHARS_SPECIES] = "";

		if( dBaseSpecies[i].lgMolecular )
		{
			fixit("should never be used if lgMolecular"); 
			/* these should never be used if lgMolecular
			 *set to dangerous values instead of unity. */
			nelem = -1;
			IonStg = -1;
			strcpy( chLabelChemical, dBaseSpecies[i].chLabel );
		}
		else
		{
			spectral_to_chemical( chLabelChemical, dBaseSpecies[i].chLabel, nelem, IonStg );
			dBaseStates[i].chLabel_set( chLabelChemical );

			dBaseSpecies[i].fmolweight = dense.AtomicWeight[nelem];

			// do not evaluate our cooling if we are using Chianti for this species
			if( dBaseSpecies[i].database == "Chianti" )
			{
				dense.lgIonChiantiOn[nelem][IonStg-1] = true;
			}
			else if( dBaseSpecies[i].database == "Stout" )
			{
				dense.lgIonStoutOn[nelem][IonStg-1] = true;
			}
			else
			{
				TotalInsanity();
			}

			if( atmdat.lgChiantiHybrid || atmdat.lgStoutHybrid )
			{
				// used in cool_dima to indicate whether to include line
				// with shorter wl than these databases
				dense.maxWN[nelem][IonStg-1] = dBaseSpecies[i].maxWN;
			}
			else
			{
				dense.maxWN[nelem][IonStg-1] = 0.;
			}

			//Store the value of ipSpecies on C-scale
			//nelem(H) = 0  IonStg(H I) = 0
			atmdat.ipSpecIon[nelem][IonStg-1] = i;

			dBaseSpecies[i].lgPrtMatrix = false;
			if( strncmp( prt.matrix.species, dBaseStates[i].chLabel().c_str(), CHARS_SPECIES ) == 0 )
			{
				dBaseSpecies[i].lgPrtMatrix = true;
			}
		}

		molecule *sp = findspecies(chLabelChemical);
		if( sp == null_mole )
		{
			dBaseSpecies[i].index = INT_MAX;
			if( nelem >= ipHYDROGEN && dense.lgElmtOn[nelem] )
				fprintf(ioQQQ," PROBLEM: could not find species %li - %s\n",i,
					chLabelChemical );
		}
		else
		{
			dBaseSpecies[i].index = sp->index;
			mole.species[ sp->index ].dbase = &dBaseSpecies[i];
			mole.species[ sp->index ].levels = &dBaseStates[i];
			mole.species[ sp->index ].lines = &dBaseTrans[i];			
		}
		
		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++ )
		{
			dBaseStates[i][j].nelem() = nelem+1;
			dBaseStates[i][j].IonStg() = IonStg;
		}
	}
	return;
}

/*This function prints the various properties of states*/
STATIC void states_propprint(void)
{
	DEBUG_ENTRY( "states_propprint()" );

	for( long i=0; i<nSpecies; i++ )
	{
		printf("The species is %s \n",dBaseSpecies[i].chLabel);
		printf("The data output is in the following format \n");
		printf("Label Energy St.wt Pop Lifetime\n");

		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++ )
		{
			printf("This is the %ld state \n",j);
			printf("%s %f %f %f %e \n",dBaseStates[i][j].chLabel().c_str(),
				dBaseStates[i][j].energy().WN(),
				dBaseStates[i][j].g(),
				dBaseStates[i][j].Pop(),
				dBaseStates[i][j].lifetime());
		}
	}
	return;
}

STATIC void database_prep(int intSpIndex)
{
	vector<realnum> fsumAs(dBaseSpecies[intSpIndex].numLevels_max,SMALLFLOAT);

	DEBUG_ENTRY( "database_prep()" );

	/*Get the lifetimes*/
	for( EmissionList::iterator em = dBaseTrans[intSpIndex].Emis().begin();
		  em != dBaseTrans[intSpIndex].Emis().end(); ++em) 
	{
		fsumAs[(*em).Tran().ipHi()] += (*em).Aul();

		// set redistribution functions for all lines
		if( intSpIndex != atmdat.ipSpecIon[ipIRON][1] )
		{
			// default for species is partial redisctribution with wings
			(*em).iRedisFun() = ipPRD;
		}
		else
		{
			// this Fe II is to be the default for cloudy post 2001
			if( em->Tran().ipLo() == 0 )
			{
				// complete redistribution, only core
				em->iRedisFun() = ipCRD;
			}
			else
			{
				/* >>chng 01 feb 27, had been -1, crd with core only,
				 * change to complete redistribution with wings as per discussion with Ivan Hubeny */
				em->iRedisFun() = ipCRDW;
			}
		}
	}
	
	dBaseStates[intSpIndex][0].lifetime()= BIGFLOAT;
	for( int ipHi=1; ipHi < dBaseSpecies[intSpIndex].numLevels_max; ipHi++ )
	{
		dBaseStates[intSpIndex][ipHi].lifetime() = 1./fsumAs[ipHi];
	}
	return;
}
