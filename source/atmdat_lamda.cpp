/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "atmdat.h"
#include "rfield.h"
#include "taulines.h"
#include "trace.h"
#include "lines_service.h"
#include "parse_species.h"

static const bool DEBUGSTATE=false;

static void check_LAMDA_comment(const char* chLine, char* chEFilename)
{
	const char* s= chLine;
	while (isspace(*s))
		++s;
	if (*s != '!')
	{
		fprintf (ioQQQ," PROBLEM: The LAMDA data file '%s' appears to be corrupt\n", 
					chEFilename);
		fprintf (ioQQQ," Expected line starting '!', found\n%s\n",chLine);
		cdEXIT(EXIT_FAILURE);
	}
}

/*Separate Routine to read in the molecules*/
void atmdat_LAMDA_readin( long intNS, char *chEFilename )
{
	DEBUG_ENTRY( "atmdat_LAMDA_readin()" );

	long int intlnct;

	FILE *ioLevData;
	realnum  fstatwt,fenergyWN,fenergy,feinsteina;

	char chLine[FILENAME_PATH_LENGTH_2];

	ASSERT( intNS >= 0 );

	static const int MAX_NUM_LEVELS = 999;

	dBaseSpecies[intNS].lgMolecular = true;
	dBaseSpecies[intNS].lgLTE = false;

	/*Load the species name in the dBaseSpecies array structure*/
	if(DEBUGSTATE)
		printf("The name of the %li species is %s \n",intNS+1,dBaseSpecies[intNS].chLabel);

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," moldat_readin opening %s:",chEFilename);

	ioLevData = open_data( chEFilename, "r" );

	intlnct = 0;
	while(intlnct < 3)
	{
		intlnct++;
		if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
	}
	/*Extracting out the molecular weight*/
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	dBaseSpecies[intNS].fmolweight = (realnum)atof(chLine);

	/*Discard this line*/
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	// sections of the file are separated by line that begin with "!"
	check_LAMDA_comment(chLine,chEFilename);

	/*Reading in the number of energy levels*/
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	long nMolLevs = (long)atoi(chLine);

	long HighestIndexInFile = nMolLevs;
	
	dBaseSpecies[intNS].numLevels_max = HighestIndexInFile;

	setProperties(dBaseSpecies[intNS]);

	/* truncate model to preset max */
	nMolLevs = MIN3( atmdat.nLamdaMaxLevels, nMolLevs, MAX_NUM_LEVELS );

	if(nMolLevs <= 0)
	{
		fprintf( ioQQQ, "The number of energy levels is non-positive in datafile %s.\n", chEFilename );
		fprintf( ioQQQ, "The file must be corrupted.\n" );
		cdEXIT( EXIT_FAILURE );
	}

	if( atmdat.lgLamdaPrint == true)
	{
		fprintf( ioQQQ,"Using LAMDA model %s with %li levels of %li available.\n",
				dBaseSpecies[intNS].chLabel , nMolLevs , HighestIndexInFile );
	}

	if (dBaseSpecies[intNS].setLevels != -1)
	{
		if (dBaseSpecies[intNS].setLevels > HighestIndexInFile)
		{
			fprintf( ioQQQ,"Using LAMDA model %s with %li requested, only %li energy levels available.\n",
				dBaseSpecies[intNS].chLabel ,dBaseSpecies[intNS].setLevels, HighestIndexInFile );
			nMolLevs = HighestIndexInFile;		  
		}
		else
		{
			nMolLevs = dBaseSpecies[intNS].setLevels;
		}
	}

	dBaseSpecies[intNS].numLevels_max = nMolLevs;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	/*malloc the States array*/
	dBaseStates[intNS].init(dBaseSpecies[intNS].chLabel,nMolLevs);
	/* allocate the Transition array*/
	ipdBaseTrans[intNS].reserve(nMolLevs);
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
		ipdBaseTrans[intNS].reserve(ipHi,ipHi);
	ipdBaseTrans[intNS].alloc();
	dBaseTrans[intNS].resize(ipdBaseTrans[intNS].size());
	dBaseTrans[intNS].states() = &dBaseStates[intNS];
	dBaseTrans[intNS].chLabel() = dBaseSpecies[intNS].chLabel;
	dBaseSpecies[intNS].database = "LAMDA";

	int ndBase = 0;
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
	{
		for( long ipLo = 0; ipLo < ipHi; ipLo++)
		{
			ipdBaseTrans[intNS][ipHi][ipLo] = ndBase;
			dBaseTrans[intNS][ndBase].Junk();
			dBaseTrans[intNS][ndBase].setLo(ipLo);
			dBaseTrans[intNS][ndBase].setHi(ipHi);
			++ndBase;
		}
	}

	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	// sections of the file are separated by line that begin with "!"
	check_LAMDA_comment(chLine,chEFilename);

	for( long ipLev=0; ipLev<HighestIndexInFile; ipLev++)
	{	
		if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		// skip these levels
		if( ipLev >= nMolLevs )
			continue;

		long i = 1;
		bool lgEOL;
		long index = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		fenergy = (realnum)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		fstatwt = (realnum)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );

		if( index != ipLev + 1 )
		{
			fprintf( ioQQQ, "internal error in file %s\n", chEFilename );
			cdEXIT(EXIT_FAILURE);
		}
		dBaseStates[intNS][ipLev].energy().set(fenergy,"cm^-1");
		dBaseStates[intNS][ipLev].g() = fstatwt;

		if( ipLev > 0 )
		{
			// Use volatile variables to ensure normalization to
			// standard precision before comparison
			volatile realnum enlev = dBaseStates[intNS][ipLev].energy().Ryd();
			volatile realnum enprev = dBaseStates[intNS][ipLev-1].energy().Ryd();
			if (enlev < enprev)
			{
				fprintf( ioQQQ, " The energy levels are not in order in species %s at index %li.\n",
					dBaseSpecies[intNS].chLabel, ipLev );
				cdEXIT(EXIT_FAILURE);
			}
		}
		if(DEBUGSTATE)
		{
			printf("The converted energy is %f \n",dBaseStates[intNS][ipLev].energy().WN());
			printf("The value of g is %f \n",dBaseStates[intNS][ipLev].g());
		}
	}

	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for(TransitionList::iterator tr=dBaseTrans[intNS].begin();
		 tr!= dBaseTrans[intNS].end(); ++tr)
	{
		int ipHi = (*tr).ipHi();
		int ipLo = (*tr).ipLo();
		//fenergyWN = (realnum)MAX2( 1.01*rfield.emm()*RYD_INF, dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN() );
		fenergyWN = (realnum)(dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN());
		
		(*tr).EnergyWN() = fenergyWN;

		/* there are OH hyperfine levels where i+1 and i have exactly
		 * the same energy.  The refractive index routine will FPE with
		 * an energy of zero - so we do this test */
		if( fenergyWN>SMALLFLOAT )
			(*tr).WLAng() = (realnum)(1e+8f/fenergyWN/RefIndex(fenergyWN));
		else
			(*tr).WLAng() = 1e30f;
	}

	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	if(chLine[0]!='!')
	{
		fprintf( ioQQQ, " The number of energy levels in file %s is not correct, expected to find line starting with!.\n",chEFilename);
		fprintf( ioQQQ , "%s\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	long intgrtct = atoi(chLine);
	/*The above gives the number of radiative transitions*/
	if(intgrtct <= 0)
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	if(DEBUGSTATE)
	{
		printf("The number of radiative transitions is %li \n",intgrtct);
	}
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	/*intrtct refers to radiative transitions count*/
	for( long intrtct=0; intrtct<intgrtct; intrtct++)
	{
		if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		long i = 1;
		bool lgEOL;

		long index = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		long ipHiInFile = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		long ipLoInFile = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		long ipHi = ipHiInFile - 1;
		long ipLo = ipLoInFile - 1;

		if( ! ( ipHi >= 0 && ipLo >= 0 && ipHi > ipLo && index == intrtct + 1 ) )
		{
			fprintf( ioQQQ, "internal error in file %s\n", chEFilename );
			cdEXIT(EXIT_FAILURE);
		}

		// skip these lines
		if( ipLo >= nMolLevs || ipHi >= nMolLevs )
			continue;

		feinsteina = (realnum)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		/* don't need the energy in GHz, so throw it away. */
		FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		fenergyWN = (realnum)(dBaseStates[intNS][ipHi].energy().WN() -dBaseStates[intNS][ipLo].energy().WN());
		fenergyWN = MAX2( fenergyWN, 1.01 * RYD_INF * rfield.emm() );

		TransitionList::iterator tr = dBaseTrans[intNS].begin()+ipdBaseTrans[intNS][ipHi][ipLo];
		ASSERT(!(*tr).hasEmis()); // Check for duplicates
		(*tr).AddLine2Stack();
		(*tr).Emis().Aul() = feinsteina;
		ASSERT( !isnan( (*tr).EnergyK() ) );
		(*tr).EnergyWN() = fenergyWN;
		(*tr).WLAng() = (realnum)(1e+8/fenergyWN/RefIndex(fenergyWN));

		(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(),(*tr).EnergyWN(), (*(*tr).Hi()).g());

		(*tr).setComment( db_comment_tran_levels( ipLoInFile, ipHiInFile ) );

		if(DEBUGSTATE)
		{
			printf("The upper level is %ld \n",ipHi+1);
			printf("The lower level is %ld \n",ipLo+1);
			printf("The Einstein A  is %E \n",(*tr).Emis().Aul());
			printf("The Energy of the transition is %E \n",(*tr).EnergyK());
		}
	}

	if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	if(chLine[0]!='!')
	{
		fprintf( ioQQQ, " The number of radiative transitions in file %s is not correct.\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	long nCollPartners = -1;
	/*Getting the number of collisional partners*/
	if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		nCollPartners = atoi(chLine);
	}
	/* if the number of colliders has a negative value, do the molecule in LTE */
	if (nCollPartners < 0)
		dBaseSpecies[intNS].lgLTE = true;
	/*Checking the number of collisional partners does not exceed 9*/
	if( nCollPartners > ipNCOLLIDER )
	{
		fprintf( ioQQQ, " The number of colliders is greater than what is expected in file %s.\n", chEFilename );
		cdEXIT(EXIT_FAILURE);
	}

	FunctPtr func = new FunctLAMDA();
	// loop over partners	
	for( long ipPartner = 0; ipPartner < nCollPartners; ++ ipPartner )
	{
		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		check_LAMDA_comment(chLine,chEFilename);

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
		/*Extract out the name of the collider*/
		/*The following are the rules expected in the datafiles to extract the names of the collider*/
		/*The line which displays the species and the collider starts with a number*/
		/*This refers to the collider in the Leiden database*/
		/*In the Leiden database 1 refers to H2,2 to para-H2,3 to ortho-H2
		4 to electrons,5 to H and 6 to He*/
		char *chCollName = strtok(chLine," ");
		/*Leiden Collider Index*/
		int intLColliderIndex = atoi(chCollName);
		int intCollIndex = -1;
		/*In Cloudy,We assign the following indices for the colliders*/
		/*electron=0,proton=1,He+=2,He++=3,atomic hydrogen=4,He=5,oH2=6,pH2=7,H2=8*/

		if(intLColliderIndex == 1)
		{
			intCollIndex = ipH2;
		}
		else if(intLColliderIndex == 2)
		{
			intCollIndex = ipH2_PARA;
		}
		else if(intLColliderIndex == 3)
		{
			intCollIndex = ipH2_ORTHO;
		}
		else if(intLColliderIndex == 4)
		{
			intCollIndex = ipELECTRON;
		}
		else if(intLColliderIndex == 5)
		{
			intCollIndex = ipATOM_H;
		}
		else if(intLColliderIndex == 6)
		{
			intCollIndex = ipATOM_HE;
		}
		else
		{
			// this happens for some LAMDA files (as of Jan 20, 2009) because there is no integer
			// in the datafile indicating which collider the subsequent table is for
			TotalInsanity();
		}

		ASSERT( intCollIndex < ipNCOLLIDER );

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		check_LAMDA_comment(chLine,chEFilename);

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
		// Get the number of transitions 
		long nCollTrans = atoi(chLine);

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		check_LAMDA_comment(chLine,chEFilename);

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
		// Get the number of temperatures
		long intCollTemp = atoi(chLine);

		// Read and store the table of rate coefficients
		ReadCollisionRateTable( AtmolCollRateCoeff[intNS][intCollIndex],
			ioLevData, func, nMolLevs, intCollTemp, nCollTrans );
	}
	delete func;

	fclose( ioLevData );

	return;
}

