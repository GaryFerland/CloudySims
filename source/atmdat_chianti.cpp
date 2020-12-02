/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "taulines.h"
#include "trace.h"
#include "thirdparty.h"
#include "atmdat.h"
#include "lines_service.h"
#include "parse_species.h"

typedef vector< pair<double,long> > DoubleLongPairVector;

/**
   AulTTError prints the error message which explains that a transition already has an Aul value for a particular
   transition type. If this occurs, the code will stop because there is probably an error in the data file.
 \param chFilename[]
 \param chLine[]
 \param TT[]
 */
static void AulTTError( const char chFilename[], const char chLine[], const char TT[] )
{
	DEBUG_ENTRY( "AulTTError()" );

	fprintf( ioQQQ, " PROBLEM File %s contains an invalid line.\n",chFilename);
	fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
	fprintf( ioQQQ, " This transition already has an Aul value set for %s.\n",TT);
	cdEXIT(EXIT_FAILURE);

}

static const bool DEBUGSTATE = false;
// minimum energy difference (wavenumbers) we will accept
const double ENERGY_MIN_WN = 1e-10;

void atmdat_STOUT_readin( long intNS, char *chPrefix )
{
	DEBUG_ENTRY( "atmdat_STOUT_readin()" );

	long int nMolLevs;
	char chLine[FILENAME_PATH_LENGTH_2],
		chNRGFilename[FILENAME_PATH_LENGTH_2],
		chTPFilename[FILENAME_PATH_LENGTH_2],
		chCOLLFilename[FILENAME_PATH_LENGTH_2];

	static const int MAX_NUM_LEVELS = 999;

	dBaseSpecies[intNS].lgMolecular = false;
	dBaseSpecies[intNS].lgLTE = false;

	strcpy( chNRGFilename , chPrefix );
	strcpy( chTPFilename , chPrefix );
	strcpy( chCOLLFilename , chPrefix );

	/*Open the energy levels file*/
	strcat( chNRGFilename , ".nrg");
	uncaps( chNRGFilename );

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_STOUT_readin opening %s:",chNRGFilename);

	FILE *ioDATA;
	ioDATA = open_data( chNRGFilename, "r" );
	bool lgEOL = false;
	long index = 0;
	double nrg = 0.0;
	double stwt = 0.0;

	/* first line is a version number - now confirm that it is valid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of %s.\n", chNRGFilename );
		cdEXIT(EXIT_FAILURE);
	}
	long int ipFFmt = 1;
	const long int iyr = 11, imo=10 , idy = 14;
	long iyrread, imoread , idyread;
	iyrread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
	imoread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
	idyread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);

	if(( iyrread != iyr ) ||
	  (  imoread != imo ) ||
	  (  idyread != idy ) )
	{
		fprintf( ioQQQ,
			" PROBLEM atmdat_STOUT_readin: the version of %s is not the current version.\n", chNRGFilename );
		fprintf( ioQQQ,
			" atmdat_STOUT_readin: I expected the magic numbers %li %li %li but found %li %li %li.\n",
			iyr, imo , idy ,iyrread, imoread , idyread  );
		cdEXIT(EXIT_FAILURE);
	}

	nMolLevs = 0;
	//Count number of levels
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* # - comment, *** ends data */
		if( chLine[0] != '#' && chLine[0] != '\n' && chLine[0] != '*' )
		{
			ipFFmt = 1;
			long n = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
			if( n < 0 )
				break;
			nMolLevs++;
		}
		else if( (chLine[0] == '*' && chLine[1] == '*' ) )
		{
			/* stop reading when field of stars encountered.*/
			break;
		}
	}

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not rewind %s.\n", chNRGFilename );
		cdEXIT(EXIT_FAILURE);
	}
	//Skip the magic numbers this time
	read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );

	long HighestIndexInFile = nMolLevs;

	dBaseSpecies[intNS].numLevels_max = HighestIndexInFile;

	setProperties(dBaseSpecies[intNS]);

	if( tolower(dBaseSpecies[intNS].chLabel[0]) == 'f' && tolower(dBaseSpecies[intNS].chLabel[1]) == 'e')
	{
		// Fe is special case with more levels
		nMolLevs = MIN3(nMolLevs, atmdat.nStoutMaxLevelsFe, MAX_NUM_LEVELS );
	}
	else
	{
		nMolLevs = MIN3(nMolLevs, atmdat.nStoutMaxLevels, MAX_NUM_LEVELS );
	}

	//Consider the masterlist specified number of levels as the min. =1 if not specified
	long numMasterlist = MIN2( dBaseSpecies[intNS].numLevels_masterlist , HighestIndexInFile );
	nMolLevs = MAX2(nMolLevs,numMasterlist);

	if (dBaseSpecies[intNS].setLevels != -1)
	{
		if (dBaseSpecies[intNS].setLevels > HighestIndexInFile)
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
			fprintf( ioQQQ,"Using STOUT spectrum %s (species: %s) with %li requested, only %li energy levels available.\n",
				dBaseSpecies[intNS].chLabel, chLabelChemical, dBaseSpecies[intNS].setLevels, HighestIndexInFile );
			nMolLevs = HighestIndexInFile;		  
		}
		else
		{
			nMolLevs = dBaseSpecies[intNS].setLevels;
		}
	}

	if( atmdat.lgStoutPrint == true)
	{
		char chLabelChemical[CHARS_SPECIES] = "";
		spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
		fprintf( ioQQQ,"Using STOUT spectrum %s (species: %s) with %li levels of %li available.\n",
				dBaseSpecies[intNS].chLabel, chLabelChemical, nMolLevs, HighestIndexInFile );
	}

	dBaseSpecies[intNS].numLevels_max = nMolLevs;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	/*Resize the States array*/
	dBaseStates[intNS].init(dBaseSpecies[intNS].chLabel,nMolLevs);
	/*Zero out the maximum wavenumber for each species */
	dBaseSpecies[intNS].maxWN = 0.;

	/* allocate the Transition array*/
	ipdBaseTrans[intNS].reserve(nMolLevs);
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
		ipdBaseTrans[intNS].reserve(ipHi,ipHi);
	ipdBaseTrans[intNS].alloc();
	dBaseTrans[intNS].resize(ipdBaseTrans[intNS].size());
	dBaseTrans[intNS].states() = &dBaseStates[intNS];
	dBaseTrans[intNS].chLabel() = dBaseSpecies[intNS].chLabel;
	dBaseSpecies[intNS].database = "Stout";

	//This is creating transitions that we don't have collision data for
	//Maybe use gbar or keep all of the Fe 2 even if it was assumed (1e-5)
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

	/* Create arrays for holding energies and statistical weights so that
	 * we can put the energies in the correct order before moving them to
	 * dBaseStates */
	DoubleLongPairVector dBaseStatesEnergy;
	vector<double> dBaseStatesStwt(HighestIndexInFile,-1.0);
	for( long ii = 0; ii < HighestIndexInFile; ii++ )
	{
		dBaseStatesEnergy.push_back(make_pair(-1.0,ii));
	}

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nStout Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Energy Level File: %s\n",chNRGFilename);
		fprintf(ioQQQ,"Number of Energy Levels in File: %li\n",HighestIndexInFile);
		fprintf(ioQQQ,"Number of Energy Levels Cloudy is Currently Using: %li\n",nMolLevs);
		fprintf(ioQQQ,"Species|File Index|Energy(WN)|StatWT\n");
	}

	/* Check for an end of file sentinel */
	bool lgSentinelReached = false;

	//Read the first line of data
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of the energy level file.\n");
		cdEXIT(EXIT_FAILURE);
	}
	//Read the remaining lines of the energy level file
	do
	{
		ipFFmt = 1;

		/* Stop on *** */
		if( chLine[0] == '*' )
		{
			lgSentinelReached = true;
			break;
		}
		//Comments start with #, skip them
		if( chLine[0] != '#' )
		{
			/* Skip blank lines */
			if( chLine[0] == '\n')
				continue;

			index = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL) -1 ;
			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"<%s>\t%li\t",dBaseSpecies[intNS].chLabel,index+1);
			}

			if( index < 0 )
			{
				fprintf( ioQQQ, " PROBLEM An energy level index was less than 1 in file %s\n",chNRGFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}

			nrg = (double)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
			stwt = (double)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"%.3f\t%.1f\n",nrg,stwt);
			}

			if( lgEOL )
			{
				fprintf( ioQQQ, " PROBLEM End of line reached prematurely in file %s\n",chNRGFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}

			/* Verify that energy levels are not overwritten */
			if( dBaseStatesEnergy.at(index).first == -1. && dBaseStatesStwt.at(index) == -1 )
			{
				dBaseStatesEnergy.at(index) = make_pair(nrg,index);
				dBaseStatesStwt.at(index) = stwt;
			}
			else
			{
				fprintf( ioQQQ, " PROBLEM Duplicate energy level index in file %s\n",chNRGFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL );
	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n",chNRGFilename);
		fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data and start in the first column of that line.\n");
		cdEXIT(EXIT_FAILURE);
	}
	fclose(ioDATA);

	//Sort levels by energy (then by the index in the file if necessary)
	sort(dBaseStatesEnergy.begin(),dBaseStatesEnergy.end());

	std::vector<long> indexold2new(dBaseStatesEnergy.size());
	for( DoubleLongPairVector::const_iterator i = dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
	{
		indexold2new[i->second] = i-dBaseStatesEnergy.begin();
	}

	
	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\n\n***Energy levels have been sorted in order of increasing energy.***\n");
		fprintf(ioQQQ,"Species|File Index|Sorted Index|Energy(WN)|StatWT\n");
	}

	/* Store sorted energies in their permanent home */
	for(DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin();i!=dBaseStatesEnergy.end();i++)
	{
		long oldindex = i->second;
		long index = i - dBaseStatesEnergy.begin();
		double nrg = i->first;
		double stwt = dBaseStatesStwt.at(oldindex);

		if( index >= nMolLevs )
			break;

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"<%s>\t%li\t%li\t%.3f\t%.1f\n",dBaseSpecies[intNS].chLabel,oldindex+1,index+1,nrg,stwt);
		}

		dBaseStates[intNS][index].energy().set(nrg,"cm^-1");
		dBaseStates[intNS][index].g() = stwt;
	}


	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	double fenergyWN = 0.;
	for(TransitionList::iterator tr=dBaseTrans[intNS].begin();
		 tr!= dBaseTrans[intNS].end(); ++tr)
	{
		int ipHi = (*tr).ipHi();
		int ipLo = (*tr).ipLo();
		ASSERT(ipHi > ipLo );
		fenergyWN = dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN();
		if( fenergyWN <= 0.)
		{
			long ipLoFile = dBaseStatesEnergy[ipLo].second;
			long ipHiFile = dBaseStatesEnergy[ipHi].second;


			if( DEBUGSTATE )
			{
				if( fenergyWN == 0. && dBaseStates[intNS][ipHi].g() != dBaseStates[intNS][ipLo].g() && abs(ipHiFile - ipLoFile) == 1)
				{
					/* This likely means that the data source lists adjacent energy levels with the same energy.
					 * Will just use the minimum energy. */
					if( atmdat.lgStoutPrint )
					{
						fprintf( ioQQQ, "\nCaution: A %s transition between adjacent energy levels has zero energy.\n",dBaseSpecies[intNS].chLabel);
						fprintf( ioQQQ, "The levels may have been sorted by energy since being read from the data files.\n");
						fprintf( ioQQQ, "The sorted lower and upper levels are %i and %i.\n",ipLo+1,ipHi+1);
						fprintf( ioQQQ, "The level indices as they appear in the Stout energy level data file, %s, are %li and %li.\n",chNRGFilename,ipLoFile+1,ipHiFile+1);
						fprintf( ioQQQ, "Verify with the data source that the correct energies are listed in the Stout data file.\n");
					}

				}
				else
				{
					fprintf( ioQQQ, "The levels may have been sorted by energy since being read from the data files.\n");
					fprintf( ioQQQ, "The sorted lower and upper levels are %i and %i.\n",ipLo+1,ipHi+1);
					fprintf( ioQQQ, "The level indices as they appear in the Stout energy level data file, %s, are %li and %li.\n",chNRGFilename,ipLoFile+1,ipHiFile+1);
					fprintf( ioQQQ, "Check the data file for missing or duplicate energies.\n");
					//cdEXIT(EXIT_FAILURE);
				}
			}
		}
		(*tr).EnergyWN() = (realnum)MAX2(ENERGY_MIN_WN,fenergyWN);
		(*tr).WLAng() = (realnum)(1e+8/(*tr).EnergyWN()/RefIndex((*tr).EnergyWN()));
		dBaseSpecies[intNS].maxWN = MAX2(dBaseSpecies[intNS].maxWN,(*tr).EnergyWN());
	}

	/******************************************************
	 ************* Transition Probability File ************
	 ******************************************************/
	strcat( chTPFilename , ".tp");
	uncaps( chTPFilename );

	ioDATA = open_data( chTPFilename, "r" );

	/* first line is a version number - now confirm that it is valid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of %s.\n", chTPFilename );
		cdEXIT(EXIT_FAILURE);
	}
	ipFFmt = 1;
	iyrread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
	imoread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
	idyread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);

	if(( iyrread != iyr ) ||
	  (  imoread != imo ) ||
	  (  idyread != idy ) )
	{
		fprintf( ioQQQ,
			" PROBLEM atmdat_STOUT_readin: the version of %s is not the current version.\n", chTPFilename );
		fprintf( ioQQQ,
			" atmdat_STOUT_readin: I expected the magic numbers %li %li %li but found %li %li %li.\n",
			iyr, imo , idy ,iyrread, imoread , idyread  );
		cdEXIT(EXIT_FAILURE);
	}

	long numtrans = 0;
	//Count number of transitions
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* # is comment, *** is end of data */
		if( chLine[0] != '#' && chLine[0] != '\n' && chLine[0] != '*' )
		{
			ipFFmt = 1;
			long n = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
			if( n < 0 )
				break;
			numtrans++;
		}
		else if( (chLine[0] == '*' && chLine[1] == '*' ) )
		{
			/* stop reading when field of stars encountered.*/
			break;
		}
	}
	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not rewind %s.\n", chTPFilename );
		cdEXIT(EXIT_FAILURE);
	}
	//Skip the magic numbers this time
	read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );


	//Read the first line of data
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of the transition probability file.\n");
		cdEXIT(EXIT_FAILURE);
	}

	double tpData = 0.0;
	lgSentinelReached = false;
	static const int intNumCols = 6;

	/* lgLineStrengthTT functions as a checklist for line strengths
	 * from the various transition types (E1,M2, etc.)
	 * When a line strength is added to the total Aul from a specific
	 * transition type, the corresponding value (based on ipLo, ipHi, and transition type)
	 * is set to true.
	 * lgLineStrengthTT[ipLo][ipHi][k]:
	 * k = 0 => E1
	 * k = 1 => E2
	 * k = 2 => E3
	 * k = 3 => M1
	 * k = 4 => M2
	 * k = 5 => M3
	 */
	bool ***lgLineStrengthTT;
	lgLineStrengthTT = (bool ***)MALLOC(nMolLevs *sizeof(bool**));
	for( int ii = 0; ii < nMolLevs; ii++ )
	{
		lgLineStrengthTT[ii] = (bool **)MALLOC(nMolLevs *sizeof(bool*));
		for( int j = 0; j < nMolLevs; j++ )
		{
			lgLineStrengthTT[ii][j] = (bool *)MALLOC(intNumCols *sizeof(bool));
		}
	}

	/* Initialize lgLineStrengthTT values to false */
	for( int ii = 0; ii < nMolLevs; ii++ )
	{
		for( int j = 0; j < nMolLevs; j++ )
		{
			for( int k = 0; k < intNumCols; k++ )
			{
				lgLineStrengthTT[ii][j][k] = false;
			}
		}
	}

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nStout Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Radiative Data File: %s\n",chTPFilename);
		fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Data Type (A,G,S)|Data\n");
	}

	//Read the remaining lines of the transition probability file
	do
	{
		if( chLine[0] == '*' )
		{
			lgSentinelReached = true;
			break;
		}

		//Comments start with #, skip them
		if( chLine[0] != '#' )
		{
			/* skip null lines */
			if( chLine[0] == '\n')
				continue;

			if( nMatch("A",chLine) || nMatch("G",chLine) || nMatch("S",chLine) )
			{
				/* reset read pointer */
				ipFFmt = 1;

				/* save original level for reference, array of levels will be sorted
				 * to become energy ordered */
				long ipLoInFile = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
				long ipHiInFile = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
				long ipLo = ipLoInFile - 1;
				long ipHi = ipHiInFile - 1;
				tpData = (double)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);

				if( tpData < atmdat.aulThreshold && nMatch("A",chLine))
				{
					// skip these lines
					continue;
				}

				/* Account for reordered energy levels */
				if (ipHi >= long(indexold2new.size()))
				{
					continue;
				}
				ipLo = indexold2new[ipLo];
				ipHi = indexold2new[ipHi];
				
				if( ipLo < 0 || ipLo >= nMolLevs || ipHi < 0 || ipHi >= nMolLevs )
				{
					// skip these lines
					continue;
				}

				ASSERT( ipLo != ipHi );

				// swap indices if energy levels were not correctly sorted
				if( ipHi < ipLo )
				{
					long swap = ipHi;
					ipHi = ipLo;
					ipLo = swap;
				}

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"<%s>\t%li:%li\t%li:%li\t%c\t%.2e\n",
						dBaseSpecies[intNS].chLabel,ipLoInFile,ipHiInFile,
						ipLo+1,ipHi+1,chLine[0],tpData);
				}

				if( lgEOL )
				{
					fprintf( ioQQQ, " PROBLEM End of line reached prematurely in file %s\n",chTPFilename);
					fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
					cdEXIT(EXIT_FAILURE);
				}

				TransitionList::iterator tr = dBaseTrans[intNS].begin()+ipdBaseTrans[intNS][ipHi][ipLo];

				/* If we don't already have the transition in the stack, add it and
				 * zero out Aul */
				if( !(*tr).hasEmis() )
				{
					(*tr).AddLine2Stack();
					(*tr).Emis().Aul() = 0.;
					(*tr).Emis().gf() = 0.;
				}

				//This means last data column has Aul.
				if( nMatch("A",chLine) )
				{
					if( (*tr).EnergyWN() > ENERGY_MIN_WN )
					{
						(*tr).Emis().Aul() += tpData;
						// use updated total Aul to get gf
						(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
					}
				}
				//This means last data column has gf.
				else if( nMatch("G",chLine) )
				{
					if( (*tr).EnergyWN() > ENERGY_MIN_WN )
					{
						(*tr).Emis().gf() += tpData;
						// use updated total gf to get Aul
						(*tr).Emis().Aul() = (realnum)eina((*tr).Emis().gf(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
					}
				}
				else if( nMatch("S",chLine) )
				{
					static const double BOHR_MAGNETON = ELEM_CHARGE_ESU*H_BAR/2/ELECTRON_MASS/SPEEDLIGHT;
					/** Bohr Magneton, 9.2740096e-21 ergs/G */

					/* Data column is composed of line strengths */
					if( nMatch("E1",chLine) )
					{
						if( lgLineStrengthTT[ipLo][ipHi][0] )
						{
							AulTTError(chTPFilename,chLine,"E1\0");
						}

						/*Convert line strength to Aul for E1 transitions
						 * Aul = 64*Pi^4*e^2*a0^2/3/h*S/gu/WLAng^3  */
						static const double E1Coeff = 64*powi(PI,4)*pow(ELEM_CHARGE_ESU,2)*pow2(BOHR_RADIUS_CM)/3/HPLANCK/pow3(1e-8);
						/* E1Coeff = 2.0261e18 */

						(*tr).Emis().Aul() += E1Coeff*tpData/(*(*tr).Hi()).g()/pow3((*tr).EnergyAng());
						(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
						lgLineStrengthTT[ipLo][ipHi][0]= true;
					}
					else if( nMatch("E2",chLine) )
					{
						if( lgLineStrengthTT[ipLo][ipHi][1] )
						{
							AulTTError(chTPFilename,chLine,"E2\0");
						}

						/*Convert line strength to Aul for E2 transitions
						 * Aul = 64*Pi^6*e^2*a0^4/15/h*S/gu/WLAng^5  */
						static const double E2Coeff = 64*powi(PI,6)*pow2(ELEM_CHARGE_ESU)*pow4(BOHR_RADIUS_CM)/15/HPLANCK/powi(1e-8,5);
						/* E2Coeff = 1.1199e18 */

						(*tr).Emis().Aul() += E2Coeff*tpData/(*(*tr).Hi()).g()/powi((*tr).EnergyAng(),5);
						(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
						lgLineStrengthTT[ipLo][ipHi][1]= true;
					}
					else if( nMatch("E3",chLine) )
					{
						if( lgLineStrengthTT[ipLo][ipHi][2] )
						{
							AulTTError(chTPFilename,chLine,"E3\0");
						}

						/*Convert line strength to Aul for E3 transitions
						 * Aul = 2048*Pi^8*e^2*a0^6/4725/h*S/gu/WLAng^7 */
						static const double E3Coeff = 2048*powi(PI,8)*pow2(ELEM_CHARGE_ESU)*powi(BOHR_RADIUS_CM,6)/4725/HPLANCK/powi(1e-8,7);
						/* E3Coeff = 3.1444165e17
						 * Atomic Transition Probabilites of Silicon. A Critical Compilation
						 * Kelleher & Podobedova 2006*/

						(*tr).Emis().Aul() += E3Coeff*tpData/(*(*tr).Hi()).g()/powi((*tr).EnergyAng(),7);
						(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
						lgLineStrengthTT[ipLo][ipHi][2]= true;
					}
					else if( nMatch("M1",chLine) )
					{
						if( lgLineStrengthTT[ipLo][ipHi][3] )
						{
							AulTTError(chTPFilename,chLine,"M1\0");
						}

						/*Convert line strength to Aul for M1 transitions
						 * Aul = 64*Pi^4*mu_B^2/3/h*S/gu/WLAng^3 */
						static const double M1Coeff = 64*powi(PI,4)*pow2(BOHR_MAGNETON)/3/HPLANCK/pow3(1e-8);
						/* M1Coeff = 2.697e13 */

						(*tr).Emis().Aul() += M1Coeff*tpData/(*(*tr).Hi()).g()/pow3((*tr).EnergyAng());
						(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
						lgLineStrengthTT[ipLo][ipHi][3]= true;
					}
					else if( nMatch("M2",chLine) )
					{
						if( lgLineStrengthTT[ipLo][ipHi][4] )
						{
							AulTTError(chTPFilename,chLine,"M2\0");
						}

						/*Convert line strength to Aul for M2 transitions
						 * Aul = 64*Pi^6*mu_B^2*a0^2/15/h*S/gu/WLAng^5 */
						static const double M2Coeff = 64*powi(PI,6)*pow2(BOHR_MAGNETON)*pow2(BOHR_RADIUS_CM)/15/HPLANCK/powi(1e-8,5);
						/* M2Coeff = 1.4909714e13
						 * Tables of Atomic Transition Probabilities for Be and B
						* Fuhr & Wiese  2010*/

						(*tr).Emis().Aul() += M2Coeff*tpData/(*(*tr).Hi()).g()/powi((*tr).EnergyAng(),5);
						(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
						lgLineStrengthTT[ipLo][ipHi][4]= true;
					}
					else if( nMatch("M3",chLine) )
					{
						if( lgLineStrengthTT[ipLo][ipHi][5] )
						{
							AulTTError(chTPFilename,chLine,"M3\0");
						}

						/*Convert line strength to Aul for M3 transitions
						 * Aul = 2048*Pi^8*mu_B^2*a0^4/4725/h*S/gu/WLAng^7 */
						static const double M3Coeff = 2048*powi(PI,8)*pow2(BOHR_MAGNETON)*pow4(BOHR_RADIUS_CM)/4725/HPLANCK/powi(1e-8,7);
						/* M2Coeff = 4.18610e12
						* Safronova & Safronova et al 2005*/


						(*tr).Emis().Aul() += M3Coeff*tpData/(*(*tr).Hi()).g()/powi((*tr).EnergyAng(),7);
						(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
						lgLineStrengthTT[ipLo][ipHi][5]= true;
					}
					else
					{
						fprintf( ioQQQ, " PROBLEM File %s contains an invalid line.\n",chTPFilename);
						fprintf( ioQQQ, " The line strength does not list a valid transition type.\n");
						fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
						cdEXIT(EXIT_FAILURE);
					}
				}
				else
				{
					/* The code has already passed a check looking for A, G or S.
					 * Now it can't find any of the three. Insanity. */
					TotalInsanity();
				}

				(*tr).setComment( db_comment_tran_levels( ipLoInFile, ipHiInFile ) );
			}
			else
			{
				fprintf( ioQQQ, " PROBLEM File %s contains an invalid line.\n",chTPFilename);
				fprintf( ioQQQ, " Data must either be in the form of Aul, gf, or S(line strength) and list"
						" the data type in the first colum as A, G, or S respectively.");
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL );
	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n",chTPFilename);
		fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data and start in the first column of that line.");
		cdEXIT(EXIT_FAILURE);
	}
	fclose(ioDATA);

	/******************************************************
	 ************* Collision Data File ********************
	 ******************************************************/

	strcat( chCOLLFilename , ".coll");
	uncaps( chCOLLFilename );

	ioDATA = open_data( chCOLLFilename, "r" );

	/* first line is a version number - now confirm that it is valid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of %s.\n", chCOLLFilename );
		cdEXIT(EXIT_FAILURE);
	}
	ipFFmt = 1;
	iyrread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
	imoread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
	idyread = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);

	if(( iyrread != iyr ) ||
	  (  imoread != imo ) ||
	  (  idyread != idy ) )
	{
		fprintf( ioQQQ,
			" PROBLEM atmdat_STOUT_readin: the version of %s is not the current version.\n", chCOLLFilename );
		fprintf( ioQQQ,
			" atmdat_STOUT_readin: I expected the magic numbers %li %li %li but found %li %li %li.\n",
			iyr, imo , idy ,iyrread, imoread , idyread  );
		cdEXIT(EXIT_FAILURE);
	}

	/****** Could add ability to count number of temperature changes, electron CS, and proton CS ****/


	//Read the first line of data
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of the collision data file.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* Malloc space for collision strengths */
	StoutCollData[intNS].alloc(nMolLevs,nMolLevs,ipNCOLLIDER);
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			for( long k=0; k<ipNCOLLIDER; k++ )
			{
				/* initialize all spline variables */
				StoutCollData[intNS].junk(ipHi,ipLo,k);
			}
		}
	}

	int numpoints = 0;
	vector<double> temps;
	long ipCollider = -1;
	lgSentinelReached = false;

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nStout Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Collision Data File: %s\n",chCOLLFilename);
		fprintf(ioQQQ,"Species|TEMP|Temperatures (K)\n");
		fprintf(ioQQQ,"Species|Data Type (CS,RATE)|Collider|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Data\n");
	}

	//Read the remaining lines of the collision data file
	do
	{
		/* Stop on *** */
		if( chLine[0] == '*' )
		{
			lgSentinelReached = true;
			break;
		}

		//Comments start with #, skip them
		if( chLine[0] != '#' )
		{
			ipFFmt = 1;

			/* Skip blank lines */
			if( chLine[0] == '\n')
				continue;

			/* Make all letters of the line upper case so they are case insensitive */
			caps( chLine );

			//This is a temperature line
			if( nMatch("TEMP",chLine) )
			{

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"<%s>\tTEMP\t",dBaseSpecies[intNS].chLabel);
				}

				// count number of temperature points
				ipFFmt = 1;
				numpoints = 0;
				while( !lgEOL )
				{
					FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
					if( lgEOL )
						break;
					numpoints++;
				}
				ASSERT( numpoints > 0 );

				temps.resize(numpoints);
				for( int j = 0; j < numpoints; j++ )
				{
					temps[j] = 0.;
				}

				ipFFmt = 1;
				for( int j = 0; j < numpoints; j++ )
				{
					temps[j] = (double)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"%.2e\t",temps[j]);
					}
					ASSERT( temps[j] > 0 );
				}
				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"\n");
				}
			}
			else if( nMatch("CS",chLine) || nMatch("RATE",chLine) )
			{

				bool isRate = false;
				ipFFmt = 1;
				if( nMatch("RATE", chLine) )
					isRate = true;

				if( nMatch( "ELECTRON",chLine ) )
				{
					ipCollider = ipELECTRON;
				}
				else if( nMatch( "PROTON",chLine ) || nMatch( "H+",chLine ) )
				{
					ipCollider = ipPROTON;
				}
				else if( nMatch( "HE+2", chLine )  )
				{
					ipCollider = ipALPHA;
					//Move the cursor past the number in the species definition
					(void)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
				}
				else if( nMatch( "HE+ ",chLine ) || nMatch( "HE+\t", chLine) )
				{
					ipCollider = ipHE_PLUS;
				}
				else if( nMatch( "H2 ",chLine ) || nMatch( "H2\t",chLine ) )
				{
					//Move the cursor past the number in the species definition
					(void)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);

					if( nMatch( "ORTHO",chLine ) )
					{
						ipCollider = ipH2_ORTHO;
					}
					else if( nMatch( "PARA",chLine ) )
					{
						ipCollider = ipH2_PARA;
					}
					else
					{
						ipCollider = ipH2;
					}
				}
				else if( nMatch( "HE ",chLine ) || nMatch( "HE\t",chLine ) )
				{
					ipCollider = ipATOM_HE;
				}
				else if( nMatch( "H ",chLine ) || nMatch( "H\t",chLine ) )
				{
					ipCollider = ipATOM_H;
				}
				else
				{
					fprintf( ioQQQ,	" PROBLEM atmdat_STOUT_readin: Each line of the collision data"
							"file must specify the collider.\n");
					fprintf( ioQQQ,	" Possible colliders are: ELECTRON, PROTON, HE, H,"
							" HE+2, HE+, H2, H2 ORTHO, H2 PARA\n");
					cdEXIT(EXIT_FAILURE);
				}

				if( temps.empty() )
				{
					fprintf( ioQQQ,	" PROBLEM atmdat_STOUT_readin: The collision "
							"file must specify temperatures before the collision strengths");
					cdEXIT(EXIT_FAILURE);
				}

				long ipLoInFile = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
				long ipHiInFile = (long)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
				long ipLo = ipLoInFile-1;
				long ipHi = ipHiInFile-1;

				/* Account for reordered energy levels */
				if (ipHi >= long(indexold2new.size()))
				{
					continue;
				}
				ipLo = indexold2new[ipLo];
				ipHi = indexold2new[ipHi];

				if( ipLo < 0 || ipLo >= nMolLevs || ipHi < 0 || ipHi >= nMolLevs )
				{
					// skip these lines
					continue;
				}

				ASSERT( ipLo != ipHi );

				// swap indices if energy levels were not correctly sorted
				if( ipHi < ipLo )
				{
					long swap = ipHi;
					ipHi = ipLo;
					ipLo = swap;
				}

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"<%s>\t%s\t%li\t%li:%li\t%li:%li",
						dBaseSpecies[intNS].chLabel,isRate?"RATE":"CS",ipCollider,
						ipLoInFile,ipHiInFile,ipLo+1,ipHi+1);
				}

				/* Set this as a collision strength not a collision rate coefficient*/
				StoutCollData[intNS].lgIsRate(ipHi,ipLo,ipCollider) = isRate;

				ASSERT( numpoints > 0 );
				StoutCollData[intNS].setpoints(ipHi,ipLo,ipCollider,numpoints);

				/* Loop over all but one CS value */
				for( int j = 0; j < numpoints; j++ )
				{
					StoutCollData[intNS].temps(ipHi,ipLo,ipCollider)[j] = temps[j];
					StoutCollData[intNS].collstrs(ipHi,ipLo,ipCollider)[j] = (double)FFmtRead(chLine,&ipFFmt,sizeof(chLine),&lgEOL);
					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"\t%.2e",StoutCollData[intNS].collstrs(ipHi,ipLo,ipCollider)[j]);
					}
					if( StoutCollData[intNS].temps(ipHi,ipLo,ipCollider)[j] <= 0 ||
						 StoutCollData[intNS].collstrs(ipHi,ipLo,ipCollider)[j] <= 0 )
					{
						fprintf(ioQQQ,"PROBLEM: A Stout temperature or collision strength is less than or equal to zero.\n");
						fprintf(ioQQQ,"Species = %s\tipLo = %li\tipHi = %li\n",dBaseSpecies[intNS].chLabel,ipLo+1,ipHi+1);
						fprintf(ioQQQ,"numpoint = %i\tTemp = %e\tCS = %e\n",j,temps[j],
								  StoutCollData[intNS].collstrs(ipHi,ipLo,ipCollider)[j]);
						cdEXIT(EXIT_FAILURE);
					}
				}
				if( DEBUGSTATE )
					fprintf(ioQQQ,"\n");
			}
			else
			{
				fprintf( ioQQQ, " PROBLEM File %s contains an invalid line.\n",chCOLLFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}

			if( lgEOL )
			{
				fprintf( ioQQQ, " PROBLEM End of line reached prematurely in file %s\n",chCOLLFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL );
	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n",chCOLLFilename);
		fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data and start in the first column.");
		cdEXIT(EXIT_FAILURE);
	}
	fclose(ioDATA);

	/* Free up memory from lgLineStrengthTT */
	for( int ii = 0; ii < nMolLevs; ii++ )
	{
		for( int j = 0; j < nMolLevs; j++ )
		{
			free( lgLineStrengthTT[ii][j] );
		}
		free( lgLineStrengthTT[ii] );
	}
	free( lgLineStrengthTT );

	return;
}

void atmdat_CHIANTI_readin( long intNS, char *chPrefix )
{
	DEBUG_ENTRY( "atmdat_CHIANTI_readin()" );

	int intsplinepts,intTranType,intxs;
	long int nMolLevs,nMolExpLevs,nElvlcLines,nTheoLevs;// number of experimental and total levels
	FILE *ioElecCollData=NULL, *ioProtCollData=NULL;
	realnum  fstatwt,fenergyWN,fWLAng,fenergy,feinsteina;
	double fScalingParam,fEnergyDiff;
	const char chCommentChianti = '#';

	char chLine[FILENAME_PATH_LENGTH_2] ,
		chEnFilename[FILENAME_PATH_LENGTH_2],
		chTraFilename[FILENAME_PATH_LENGTH_2],
		chEleColFilename[FILENAME_PATH_LENGTH_2],
		chProColFilename[FILENAME_PATH_LENGTH_2];

	bool lgProtonData=false;

	// this is the largest number of levels allowed by the chianti format, I3
	static const int MAX_NUM_LEVELS = 999;

	dBaseSpecies[intNS].lgMolecular = false;
	dBaseSpecies[intNS].lgLTE = false;

	strcpy( chEnFilename , chPrefix );
	strcpy( chTraFilename , chPrefix );	
	strcpy( chEleColFilename , chPrefix );		
	strcpy( chProColFilename , chPrefix );			

	/*For the CHIANTI DATABASE*/
	/*Open the energy levels file*/
	strcat( chEnFilename , ".elvlc");
	uncaps( chEnFilename );

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEnFilename);

	fstream elvlcstream;
	open_data(elvlcstream, chEnFilename,mode_r);

	/*Open the transition probabilities file*/
	strcat( chTraFilename , ".wgfa");
	uncaps( chTraFilename );

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chTraFilename);

	fstream wgfastream;
	open_data(wgfastream, chTraFilename,mode_r);

	/*Open the electron collision data*/
	strcat( chEleColFilename , ".splups");
	uncaps( chEleColFilename );

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEleColFilename);

	ioElecCollData = open_data( chEleColFilename, "r" );

	/*Open the proton collision data*/
	strcat( chProColFilename , ".psplups");
	uncaps( chProColFilename );

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chProColFilename);

	/*We will set a flag here to indicate if the proton collision strengths are available */
	if( ( ioProtCollData = open_data( chProColFilename, "r", AS_DATA_ONLY_TRY ) ) != NULL )
	{
		lgProtonData = true;
	}
	else
	{
		lgProtonData = false;
	}

	/*Loop finds how many theoretical and experimental levels are in the elvlc file */
	//eof_col is used get the first 4 charcters per line to find end of file
	const int eof_col = 5;
	//length (+1) of the nrg in the elvlc file
	const int lvl_nrg_col=16;
	//# of columns skipped from the left to get to nrg start
	const int lvl_skipto_nrg = 40;
	/* # of columns to skip from eof check to nrg start */
	const int lvl_eof_to_nrg = lvl_skipto_nrg - eof_col + 1;
	//# of columns to skip over the rydberg energy, we don't use it
	const int lvl_skip_ryd = 15;
	nElvlcLines = 0;
	nMolExpLevs = 1;
	nTheoLevs = 1;
	if (elvlcstream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		char exptemp[lvl_nrg_col],theotemp[lvl_nrg_col];
		double tempexpenergy = 0.,theoenergy = 0.;
		/*This loop counts the number of valid rows within the elvlc file
		  as well as the number of experimental energy levels.*/
		while(nj != -1)
		{
			elvlcstream.get(otemp,eof_col);
			nj = atoi(otemp);
			if( nj == -1)
				break;
			nElvlcLines++;

			elvlcstream.seekg(lvl_eof_to_nrg,ios::cur);
			elvlcstream.get(exptemp,lvl_nrg_col);
			tempexpenergy = (double) atof(exptemp);
			if( tempexpenergy != 0.)
				nMolExpLevs++;

			elvlcstream.seekg(lvl_skip_ryd,ios::cur);
			elvlcstream.get(theotemp,lvl_nrg_col);
			theoenergy = (double) atof(theotemp);
			if( theoenergy != 0. )
				nTheoLevs++;

			elvlcstream.ignore(INT_MAX,'\n');

		}
		elvlcstream.seekg(0,ios::beg);
	}

	//Sometimes the theoretical chianti level data is incomplete.
	//If it is bad use experimental
	bool lgChiaBadTheo = false;
	if( !atmdat.lgChiantiExp && nTheoLevs < nElvlcLines )
	{
		lgChiaBadTheo = true;
		atmdat.lgChiantiExp = true;
		fprintf(ioQQQ,"Warning: The theoretical energy levels for %s are incomplete.",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Switching to the experimental levels for this species.");
	}

	long HighestIndexInFile = -1;

	/* The total number of levels depends on the experimental Chianti switch */
	if( atmdat.lgChiantiExp )
	{
		HighestIndexInFile = nMolExpLevs;
	}
	else
	{
		HighestIndexInFile = nElvlcLines;
	}

	dBaseSpecies[intNS].numLevels_max = HighestIndexInFile;

	setProperties(dBaseSpecies[intNS]);
	
	if( tolower(dBaseSpecies[intNS].chLabel[0]) == 'f' && tolower(dBaseSpecies[intNS].chLabel[1]) == 'e')
	{
		// Fe is special case with more levels
		nMolLevs = MIN3(HighestIndexInFile, atmdat.nChiantiMaxLevelsFe,MAX_NUM_LEVELS );
	}
	else
	{
		nMolLevs = MIN3(HighestIndexInFile, atmdat.nChiantiMaxLevels,MAX_NUM_LEVELS );
	}

	if( nMolLevs <= 0 )
	{
		fprintf( ioQQQ, "The number of energy levels is non-positive in datafile %s.\n", chEnFilename );
		fprintf( ioQQQ, "The file must be corrupted.\n" );
		cdEXIT( EXIT_FAILURE );
	}

	//Consider the masterlist specified number of levels as the min. =1 if not specified
	long numMasterlist = MIN2( dBaseSpecies[intNS].numLevels_masterlist , HighestIndexInFile );
	nMolLevs = MAX2(nMolLevs,numMasterlist);

	if (dBaseSpecies[intNS].setLevels != -1)
	{
		if (dBaseSpecies[intNS].setLevels > HighestIndexInFile)
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
			fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li requested, only %li energy levels available.\n",
				dBaseSpecies[intNS].chLabel, chLabelChemical, dBaseSpecies[intNS].setLevels, HighestIndexInFile );
			nMolLevs = HighestIndexInFile;		  
		}
		else
		{
			nMolLevs = dBaseSpecies[intNS].setLevels;
		}
	}

	dBaseSpecies[intNS].numLevels_max = nMolLevs;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	if( atmdat.lgChiantiPrint == true)
	{
		if( atmdat.lgChiantiExp )
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
			fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li experimental energy levels of %li available.\n",
				dBaseSpecies[intNS].chLabel, chLabelChemical, nMolLevs , nMolExpLevs );
		}
		else
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
			fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li theoretical energy levels of %li available.\n",
				dBaseSpecies[intNS].chLabel, chLabelChemical, nMolLevs , nElvlcLines );
		}
	}

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
	dBaseSpecies[intNS].database = "Chianti";

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

	/*Keep track of which levels have experimental data and then create a vector
	which relates their indices to the default chianti energy indices.
	 */
	long ncounter = 0;

	//Relate Chianti level indices to a set that only include experimental levels
	vector<long> intExperIndex(nElvlcLines,-1);

	DoubleLongPairVector dBaseStatesEnergy;
	vector<double> dBaseStatesStwt(HighestIndexInFile,-1.0);
	for( long ii = 0; ii < HighestIndexInFile; ii++ )
	{
		dBaseStatesEnergy.push_back(make_pair(-1.0,ii));
	}

	//lvl_skipto_statwt is the # of columns to skip to statwt from left
	const int lvl_skipto_statwt = 37;
	//lvl_statwt_col is the length (+1) of statwt
	const int lvl_statwt_col = 4;
	//Read in stat weight and energy

	//Read in nrg levels to see if they are in order
	for( long ipLev=0; ipLev<nElvlcLines; ipLev++ )
	{
		if(elvlcstream.is_open())
		{
			char gtemp[lvl_statwt_col],thtemp[lvl_nrg_col],obtemp[lvl_nrg_col];
			elvlcstream.seekg(lvl_skipto_statwt,ios::cur);
			elvlcstream.get(gtemp,lvl_statwt_col);
			fstatwt = (realnum)atof(gtemp);
			elvlcstream.get(thtemp,lvl_nrg_col);
			fenergy = (double) atof(thtemp);

			if(fstatwt <= 0.)
			{
				fprintf( ioQQQ, " WARNING: A positive non zero value is expected for the "
						"statistical weight but was not found in %s"
						" level %li\n", chEnFilename,ipLev);
				cdEXIT(EXIT_FAILURE);
			}

			if( atmdat.lgChiantiExp )
			{
				/* Go through the entire level list selectively choosing only experimental level energies.
				 * Store them, not zeroes, in order using ncounter to count the index.
				 * Any row on the level list where there is no experimental energy, put a -1 in the relational vector.
				 * If it is a valid experimental energy level store the new ncounter index.
				 */

				if( fenergy != 0. || ipLev == 0 )
				{
					dBaseStatesEnergy.at(ncounter).first = fenergy;
					dBaseStatesEnergy.at(ncounter).second = ncounter;
					dBaseStatesStwt.at(ncounter) = fstatwt;
					intExperIndex.at(ipLev) = ncounter;
					ncounter++;
				}
				else
				{
					intExperIndex.at(ipLev) = -1;
				}
			}
			else
			{
				elvlcstream.seekg(lvl_skip_ryd,ios::cur);
				elvlcstream.get(obtemp,lvl_nrg_col);
				fenergy = (double) atof(obtemp);
				if(fenergy != 0. || ipLev == 0)
				{
					dBaseStatesEnergy.at(ipLev).first = fenergy;
					dBaseStatesEnergy.at(ipLev).second = ipLev;
					dBaseStatesStwt.at(ipLev) = fstatwt;
				}
				else
				{
					dBaseStatesEnergy.at(ipLev).first = -1.;
					dBaseStatesEnergy.at(ipLev).second = ipLev;
					dBaseStatesStwt.at(ipLev) = -1.;
				}
			}

			elvlcstream.ignore(INT_MAX,'\n');
		}
		else
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEnFilename);
			fclose( ioProtCollData );
			cdEXIT(EXIT_FAILURE);
		}
	}

	elvlcstream.close();

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nintExperIndex Vector:\n");
		fprintf(ioQQQ,"File Index|Exper Index\n");
		for( vector<long>::iterator i = intExperIndex.begin(); i != intExperIndex.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-intExperIndex.begin())+1;
			fprintf(ioQQQ,"%li\t%li\n",iPrt,(*i)+1);
		}

		for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-dBaseStatesEnergy.begin())+1;
			fprintf(ioQQQ,"PreSort:%li\t%li\t%f\t%f\n",iPrt,
					(i->second)+1,i->first,dBaseStatesStwt.at(i->second));
		}
	}

	//Sort energy levels
	sort(dBaseStatesEnergy.begin(),dBaseStatesEnergy.end());

	std::vector<long> indexold2new(dBaseStatesEnergy.size());
	for( DoubleLongPairVector::const_iterator i = dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
	{
		indexold2new[i->second] = i-dBaseStatesEnergy.begin();
	}

	if( DEBUGSTATE )
	{
		for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-dBaseStatesEnergy.begin())+1;
			if( iPrt > nMolLevs )
				break;
			fprintf(ioQQQ,"PostSort:%li\t%li\t%f\t%f\n",iPrt,
					(i->second)+1,i->first,dBaseStatesStwt.at(i->second));
		}

		fprintf(ioQQQ,"\nChianti Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Energy Level File: %s\n",chEnFilename);
		if( atmdat.lgChiantiExp )
		{
			fprintf(ioQQQ,"Number of Experimental Energy Levels in File: %li\n",nMolExpLevs);
		}
		else
		{
			fprintf(ioQQQ,"Number of Theoretical Energy Levels in File: %li\n",nElvlcLines);
		}
		fprintf(ioQQQ,"Number of Energy Levels Cloudy is Currently Using: %li\n",nMolLevs);
		fprintf(ioQQQ,"Species|File Index|Cloudy Index|StatWT|Energy(WN)\n");
	}

	vector<long> revIntExperIndex;
	if ( atmdat.lgChiantiExp )
	{
		revIntExperIndex.resize(dBaseStatesEnergy.size());
		for (size_t i = 0; i<dBaseStatesEnergy.size(); ++i)
			revIntExperIndex[i] = -1;
		for ( vector<long>::const_iterator i= intExperIndex.begin();
		      i != intExperIndex.end(); ++i )
		{
			long ipos = intExperIndex[i-intExperIndex.begin()];
			if (ipos >= 0 && ipos < long(dBaseStatesEnergy.size()))
			    revIntExperIndex[ipos] = i-intExperIndex.begin();
		}
	}
	
	for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
	{

		long ipLevNew = i - dBaseStatesEnergy.begin();
		long ipLevFile = -1;

		if( ipLevNew >= nMolLevs )
			break;

		if( atmdat.lgChiantiExp )
		{
			ipLevFile = revIntExperIndex[ipLevNew];
		}
		else
		{
			ipLevFile = i->second;
		}

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"<%s>\t%li\t%li\t",dBaseSpecies[intNS].chLabel,ipLevFile+1,ipLevNew+1);
		}

		dBaseStates[intNS][ipLevNew].g() = dBaseStatesStwt.at(i->second);
		dBaseStates[intNS][ipLevNew].energy().set(i->first,"cm^-1");

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"%.1f\t",dBaseStatesStwt.at(i->second));
			fprintf(ioQQQ,"%.3f\n",i->first);
		}
	}

	// highest energy transition in chianti
	dBaseSpecies[intNS].maxWN = 0.;
	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for(TransitionList::iterator tr=dBaseTrans[intNS].begin();
		 tr!= dBaseTrans[intNS].end(); ++tr)
	{
		int ipHi = (*tr).ipHi();
		int ipLo = (*tr).ipLo();
		fenergyWN = (realnum)MAX2( ENERGY_MIN_WN , dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN() );

		(*tr).EnergyWN() = fenergyWN;

		(*tr).WLAng() = (realnum)(1e+8/fenergyWN/RefIndex(fenergyWN));
		dBaseSpecies[intNS].maxWN = MAX2(dBaseSpecies[intNS].maxWN,fenergyWN);
	}

	/************************************************************************/
	/*** Read in the level numbers, Einstein A and transition wavelength  ***/
	/************************************************************************/

	//Count the number of rows first
	long wgfarows = -1;
	if (wgfastream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		while(nj != -1)
		{
			wgfastream.get(otemp,eof_col);
			wgfastream.ignore(INT_MAX,'\n');
			if( otemp[0] == chCommentChianti ) continue;
			nj = atoi(otemp);
			wgfarows++;
		}
		wgfastream.seekg(0,ios::beg);
	}
	else 
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chTraFilename);


	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\n\nTransition Probability File: %s\n",chTraFilename);
		fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Wavelength(A)|Ein A\n");
	}


	//line_index_col is the length(+1) of the level indexes in the WGFA file
	const int line_index_col = 6;
	//line_nrg_to_eina is the # of columns to skip from wavelength to eina in WGFA file
	const int line_nrg_to_eina =15;
	//line_eina_col is the length(+1) of einsteinA in WGFA
	const int line_eina_col = 16;
	char lvltemp[line_index_col];
	//Start reading WGFA file
	if (wgfastream.is_open())
	{
		for (long ii = 0;ii<wgfarows;ii++)
		{
			wgfastream.get(lvltemp,line_index_col);
			if( lvltemp[0] == chCommentChianti )
			{
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			long ipLoInFile = atoi(lvltemp);
			wgfastream.get(lvltemp,line_index_col);
			long ipHiInFile = atoi(lvltemp);

			// ipLo and ipHi will be manipulated below, want to retain original vals for prints
			long ipLo = ipLoInFile - 1;
			long ipHi = ipHiInFile - 1;

			if( atmdat.lgChiantiExp )
			{
				/* If either upper or lower index is -1 in the relational vector,
				 * skip that line in the wgfa file.
				 * Otherwise translate the level indices.*/
				if( intExperIndex[ipLo] == -1 || intExperIndex[ipHi] == -1 )
				{
					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}
				else
				{
					ipHi = intExperIndex.at(ipHi);
					if (ipHi < long(indexold2new.size()))
					{
						ipHi = indexold2new[ipHi];
					}
					else
					{
						ipHi = -1;
					}
					ipLo = intExperIndex.at(ipLo);
					if (ipLo < long(indexold2new.size()))
					{
						ipLo = indexold2new[ipLo];
					}
					else
					{
						ipLo = -1;
					}
				}
			}
			else
			{
				long testlo = -1, testhi = -1;

				try
				{
					testlo = indexold2new[ipLo];
					testhi = indexold2new[ipHi];
				}
				catch ( out_of_range& /* e */ )
				{
					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"NOTE: An out of range exception has occurred"
								" reading in data from %s\n",chTraFilename);
						fprintf(ioQQQ," The line in the file containing the unidentifiable"
								" levels has been ignored.\n");
						fprintf(ioQQQ,"There is no reason for alarm."
								" This message is just for documentation.\n");
					}

					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}

				if(  testlo == -1 || testhi == -1 )
				{
					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}
				else
				{
					ipLo = testlo;
					ipHi = testhi;
				}
			}

			if( ipLo >= nMolLevs || ipHi >= nMolLevs )
			{
				// skip these lines
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}
	
			if( ipHi == ipLo )
			{
				fprintf( ioQQQ," WARNING: Upper level = lower for a radiative transition in %s. Ignoring.\n", chTraFilename );
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"<%s>\t%li:%li\t%li:%li\t",dBaseSpecies[intNS].chLabel,ipLoInFile,ipHiInFile,ipLo+1,ipHi+1);
			}
	
			ASSERT( ipHi != ipLo );
			ASSERT( ipHi >= 0 );
			ASSERT( ipLo >= 0 );

			// sometimes the CHIANTI datafiles list the highest index first as in the middle of these five lines in ne_10.wgfa:
			//    ...
			//    8   10       187.5299      0.000e+00      4.127e+05                 3d 2D1.5 -                  4s 2S0.5           E2
			//    9   10       187.6573      0.000e+00      6.197e+05                 3d 2D2.5 -                  4s 2S0.5           E2
			//   11   10   4842624.0000      1.499e-05      9.423e-06                 4p 2P0.5 -                  4s 2S0.5           E1
			//    1   11         9.7085      1.892e-02      6.695e+11                 1s 2S0.5 -                  4p 2P0.5           E1
			//    2   11        48.5157      6.787e-02      9.618e+10                 2s 2S0.5 -                  4p 2P0.5           E1
			//    ...
			// so, just set ipHi (ipLo) equal to the max (min) of the two indices.
			// NB NB NB it looks like this may depend upon whether one uses observed or theoretical energies.

			//Read in wavelengths
			char trantemp[lvl_nrg_col];
			wgfastream.get(trantemp,lvl_nrg_col);
			fWLAng = (realnum)atof(trantemp);
			if( DEBUGSTATE && atmdat.lgChiantiExp)
			{
				fprintf(ioQQQ,"%.4f\t",fWLAng);
			}

			/* \todo 2 CHIANTI labels the H 1 2-photon transition as z wavelength of zero.
			 * Should we just ignore all of the wavelengths in this file and use the
			 * difference of level energies instead. */

			if( ipHi < ipLo )
			{
				long swap = ipHi;
				ipHi = ipLo;
				ipLo = swap;
			}

			/* If the given wavelength is negative, then theoretical energies are being used.
			 * Take the difference in stored theoretical energies.
			 * It should equal the absolute value of the wavelength in the wgfa file. */
			if( fWLAng <= 0. ) // && !atmdat.lgChiantiExp )
			{
				//if( fWLAng < 0.)
					//fprintf( ioQQQ," WARNING: Negative wavelength for species %6s, indices %3li %3li \n", dBaseSpecies[intNS].chLabel, ipLo, ipHi);
				fWLAng = (realnum)(1e8/abs(dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN()));
			}

			if( DEBUGSTATE && !atmdat.lgChiantiExp)
			{
				fprintf(ioQQQ,"%.4f\t",fWLAng);
			}
			//Skip from end of Wavelength to Einstein A and read in
			wgfastream.seekg(line_nrg_to_eina,ios::cur);
			wgfastream.get(trantemp,line_eina_col);
			feinsteina = (realnum)atof(trantemp);
			if( feinsteina == 0. )
			{
				static bool notPrintedYet = true;
				if( notPrintedYet && atmdat.lgChiantiPrint)
				{
					fprintf( ioQQQ," CAUTION: Radiative rate(s) equal to zero in %s.\n", chTraFilename );
					notPrintedYet = false;
				}
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}
			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"%.3e\n",feinsteina);
			}

			fixit("may need to do something with these rates");
			//Read in the rest of the line and look for auto
			wgfastream.getline(chLine,INT_MAX);
			TransitionList::iterator tr = dBaseTrans[intNS].begin()+ipdBaseTrans[intNS][ipHi][ipLo];
			if( nMatch("auto", chLine) )
			{
				if( (*tr).hasEmis() )
				{
					(*tr).Emis().AutoIonizFrac() =
						feinsteina/((*tr).Emis().Aul() + feinsteina);
					ASSERT( (*tr).Emis().AutoIonizFrac() >= 0. );
					ASSERT( (*tr).Emis().AutoIonizFrac() <= 1. );
				}
				continue;
			}

			if( (*tr).hasEmis() )
			{
				fprintf(ioQQQ," PROBLEM duplicate transition found by atmdat_chianti in %s, "
						"wavelength=%f\n", chTraFilename,fWLAng);
				wgfastream.close();
				cdEXIT(EXIT_FAILURE);
			}
			(*tr).AddLine2Stack();
			(*tr).Emis().Aul() = feinsteina;

			fenergyWN = (realnum)(1e+8/fWLAng);

			// TODO::Check the wavelength in the file with the difference in energy levels

			(*tr).EnergyWN() = fenergyWN;
			(*tr).WLAng() = (realnum)(1e+8/fenergyWN/RefIndex(fenergyWN));
			(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());

			(*tr).setComment( db_comment_tran_levels( ipLoInFile, ipHiInFile ) );
		}
	}
	else fprintf( ioQQQ, " The data file %s is corrupted .\n",chTraFilename);
	wgfastream.close();

	/* Malloc space for splines */
	AtmolCollSplines[intNS].reserve(nMolLevs);
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		AtmolCollSplines[intNS].reserve(ipHi,nMolLevs);
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			AtmolCollSplines[intNS].reserve(ipHi,ipLo,ipNCOLLIDER);
		}
	}
	AtmolCollSplines[intNS].alloc();

	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			for( long k=0; k<ipNCOLLIDER; k++ )
			{
				/* initialize all spline variables */
				AtmolCollSplines[intNS][ipHi][ipLo][k].collspline = NULL;
				AtmolCollSplines[intNS][ipHi][ipLo][k].SplineSecDer = NULL;
				AtmolCollSplines[intNS][ipHi][ipLo][k].nSplinePts = -1; 
				AtmolCollSplines[intNS][ipHi][ipLo][k].intTranType = -1;
				AtmolCollSplines[intNS][ipHi][ipLo][k].EnergyDiff = BIGDOUBLE;
				AtmolCollSplines[intNS][ipHi][ipLo][k].ScalingParam = BIGDOUBLE;
			}
		}
	}

	/************************************/
	/*** Read in the collisional data ***/
	/************************************/

	// ipCollider 0 is electrons, 1 is protons
	for( long ipCollider=0; ipCollider<=1; ipCollider++ )
	{
		char chFilename[FILENAME_PATH_LENGTH_2];

		if( ipCollider == ipELECTRON )
		{
			strcpy( chFilename, chEleColFilename );
		}
		else if( ipCollider == ipPROTON )
		{
			if( !lgProtonData )
				break;
			strcpy( chFilename, chProColFilename );
		}
		else
			TotalInsanity();

		/*Dummy string used for convenience*/
		strcpy(chLine,"A");

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"\n\nCollision Data File: %s\n",chTraFilename);
			fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Spline Points\n");
		}

		fstream splupsstream;
		open_data(splupsstream, chFilename,mode_r);

		//cs_eof_col is the length(+1) of the first column used for finding the end of file
		const int cs_eof_col = 4;
		//cs_index_col is the length(+1) of the indexes in the CS file
		const int cs_index_col = 4;
		//cs_trantype_col is the length(+1) of the transition type in the CS file
		const int cs_trantype_col = 4;
		//cs_values_col is the length(+1) of the other values in the CS file
		//including: GF, nrg diff, scaling parameter, and spline points
		const int cs_values_col = 11;
		//Determine the number of rows in the CS file
		if (splupsstream.is_open())
		{
			int nj = 0;
			//splupslines is -1 since the loop runs 1 extra time
			long splupslines = -1;
			char otemp[cs_eof_col];
			while(nj != -1)
			{
				splupsstream.get(otemp,cs_eof_col);
				splupsstream.ignore(INT_MAX,'\n');
				nj = atoi(otemp);
				splupslines++;
			}
			splupsstream.seekg(0,ios::beg);
	
			for (int m = 0;m<splupslines;m++)
			{
				if( ipCollider == ipELECTRON )
				{
					splupsstream.seekg(6,ios::cur);
				}

				/* level indices */
				splupsstream.get(otemp,cs_index_col);
				long ipLo = atoi(otemp)-1;
				splupsstream.get(otemp,cs_index_col);
				long ipHi = atoi(otemp)-1;

				long ipLoFile = ipLo;
				long ipHiFile = ipHi;

				/* If either upper or lower index is -1 in the relational vector,
				* skip that line in the splups file.
				* Otherwise translate the level indices.*/
				if( atmdat.lgChiantiExp )
				{
					if( intExperIndex[ipLo] == - 1 || intExperIndex[ipHi] == -1 )
					{
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}
					else
					{
						ipHi = intExperIndex.at(ipHi);
						if (ipHi < long(indexold2new.size()))
						{
							ipHi = indexold2new[ipHi];
						}
						else
						{
							ipHi = -1;
						}
						ipLo = intExperIndex.at(ipLo);
						if (ipLo < long(indexold2new.size()))
						{
							ipLo = indexold2new[ipLo];
						}
						else
						{
							ipLo = -1;
						}
					}
				}
				else
				{
					long testlo = -1, testhi = -1;

					/* With level trimming on it is possible that there can be rows that
					 * have to be skipped when using theoretical
					 * since the levels no longer exist */
					try
					{
						testlo = indexold2new[ipLo];
						testhi = indexold2new[ipHi];
					}
					catch ( out_of_range& /* e */ )
					{
						if( DEBUGSTATE )
						{
							fprintf(ioQQQ,"NOTE: An out of range exception has occurred"
									" reading in data from %s\n",chEleColFilename);
							fprintf(ioQQQ," The line in the file containing the unidentifiable"
									" levels has been ignored.\n");
							fprintf(ioQQQ,"There is no reason for alarm."
									" This message is for documentation.\n");
						}
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}

					if( testlo == -1 || testhi == -1 )
					{
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}
					else
					{
						ipLo = testlo;
						ipHi = testhi;
					}
				}

				if( ipLo >= nMolLevs || ipHi >= nMolLevs )
				{
					// skip these transitions
					splupsstream.ignore(INT_MAX,'\n');
					continue;
				}

				if( ipHi < ipLo )
				{
					long swap = ipHi;
					ipHi = ipLo;
					ipLo = swap;
				}

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"<%s>\t%li:%li\t%li:%li",dBaseSpecies[intNS].chLabel,ipLoFile+1,ipHiFile+1,ipLo+1,ipHi+1);
				}

				/*Transition Type*/
				splupsstream.get(otemp,cs_trantype_col);
				intTranType = atoi(otemp);
				char qtemp[cs_values_col];
				splupsstream.get(qtemp,cs_values_col);
				/*Energy Difference*/
				splupsstream.get(qtemp,cs_values_col);
				fEnergyDiff = atof(qtemp);
				/*Scaling Parameter*/
				splupsstream.get(qtemp,cs_values_col);
				fScalingParam = atof(qtemp);

				ASSERT( ipLo != ipHi );
				ASSERT( ipLo >= 0 && ipLo < nMolLevs );
				ASSERT( ipHi >= 0 && ipHi < nMolLevs );
				ASSERT( AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline == NULL );
				ASSERT( AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer == NULL );

				const int CHIANTI_SPLINE_MAX=9, CHIANTI_SPLINE_MIN=5;
				STATIC_ASSERT(CHIANTI_SPLINE_MAX > CHIANTI_SPLINE_MIN);

				/*We malloc the space here*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline =
					(double *)MALLOC((unsigned long)(CHIANTI_SPLINE_MAX)*sizeof(double));
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer =
					(double *)MALLOC((unsigned long)(CHIANTI_SPLINE_MAX)*sizeof(double));

				/* always read at least CHIANTI_SPLINE_MIN */
				for( intsplinepts=0; intsplinepts<=CHIANTI_SPLINE_MAX; intsplinepts++ )
				{
					//Look at the next character to see if it is the end of line.
					char p = splupsstream.peek();
					if( p == '\n' )
					{
						break;
					}
					else
					{
						if( intsplinepts >= CHIANTI_SPLINE_MAX )
						{
							fprintf( ioQQQ, " WARNING: More spline points than expected in %s, indices %3li %3li.  Ignoring extras.\n", chFilename, ipHi, ipLo );
							break;
						}
						ASSERT( intsplinepts < CHIANTI_SPLINE_MAX );
						double temp;
						//Store a single spline point then look for more
						splupsstream.get(qtemp,cs_values_col);
						temp = atof(qtemp);
						if( DEBUGSTATE )
						{
							fprintf(ioQQQ,"\t%.3e",temp);
						}
						// intTranType == 6 means log10 of numbers have been fit => allow negative numbers
						// intTranType < 6 means linear numbers have been fit => negative numbers are unphysical
						if( intTranType < 6 )
							temp = max( temp, 0. );
						AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intsplinepts] = temp;
					}
				}

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"\n");
				}

				ASSERT( intsplinepts > 2 );

				/*The zeroth element contains the number of spline points*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].nSplinePts = intsplinepts;
				/*Transition type*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].intTranType = intTranType;
				/*Energy difference between two levels in Rydbergs*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].EnergyDiff = fEnergyDiff;
				/*Scaling parameter C*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].ScalingParam = fScalingParam;

				/*Once the spline points have been filled,fill the second derivatives structure*/
				/*Creating spline points array*/
				vector<double> xs (intsplinepts),
					spl(intsplinepts),
					spl2(intsplinepts);

				for(intxs=0;intxs<intsplinepts;intxs++)
				{
					double coeff = (double)1/(intsplinepts-1);
					xs[intxs] = coeff*intxs;
					spl[intxs] = AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intxs];
				}

				spline(&xs[0], &spl[0],intsplinepts,2e31,2e31,&spl2[0]);

				/*Filling the second derivative structure*/
				for( long ii=0; ii<intsplinepts; ii++)
				{
					AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer[ii] = spl2[ii];
				}

				splupsstream.ignore(INT_MAX,'\n');
			}
			splupsstream.close();
		}
	}

	// close open file handles
	fclose( ioElecCollData );
	if( lgProtonData )
		fclose( ioProtCollData );

	//Chianti had bad theo level data so we used experimental
	//Changing lgChiantiExp back to false so next speices will use theoretical
	if( lgChiaBadTheo )
	{
		atmdat.lgChiantiExp = false;
	}

	return;
}
