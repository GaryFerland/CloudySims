/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAbundances parse and read in composition as set by abundances command */
#include "cddefines.h"
#include "abund.h"
#include "dense.h"
#include "elementnames.h"
#include "input.h"
#include "parser.h"


STATIC FILE *open_abn_file( string chFile )
{
	FILE *ioDATA;

	// first try local directory, then data/abundances
	if( (ioDATA = open_data( chFile.c_str(), "r", AS_LOCAL_ONLY_TRY ) ) == NULL )
	{
		char chPath[FILENAME_PATH_LENGTH_2] = { 0 };	/*entire path to file including file name*/

		/* change chPath to the abundances/chFile, using the proper delimiter */
		strcpy( chPath, "abundances" );
		strcat( chPath, input.chDelimiter );
		strcat( chPath, chFile.c_str() );
		ioDATA = open_data( chPath, "r" );	// will abort if not found
	}

	return	ioDATA;
}



void ParseAbundances(Parser &p)
			/* following set true by grains command,
			 * so this will not set more grains if grains already on. */
{
	bool lgLog;
	long int i;
	double absav[LIMELM], 
	  chk;

	DEBUG_ENTRY( "ParseAbundances()" );

	/* abundances no longer solar */
	abund.lgAbnSolar = false;

	if( p.nMatch("STAR") )
	{
		/* Fred Hamann's star burst galaxy mixture -- includes a number which isn't an abundance */
		abund_starburst(p);
		return;
	}

	/* GetQuote should be above the other options so other commands don't trip
	 * if there is any overlapping text in the quotes */
	string chFile;	/*file name for table read */
	//records whether or not quotes were found and stores what is between them in chFile
	bool lgMatchFound = true;
	if( p.GetQuote( chFile ) )
		lgMatchFound = false;
	
	bool lgPrint = p.nMatch( "PRINT");
	bool lgIsotp = p.nMatch( "ISOT" );

	if(!lgMatchFound && !lgIsotp)
	{
		lgMatchFound = true;
		if(p.nMatch("CAME"))
			chFile = "Cameron.abn";
		else if (p.nMatch("CRAB"))
			chFile = "Crab.abn";
		else if (p.nMatch("GASS"))
			chFile = "solar_GASS10.abn";
		else if (p.nMatch("H II") || p.nMatch("HII ") || p.nMatch("ORIO"))
			chFile = "HII.abn";
		else if (p.nMatch("ISM"))
			chFile = "ISM.abn";
		else if (p.nMatch("NOVA"))
			chFile = "nova.abn";
		else if (p.nMatch(" AGB") || p.nMatch("AGB ") || p.nMatch("PLAN"))
			chFile = "PN.abn";
		else if (p.nMatch("PRIM"))
			chFile = "primordial.abn";
		else if (p.nMatch("OLD ") && p.nMatch("SOLA"))
			chFile = "solar84.abn";
		else if (p.nMatch("ALLE"))
			chFile = "allen73.abn";
		else
			lgMatchFound = false;
	}
	else if( lgIsotp && chFile.length() == 0 )
	{
		if( p.nMatch( "ASPL" ) )
			chFile = "Asplund09-iso.abn";
		else if( p.nMatch( "LODDERS03" ) )
			chFile = "Lodders03-iso.abn";
		else if( p.nMatch( "LODDERS09" ) )
			chFile = "Lodders09-iso.abn";
		else if( p.nMatch( "ROSM" ) )
			chFile = "Rosman98-iso.abn";
		else
		{
			fprintf(ioQQQ, "Unknown isotope abundances file:\t %s\n", chFile.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
	}

	if(lgMatchFound && !lgIsotp)
	{
		// we have a file name - are other parameters present?
		bool lgGrainsON = true;
		if( p.nMatch("NO GR") != 0 )
			lgGrainsON = false;

		bool lgQHeat = true;
		if( p.nMatch("NO QH") != 0 )
			lgQHeat = false;

		abund.lgAbundancesSet = true;

		// initialization
		for( int nelem=1; nelem < LIMELM; nelem++ )
		{
			/* turn off all elements except Hydrogen,
			 * then turn on each element when found in the file */
			char chDUMMY[INPUT_LINE_LENGTH];
			sprintf(chDUMMY,"element %s off ", elementnames.chElementName[nelem] );
			p.setline(chDUMMY);
			// need to retain this flag, set by user to explicitly turn off element,
			// overriding our default to turn on elements that appear in abundances list
			bool lgSave = dense.lgElmtSetOff[nelem];
			ParseElement( p );
			dense.lgElmtSetOff[nelem] = lgSave;
		}

		FILE *ioDATA = open_abn_file( chFile );

		/* Sets abundance of Hydrogen to 1.0 in the expectation
		 * that other values are abundances relative to Hydrogen,
		 * although Hydrogen's abundance will still be set to
		 * whatever value is found in the code. This is just in
		 * case it is assumed */
		abund.solar[ipHYDROGEN] = 1.0;//The value found in abund.solar[0] before this point is found to be 1.0 as well

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
				fprintf(ioQQQ, "PROBLEM in ABUNDANCES: Encountered unexpected empty line.\n");
				cdEXIT(EXIT_FAILURE);
			}

			// If keyword "grains" is found on a line calls grain command
			// NOQHEAT may have been on the line in the input deck, or in the abn file
			/* copy input to CAPS to make sure hide not on it */
			char chLineCAPS[INPUT_LINE_LENGTH];
			strcpy( chLineCAPS , chLine );
			caps( chLineCAPS );
			if( nMatch("GRAINS", chLineCAPS) != 0 )
			{
				if( lgPrint )
					fprintf(ioQQQ,"%s\n",chLine);
				//Makes sure grains have not already been set and are on
				//Either way it skips to the next loop iteration and does not check the line for elements
				if( !p.m_lgDSet && lgGrainsON )
				{
					if( !lgQHeat )
						strcat( chLineCAPS, " NO QHEAT");

					p.setline(chLineCAPS);
					ParseGrain(p);
					continue;
				}
				else
					continue;
			}

			i = 1;
			bool lgFound = false;
			for(int nelem=0; nelem<LIMELM; nelem++)
			{
				if( nMatch( elementnames.chElementNameShort[nelem], chLineCAPS ) != 0 )
				{
					lgFound = true;
					abund.solar[nelem] = FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
					if( abund.solar[nelem] <=0. )
					{
						fprintf(ioQQQ, "PROBLEM in ABUNDANCES: negative abundance not allowed.\n");
						fprintf(ioQQQ, "Non-positive abundance found on this line: %s\n", chLine);
						cdEXIT(EXIT_FAILURE);
					}

					// do not implicitly turn on an element that has explicitly been turned off
					if( dense.lgElmtSetOff[nelem] )
						break;

					/* turn on any element found while parsing the file */
					char chDUMMY[INPUT_LINE_LENGTH];
					sprintf(chDUMMY,"element %s on ", elementnames.chElementName[nelem] );
					p.setline(chDUMMY);
					ParseElement( p );

					//we shouldn't need to continue once an element name is found on the line...
					break;
				}
			}
			if( !lgFound )
			{
				fprintf(ioQQQ, "PROBLEM in ABUNDANCES: did not identify element name on this line: %s\n",
					chLine);
				cdEXIT(EXIT_FAILURE);
			}
		}

		//Normalizes all elements in abund.solar relative to the quantity of Hydrogen
		ASSERT( abund.solar[0]>0. );
		for(int nelem=0; nelem<LIMELM; nelem++)
		{
			abund.solar[nelem] /= abund.solar[ipHYDROGEN];
			if( lgPrint && dense.lgElmtOn[nelem] )
				fprintf(ioQQQ,"%s\t%.3e\t%.3f\n",elementnames.chElementName[nelem],
						abund.solar[nelem] , log10(SDIV(abund.solar[nelem])) );
		}
		fclose( ioDATA );
		return;
	}
	else if( !lgIsotp )
	{
		abund.lgAbundancesSet = true;

		/* Looks for a number on the parser*/
		absav[0] = p.FFmtRead();
		/* this branch at least one number on the line - an abundance */
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " PROBLEM in ABUNDANCES: I did not find a keyword, file name, or any numbers. Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
		absav[1] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* this branch, we found one, but not two, numbers -
			 * must be a few special case */
			if( p.nMatch(" ALL") )
			{
				/* special option, all abundances will be this number */
				if( absav[0] <= 0. )
				{
					absav[0] = exp10(absav[0]);
				}
				for( i=1; i < LIMELM; i++ )
				{
					abund.solar[i] = (realnum)absav[0];
				}

			}
			else
			{
				fprintf( ioQQQ, 
					" Only one number was found. Did you include the keyword ALL? Sorry.\n");
				cdEXIT(EXIT_FAILURE);
			}

			/* normal return */
			return;
		}

		/* we get here if there is a second number - read in all abundances */
		/* abund.npSolar is changed by ParseElement, which was called previously
		 * in this routine.  It should not be changed by anything in this loop.
		 * Make sure this is the case.
		 */
		long nElem = abund.npSolar;
		for( i=2; i < nElem; i++ )
		{
			absav[i] = p.FFmtRead();
			if( p.lgEOL() )
			{
				/* read CONTINUE line if J scanned off end before reading all abundances */
				do
				{
					p.getline();
					if( p.m_lgEOF )
					{
						fprintf( ioQQQ, " There MUST be%3ld abundances entered, there were only%3ld.  Sorry.\n",
									abund.npSolar, i );
						cdEXIT(EXIT_FAILURE);
					}
				} while( p.isComment() );

				p.echo();

				if( ! p.hasCommand("CONT") )
				{
					fprintf( ioQQQ, " There MUST be%3ld abundances entered, there were only%3ld.  Sorry.\n",
						abund.npSolar, i );
					cdEXIT(EXIT_FAILURE);
				}
				else
				{
					absav[i] = p.FFmtRead();
					if( p.lgEOL() )
					{
						fprintf( ioQQQ, " There MUST be%3ld abundances entered, there were only%3ld.  Sorry.\n", 
							abund.npSolar, i);
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
		}
		/* would only fail if npSolar is changed within this loop */
		ASSERT( nElem == abund.npSolar);

		/* fell through to here after reading in N abundances for N elements
		 * check that there are no more abundances on the line - that would
		 * be an error - a typo */
		chk = p.FFmtRead();
		if( !p.lgEOL() || (chk!=0.) )
		{
			/* got another number, not lgEOL, so too many numbers */
			fprintf( ioQQQ, " There were more than %ld abundances entered\n",
				abund.npSolar );
			fprintf( ioQQQ, " Could there have been a typo somewhere?\n" );
		}

		/* are numbers scale factors, or log of abund rel to H?? */
		lgLog = false;
		for( i=0; i < abund.npSolar; i++ )
			if( absav[i] < 0. )
				lgLog = true;

		if( lgLog )
		{
			/* entered as log of number rel to hydrogen */
			for( i=0; i < abund.npSolar; i++ )
				abund.solar[abund.ipSolar[i]-1] = (realnum)exp10(absav[i]);
		}
		else
		{
			/* scale factors relative to solar */
			for( i=0; i < abund.npSolar; i++ )
				abund.solar[abund.ipSolar[i]-1] *= (realnum)absav[i];
		}

		/* check whether the abundances are reasonable */
		for( i=1; i < LIMELM; i++ )
		{
			if( abund.solar[i] > 0.2 )
			{
				fprintf( ioQQQ, " Is an abundance of %.3e relative to H reasonable for %2.2s?\n",
					abund.solar[i], elementnames.chElementSym[i] );
			}
		}
		return;
	}
	else if( lgIsotp )
	{
		FILE *ioDATA = open_abn_file( chFile );

		char chLine[INPUT_LINE_LENGTH];			/*to store file lines*/
		bool lgEOL;
		int nRead[LIMELM] = { 0 };

		while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
		{
			if( chLine[0]=='*' )
				break;

			/* skip comment */
			if( chLine[0]=='#' )
				continue;
			if( chLine[0]=='\n' || chLine[0]=='\0' )
			{
				fprintf(ioQQQ, "PROBLEM in ABUNDANCES ISOTOPES: Encountered unexpected empty line.\n");
				cdEXIT(EXIT_FAILURE);
			}

			long i = 1;
			int ielem = (int)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL) - 1;
			int Aiso  = (int)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
			realnum Fiso = (realnum)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
			int j = abund.IsoAbn[ielem].setAbn( Aiso, Fiso );
			if( j == -1 )
			{
				fprintf(ioQQQ,
					"PROBLEM in ABUNDANCES ISOTOPES: Could not store isotope fraction (%7.4f) for ^%d %s\n",
					Fiso, Aiso, elementnames.chElementSym[i]);
					cdEXIT( EXIT_FAILURE );
			}
			++nRead[ielem];
		}
		fclose( ioDATA );

		/* Express all isotope fractions in terms of the most abundant isotope. */
		for( int i = 0; i < LIMELM; i++ )
		{
			if( nRead[i] != abund.IsoAbn[i].getNiso() )
			{
				fprintf(ioQQQ, "Abundaces Isotopes: %s requires %d isotope pairs to be specified, but found %d\n",
					elementnames.chElementName[i], abund.IsoAbn[i].getNiso(), nRead[i]);
				cdEXIT(EXIT_FAILURE);
			}
			if( abund.IsoAbn[i].isAnyIllegal() )
			{
				fprintf(ioQQQ, "Abundaces Isotopes: Non-positive isotope fractions are illegal.\n");
				fprintf(ioQQQ, "File: %s\t Read: %s\t iso:", chFile.c_str(), elementnames.chElementName[i]);
				abund.IsoAbn[i].prtIsoPairs( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}
			abund.IsoAbn[i].normAbn( );
			//	abund.IsoAbn[i].prtIsoPairs( stdout );
		}

		return;
	}
	//return;
	//Should never reach this return, because if & else branch return
}
