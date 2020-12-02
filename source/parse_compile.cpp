/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCompile compile Werner or kurucz model atmospheres into cloudy format, originally by K Volk,
 * also compile opacity and grains  */
#include "cddefines.h"
#include "continuum.h"
#include "dense.h"
#include "iso.h"
#include "grains.h"
#include "rfield.h"
#include "stars.h"
#include "parser.h"

void ParseCompile(Parser &p)
{
	long int ncell;
	char 
	  chRFI[FILENAME_PATH_LENGTH_2],
	  chSZD[FILENAME_PATH_LENGTH_2];

	DEBUG_ENTRY( "ParseCompile()" );

	/* >>chng 01 aug 24, remove compile opacity command */
	/* this option to compile opacities into file for later use */
	if( p.nMatch("OPAC") )
	{
		fprintf( ioQQQ, "The COMPILE OPACITIES command is currently not supported\n" );
		cdEXIT(EXIT_FAILURE);

#	if 0
		/* calls fill to set up continuum energy mesh if first call, 
		 * otherwise reset to original mesh */
		ContCreateMesh();

		/* read in some external data files, but only if this is first call */
		atmdat_readin();

		/* first generate the frequency array */
		ContCreatePointers();

		Badnell_rec_init();

		/* say that we want to compile the opacities */
		opac.lgCompileOpac = true;

		/* generate initial set of opacities but only if this is the first call 
		 * in this coreload */
		OpacityCreateAll();

		fprintf(ioQQQ, 
			"Success!!  Created file opacity.opc\nMake sure this is on the path.\n" );
		cdEXIT(EXIT_SUCCESS);
#	endif
	}

	/* >>chng 00 apr 27, modified for arbitrary file names by PvH
	 *
	 * this option to compile grains into file for later use
	 *
	 * the command supports the following syntax:
	 *
	 * COMPILE GRAINS
	 *     compile a standard set of opacity files
	 *
	 * COMPILE GRAINS <refr-ind-file> <size-distr-file> [ <no-bins> ]
	 *     compile a single opacity file
	 *
	 * Remarks:
	 * - the parameters of this command can be supplied in arbitrary order.
	 * - the file names can either be supplied as names between quotes or
	 *   as keywords; it is allowed to use a filename between quotes for
	 *   one file and a keyword for the other file; both names have to be
	 *   present in either form, there are no defaults.
	 * - filenames are recognized by their extension: .rfi or .mix for
	 *   refractive index files, and .szd for size distribution files,
	 *   this allows their sequence to be arbitrary.
	 * - the number-of-bins parameter is optional, it is defaulted to 10.
	 *
	 * NB NB NB NB NB NB NB NB NB NB NB NB NB
	 *
	 * - in order not to upset FFmtRead for reading the number-of-bins
	 *   parameter, all file names and keywords should be read first and
	 *   erased after being read ! to assure that all digits are erased,
	 *   the keywords 0M010, 0M100 and 1M000 are matched on all 5 characters.
	 *   if keywords are known not to contain digits or minus signs, erasing
	 *   is not necessary of course....
	 */
	if( p.nMatch("GRAI") )
	{

		/* calls fill to set up continuum energy mesh if first call, 
		 * otherwise reset to original mesh */
		ContCreateMesh();
		/* >>chng 06 dec 13, this had been followed by calls to atmdat_readin &
		 * ConCreatePointer which had problems because the code was not
		 * fully initialized yet.  Compile stars and compile grains would
		 * fail on mallocing an array of length zero */

		chRFI[0] = '\0';
		chSZD[0] = '\0';

		/* get first filename (either .rfi or .szd file) */
		if( p.nMatch( "\"" ) )
		{
			string chRead;
			if (p.GetQuote(chRead ))
				p.StringError();
			if( strstr_s(chRead.c_str(),".rfi") != NULL || strstr_s(chRead.c_str(),".mix") != NULL ) 
			{
				strcpy(chRFI,chRead.c_str());
			}
			else if( strstr_s(chRead.c_str(),".szd") != NULL ) 
			{
				strcpy(chSZD,chRead.c_str());
			}
			else 
			{
				fprintf( ioQQQ, " filename %s has unknown extension, sorry\n" , chRead.c_str() );
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* get second filename (either .rfi or .szd file) */
		if( p.nMatch( "\"" ) )
		{
			string chRead;
			if (p.GetQuote(chRead ))
				p.StringError();
			if( strstr_s(chRead.c_str(),".rfi") != NULL || strstr_s(chRead.c_str(),".mix") != NULL ) 
			{
				strcpy(chRFI,chRead.c_str());
			}
			else if( strstr_s(chRead.c_str(),".szd") != NULL ) 
			{
				strcpy(chSZD,chRead.c_str());
			}
			else 
			{
				fprintf( ioQQQ, " filename %s has unknown extension, sorry\n" , chRead.c_str() );
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* if no .rfi file was supplied between quotes, check for keywords */
		if( chRFI[0] == '\0' ) 
		{
			/* check on index of refraction names */
			if( p.nMatchErase("AC1-") )
			{
				/* amorphous carbon from Rouleau & Martin 1991 */
				strcpy(chRFI , "ac1-amcarb.rfi" );
				/* erase this keyword, it upsets FFmtRead */
			}
			else if( p.nMatchErase("BE1-"))
			{
				/* amorphous carbon from Rouleau & Martin 1991 */
				strcpy(chRFI , "be1-amcarb.rfi" );
				/* erase this keyword, it upsets FFmtRead */				
			}
			else if( p.nMatch( "GRAP") )
			{
				/* graphite */
				strcpy(chRFI , "graphite.rfi" );
			}
			else if( p.nMatch( "SILI" ) )
			{
				/* astronomical silicate */
				strcpy(chRFI , "silicate.rfi" );
			}
			else if( p.nMatch( " SIC" ) )
			{
				/* alpha-SiC */
				strcpy(chRFI , "sic.rfi" );
			}
			else if( p.nMatch( " PAH" ) )
			{
				/* original PAHs */
				strcpy(chRFI , "pah1.rfi" );
			}
			else if( p.nMatch( "GREY" ) || p.nMatch( "GRAY" ))
			{
				strcpy(chRFI , "grey.rfi" );
			}
		}

		/* if no .szd file was supplied between quotes, check for keywords */
		if( chSZD[0] == '\0' ) 
		{
			/* check on size distribution */
			if( p.nMatchErase("0M010") )
			{
				strcpy(chSZD , "0m010.szd" );
			}
			else if( p.nMatchErase("0M100") )
			{
				strcpy(chSZD , "0m100.szd" );
			}
			else if( p.nMatchErase("1M000") )
			{
				strcpy(chSZD , "1m000.szd" );
			}
			else if( p.nMatch( "ORIO" ) )
			{
				strcpy(chSZD , "orion.szd" );
			}
			else if( p.nMatch( " ISM" ) )
			{
				strcpy(chSZD , "ism.szd" );
			}
			else if( p.nMatchErase("AB08") )
			{
				/* Abel et al., 2008 size distribution */
				strcpy(chSZD , "ab08.szd" );
			}
			else if( p.nMatchErase("C15") )
			{
				/* small PAH, 15 C atoms */
				strcpy(chSZD , "c15.szd" );
			}
			else if( p.nMatchErase("C120") )
			{
				/* large PAH, 120 C atoms */
				strcpy(chSZD , "c120.szd" );
			}
		}

		/* the user has to supply either both the .rfi and .szd files, or neither
		 * (to compile the complete standard set of files); anything else is illegal */
		if( chRFI[0] == '\0' && chSZD[0] != '\0' )
		{
			fprintf(ioQQQ,"Sorry, but I did not recognize a refractive index file.\n");
			fprintf(ioQQQ,"Supply a file name between quotes or one of the following ");
			fprintf(ioQQQ,"keywords: ac1-amcarb, be1-amcarb, graphite, silicate, grey, pah\n");
			cdEXIT(EXIT_FAILURE);
		}

		if( chSZD[0] == '\0' && chRFI[0] != '\0' )
		{
			fprintf(ioQQQ,"Sorry, but I did not recognize a size distribution file.\n");
			fprintf(ioQQQ,"Supply a file name between quotes or one of the following ");
			fprintf(ioQQQ,"keywords: 0m010, 0m100, 1m000, ism, orion, c15, c120, ab08\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* compile the complete standard set of files */
		if( chRFI[0] == '\0' && chSZD[0] == '\0' )
		{
			/* ism graphite, single bin */
			mie_write_opc( "graphite.rfi" , "ism.szd" , 1 );

			/* ism silicate, single bin */
			mie_write_opc( "silicate.rfi" , "ism.szd" , 1 );

			/* ism graphite, 10 bins */
			mie_write_opc( "graphite.rfi" , "ism.szd" , 10 );

			/* ism silicate, 10 bins */
			mie_write_opc( "silicate.rfi" , "ism.szd" , 10 );

			/* orion graphite, single bin */
			mie_write_opc( "graphite.rfi" , "orion.szd" , 1 );

			/* orion silicate, single bin */
			mie_write_opc( "silicate.rfi" , "orion.szd" , 1 );

			/* orion graphite, 10 bins */
			mie_write_opc( "graphite.rfi" , "orion.szd" , 10 );

			/* orion silicate, 10 bins */
			mie_write_opc( "silicate.rfi" , "orion.szd" , 10 );

			/* 0.01 micron silicate */
			mie_write_opc( "silicate.rfi" , "0m010.szd" , 1 );

			/* 0.1 micron silicate */
			mie_write_opc( "silicate.rfi" , "0m100.szd" , 1 );

			/* 1 micron silicate */
			mie_write_opc( "silicate.rfi" , "1m000.szd" , 1 );

			/* 0.01 micron graphite */
			mie_write_opc( "graphite.rfi" , "0m010.szd" , 1 );

			/* 0.1 micron graphite */
			mie_write_opc( "graphite.rfi" , "0m100.szd" , 1 );

			/* 1 micron graphite */
			mie_write_opc( "graphite.rfi" , "1m000.szd" , 1 );

			/* grey single bin */
			mie_write_opc( "grey.rfi" , "ism.szd" , 1 );

			/* grey resolved distribution */
			mie_write_opc( "grey.rfi" , "ism.szd" , 10 );

			/* small pah */
			mie_write_opc( "pah1.rfi" , "c15.szd" , 1 );

			/* large pah */
			mie_write_opc( "pah1.rfi" , "c120.szd" , 1 );

			/* distributed pah */
			mie_write_opc( "pah1.rfi" , "ab08.szd" , 10 );

			/* single pah */
			mie_write_opc( "pah1.rfi" , "ab08.szd" , 1 );
		}
		/* this option is to compile a single type of grain */
		else
		{
			ncell = (long)p.FFmtRead();
			if( p.lgEOL() )
			{
				/* the default, 10 cells */
				ncell = 10;
			}
			if( ncell <= 0 )
			{
				fprintf(ioQQQ,"Number of bins must be positive.  Sorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
			/* this actually does the work */
			mie_write_opc( chRFI , chSZD , ncell );
		}

		fprintf(ioQQQ, 
			"Success!!  Created grain opacity file(s).\nMake sure this directory is on the path.\n" );
		cdEXIT(EXIT_SUCCESS);
	}

	/* compile recombination coefficients command */
	else if( p.nMatch("RECO") && p.nMatch("COEF") )
	{
		long ipISO;
		int nelem;

		if( p.nMatch("H-LI") )
			ipISO = ipH_LIKE;
		else if( p.nMatch("HE-L") )
			ipISO = ipHE_LIKE;
		else
		{
			fprintf(ioQQQ,"Sorry, but I did not recognize an iso sequence.\n");
			fprintf(ioQQQ,"The available options are H-like and He-like.\nSorry.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* compile he-like command - compiles table of recombination coeficients */
		iso_ctrl.lgCompileRecomb[ipISO] = true;

		/* we will want to create the rec coefficient for a large number of levels, then stop.
		 * this sets the number of levels to a large number.  macro is in helike.h */
		for( nelem = ipISO; nelem < LIMELM; nelem++)
		{
			long maxN;
			dense.lgElmtOn[nelem] = true;
			iso_sp[ipISO][nelem].nCollapsed_max = 0;

			if( nelem == ipISO )
				maxN = RREC_MAXN;
			else
				maxN = LIKE_RREC_MAXN( nelem );

			iso_sp[ipISO][nelem].n_HighestResolved_max = maxN;

			iso_update_num_levels( ipISO, nelem );
		}
	}

	else if( p.nMatch("STAR") )
	{
		bool lgProblems = false;

		/* check that frequency mesh is set up */
		ASSERT( rfield.lgMeshSetUp() );

		if( p.nMatch( "\"" ))
		{
			/* this is branch for for user-supplied *.ascii file */

			/* this will both scan in whatever label is inside the quotes in OrgCard, 
			 * but also remove the contents there and in chCard,
			 * so that following keywords will not trigger off it */
			string chRead;
			if (p.GetQuote( chRead ))
				p.StringError();

			string::size_type ptr = chRead.rfind('.');
			if( ptr != string::npos )
			{
				if( chRead.substr(ptr) == ".ascii" )
				{
					lgProblems = GridCompile( chRead.c_str() );
				}
				else if( chRead.substr(ptr) == ".stb99" )
				{
					// keyword to only include the stellar component
					// default is to include both the stellar and nebular component
					sb_mode mode;
					if( p.nMatch( "STEL" ) )
						mode = SB_STELLAR;
					else if( p.nMatch( "NEBU" ) )
						mode = SB_NEBULAR;
					else
						mode = SB_TOTAL;
					string chASC(chRead.c_str(), ptr);
					chASC += ".ascii";
					fprintf( ioQQQ, " Creating %s....\n", chASC.c_str() );
					lgProblems = StarburstInitialize( chRead.c_str(), chASC.c_str(), mode );
					lgProblems = lgProblems || GridCompile( chASC.c_str() );
				}
				else
				{
					fprintf( ioQQQ, " I did not recognize this file extension: %s\n",
						 chRead.substr(ptr).c_str() );
					lgProblems = true;
				}
			}
			else
			{
				fprintf( ioQQQ, " I did not find any file extension: %s\n", chRead.c_str() );
				lgProblems = true;
			}
		}
		else
		{
			/* this branch is intended to convert ascii versions of stellar
			 * atmosphere grids into a direct access version for faster access.
			 * the original file is usually named *.ascii, and the new direct
			 * access file will always be named *.mod.
			 * - if the *.ascii file does not exist, the grid will be skipped
			 * - if the *.ascii file exists, but the *.mod file does not, or is
			 *   out of date, a new *.mod file will be generated
			 * - if the *.mod file is up to date, it will not be touched. */

			process_counter pc;

			/* These are the current Atlas grids */
			lgProblems = lgProblems || AtlasCompile(pc);
			/* do the costar OB stars */
			lgProblems = lgProblems || CoStarCompile(pc);
			/* legacy Atlas grid - for backward compatibility only */
			lgProblems = lgProblems || Kurucz79Compile(pc);
			/* Mihalas grid - for backward compatibility only */
			lgProblems = lgProblems || MihalasCompile(pc);
			/* do the rauch PN central stars */
			lgProblems = lgProblems || RauchCompile(pc);
			/* do the Starburst99 sample output */
			lgProblems = lgProblems || StarburstCompile(pc);
			/* do the Tlusty OSTAR2002 grid */
			lgProblems = lgProblems || TlustyCompile(pc);
			/* do the Werner PN central stars - for backward compatibility only */
			lgProblems = lgProblems || WernerCompile(pc);
			/* WMBASIC O-star grid by Pauldrach */
			lgProblems = lgProblems || WMBASICCompile(pc);

			if( pc.nFound == 0 )
			{
				fprintf( ioQQQ, "\n PROBLEM - No ascii files were found!\n" );
				fprintf( ioQQQ, " Did you change directory to where the stellar atmosphere files are?\n" );
				fprintf( ioQQQ, " This command will only work on files in the local directory. Sorry.\n" );
				lgProblems = true;
			}
			else
			{
				fprintf( ioQQQ, "\n %d ascii file(s) found", pc.nFound );
				if( pc.notProcessed > 0 )
					fprintf( ioQQQ, ", %d file(s) up to date", pc.notProcessed );
				if( pc.nOK > 0 )
					fprintf( ioQQQ, ", %d update(s) OK", pc.nOK );
				if( pc.nFail > 0 )
					fprintf( ioQQQ, ", %d update(s) failed", pc.nFail );
				int nSkip = pc.nFound - pc.notProcessed - pc.nOK - pc.nFail;
				if( nSkip > 0 )
					fprintf( ioQQQ, ", %d file(s) skipped after failure", nSkip );
				fprintf( ioQQQ, ".\n" );
			}
		}

		if( lgProblems )
		{
			fprintf( ioQQQ, "\n Problems occurred during the compilation - check output.\n" );
		}
		else
		{
			fprintf( ioQQQ, "\n The compilation was successful!\n" );
			fprintf( ioQQQ, "\n Good Luck!!\n\n\n" );
		}

		cdEXIT( lgProblems ? ES_FAILURE : ES_SUCCESS );
	}
	else
	{
		fprintf( ioQQQ, " One of the keywords, GRAINS, RECO COEF, or STARS, must appear.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}
