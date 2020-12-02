/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*grid_do called by cdDrive, this returns 0 if things went ok, 1 for disaster */
#include "cddefines.h"
#include "conv.h"
#include "input.h"
#include "called.h"
#include "version.h"
#include "init.h"
#include "prt.h"
#include "trace.h"
#include "grainvar.h"
#include "parse.h"
#include "grid.h"
#include "atmdat.h"
#include "flux.h"

/* grid_do called by cdDrive, calls gridXspec or lgOptimize_do, returns false if ok, true for disaster */
bool grid_do()
{
	char chLine[INPUT_LINE_LENGTH], 
	  chNote[8];
	long int i, 
	  ii, 
	  j;
	realnum ptem[LIMPAR]; 

	DEBUG_ENTRY( "grid_do()" );

	/* main driver for optimization runs
	 * Drives cloudy to grid variables;*/

	/* code originally written by R.F. Carswell, IOA Cambridge */

	/* this will be number of times grid calls cloudy */
	optimize.nOptimiz = 0;

	/* variables with optimizer */
	for( i=0; i < LIMPAR; i++ )
	{
		optimize.OptIncrm[i] = 0.;
		optimize.varang[i][0] = -FLT_MAX;
		optimize.varang[i][1] = FLT_MAX;
		/* this should be overwritten by format of vary line */
		strcpy( optimize.chVarFmt[i], "error - no optimizer line image was set" );
	}

	/* necessary to do this to keep all lines in */
	prt.lgFaintOn = false;
	conv.LimFail = 1000;

	/* this initializes variables at the start of each simulation
	* in a grid, before the parser is called - this must set any values
	* that may be changed by the command parser */
	InitDefaultsPreparse();

	/* Read the isotope data and allocate the required space */
	LoadIsotopes();

	optimize.lgInitialParse = true;

	/* call READR the first time to scan off all variable options */
	/* this is just an initial parsing to get the number of iterations and
	 * the number of varied parameters.  The other Init* routines are not 
	 * called after this because this is all done again later for each grid point */
	ParseCommands();

	optimize.lgInitialParse = false;

	/* >>chng 00 aug 09, return memory allocated for grains, they are not used, PvH */
	gv.clear();

	optimize.nvary = optimize.nparm;

	/* option to change default increments; if zero then leave as is */
	for( i=0; i < LIMPAR; i++ )
	{
		if( optimize.OptIncrm[i] != 0. )
		{
			optimize.vincr[i] = optimize.OptIncrm[i];
		}
	}

	if( called.lgTalk )
	{
		/* check that at least 1 observed quantity was entered */
		unsigned long nObsQuant = optimize.xLineInt_Obs.size() + optimize.ContNFnu.size() +
			optimize.temp_obs.size() + optimize.ColDen_Obs.size();
		if( optimize.lgOptLum )
			nObsQuant++;
		if( optimize.lgOptDiam )
			nObsQuant++;
		if( nObsQuant == 0 && !grid.lgGrid )
		{
			fprintf( ioQQQ, " The input stream has vary commands, but\n" );
			fprintf( ioQQQ, " no observed quantities were entered.  Whats up?\n" );
			fprintf( ioQQQ, " Use the NO VARY command if you intended to disable optimization.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that the total number of parameters to vary is greater than 1 */
		if( optimize.nvary < 1 )
		{
			fprintf( ioQQQ, " No parameters to vary were entered. Whats up?\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( optimize.nvary > long(nObsQuant) && !grid.lgGrid )
		{
			fprintf( ioQQQ, " PROBLEM - More parameters are varied then there are observables.\n" );
			fprintf( ioQQQ, " PROBLEM - This run may not converge as a result.\n" );
			fprintf( ioQQQ, " PROBLEM - Please reduce the number of free parameters,"
				 " or add more observables.\n" );
		}

		if( strcmp(optimize.chOptRtn,"XSPE") == 0 && optimize.nRangeSet != optimize.nvary )
		{
			fprintf( ioQQQ, " Every parameter with a VARY option must have a GRID specified,\n" );
			fprintf( ioQQQ, " and the GRID must be specified after the VARY option.\n" );
			fprintf( ioQQQ, " These requirements were not satisfied for %ld parameter(s).\n",
				 abs(optimize.nvary - optimize.nRangeSet) );
			cdEXIT(EXIT_FAILURE);
		}

		/* lgTrOptm set with trace grid command */
		if( trace.lgTrOptm )
		{
			for( i=0; i < optimize.nvary; i++ )
			{
				/*print the command format as debugging aid */
				fprintf( ioQQQ, "%s\n", optimize.chVarFmt[i]);

				/* now generate the actual command with parameter,
				 * there will be from 1 to 3 numbers on the line */
				if( optimize.nvarxt[i] == 1 )
				{
					/* case with 1 parameter */
					sprintf( chLine, optimize.chVarFmt[i], optimize.vparm[0][i] );
				}

				else if( optimize.nvarxt[i] == 2 )
				{
					/* case with 2 parameter */
					sprintf( chLine, optimize.chVarFmt[i], optimize.vparm[0][i],
						 optimize.vparm[1][i]);
				}

				else if( optimize.nvarxt[i] == 3 )
				{
					/* case with 3 parameter */
					sprintf( chLine, optimize.chVarFmt[i], 
						 optimize.vparm[0][i], optimize.vparm[1][i], optimize.vparm[2][i] );
				}

				else if( optimize.nvarxt[i] == 4 )
				{
					/* case with 4 parameter */
					sprintf( chLine, optimize.chVarFmt[i], 
						 optimize.vparm[0][i], optimize.vparm[1][i], optimize.vparm[2][i],
						 optimize.vparm[3][i] );
				}

				else if( optimize.nvarxt[i] == 5 )
				{
					/* case with 5 parameter */
					sprintf( chLine, optimize.chVarFmt[i], 
						 optimize.vparm[0][i], optimize.vparm[1][i], optimize.vparm[2][i],
						 optimize.vparm[3][i], optimize.vparm[4][i]);
				}

				else
				{
					fprintf(ioQQQ,"The number of variable options on this line makes no sense to me1\n");
					cdEXIT(EXIT_FAILURE);
				}

				/* print the resulting command line*/
				fprintf( ioQQQ, "%s\n", chLine );
			}
		}

		/* say who we are */
		if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
			fprintf( ioQQQ, "%58cGrid  Driver\n", ' ' );
		else
			fprintf( ioQQQ, "%54cOptimization  Driver\n", ' ' );
		int indent = (int)((122 - strlen(t_version::Inst().chVersion))/2);
		fprintf( ioQQQ, "%*cCloudy %s\n\n",indent,' ',t_version::Inst().chVersion);
		fprintf( ioQQQ, "%23c**************************************%7.7s**************************************\n",
			 ' ', t_version::Inst().chDate );
		fprintf( ioQQQ, "%23c*%81c*\n", ' ', ' ' );

		/* now echo initial input quantities with flag for vary */
		/* first loop steps over all command lines entered */
		for( i=0; i <= input.nSave; i++ )
		{
			/* put space to start line, overwrite if vary found */
			strcpy( chNote, "       " );
			/* loop over all vary commands, see if this is one */
			for( j=0; j < optimize.nvary; j++ )
			{
				/* input.nSave is on C array counting, rest are on fortran */
				if( i == optimize.nvfpnt[j] )
				{
					/* this is a vary command, put keyword at start */
					strcpy( chNote, "VARY>>>" );
				}
			}

			// exclude lines from init files
			if( input.InclLevel[i] == 0 )
				fprintf( ioQQQ, "%22.7s * %-80s*\n", chNote, input.chCardSav[i] );
		}
		fprintf( ioQQQ, "%23c*%81c*\n", ' ', ' ' );
		fprintf( ioQQQ, "%23c***********************************************************************************\n\n\n", ' ' );

		/* option to trace logical flow within this sub */
		if( optimize.lgOptimFlow )
		{
			for( j=0; j < optimize.nvary; j++ )
			{
				i = optimize.nvfpnt[j];
				fprintf( ioQQQ, " trace:%80.80s\n", input.chCardSav[i]);
				fprintf( ioQQQ, "%80.80s\n", optimize.chVarFmt[j]);
				fprintf( ioQQQ, " number of variables on line:%4ld\n", 
					 optimize.nvarxt[j] );
				fprintf( ioQQQ, " Values:" );
				for( ii=1; ii <= optimize.nvarxt[j]; ii++ )
				{
					fprintf( ioQQQ, "%10.2e", optimize.vparm[ii-1][j] );
				}
				fprintf( ioQQQ, "\n" );
			}
		}

		if( strcmp(optimize.chOptRtn,"PHYM") == 0 )
		{
			fprintf( ioQQQ, " Up to %ld iterations will be performed,\n", 
				 optimize.nIterOptim );
			fprintf( ioQQQ, " and the final version of the input file will be written to the file %s\n", 
				 chOptimFileName );

			fprintf( ioQQQ, " The Phymir method will be used" );
			if( optimize.lgParallel )
			{
				if( cpu.i().lgMPI() )
					fprintf( ioQQQ, " in MPI mode.\n" );
				else
					fprintf( ioQQQ, " in parallel mode.\n" );

				fprintf( ioQQQ, " The maximum no. of CPU's to be used is %ld.\n",
					 optimize.useCPU );
			}
			else
				fprintf( ioQQQ, " in sequential mode.\n" );
		}

		else if( strcmp(optimize.chOptRtn,"SUBP") == 0 )
		{
			fprintf( ioQQQ, " Up to %ld iterations will be performed,\n", 
				 optimize.nIterOptim );
			fprintf( ioQQQ, " and the final version of the input file will be written to the file %s\n", 
				 chOptimFileName );

			fprintf( ioQQQ, " The Subplex method will be used.\n" );
		}

		else if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
		{
			fprintf( ioQQQ, " Producing grid output.\n" );
		}

		else
		{
			fprintf( ioQQQ, " I do not understand what method to use.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		fprintf( ioQQQ, "\n   %ld parameter(s) will be varied.  The first lines, and the increments are:\n", 
			 optimize.nvary );

		for( i=0; i < optimize.nvary; i++ )
		{
			optimize.varmax[i] = -FLT_MAX;
			optimize.varmin[i] = FLT_MAX;
			/*  write formatted to output using the format held in chVarFmt(np) */

			/* now generate the actual command with parameter,
			 * there will be from 1 to 3 numbers on the line */
			if( optimize.nvarxt[i] == 1 )
			{
				/* case with 1 parameter */
				sprintf( chLine, optimize.chVarFmt[i], optimize.vparm[0][i] );
			}

			else if( optimize.nvarxt[i] == 2 )
			{
				/* case with 2 parameter */
				sprintf( chLine, optimize.chVarFmt[i], optimize.vparm[0][i], optimize.vparm[1][i]);
			}

			else if( optimize.nvarxt[i] == 3 )
			{
				/* case with 3 parameter */
				sprintf( chLine, optimize.chVarFmt[i], 
					 optimize.vparm[0][i], optimize.vparm[1][i], optimize.vparm[2][i] );
			}

			else if( optimize.nvarxt[i] == 4 )
			{
				/* case with 4 parameter */
				sprintf( chLine, optimize.chVarFmt[i], 
					 optimize.vparm[0][i], optimize.vparm[1][i], optimize.vparm[2][i],
					 optimize.vparm[3][i] );
			}

			else if( optimize.nvarxt[i] == 5 )
			{
				/* case with 5 parameter */
				sprintf( chLine, optimize.chVarFmt[i], 
					 optimize.vparm[0][i], optimize.vparm[1][i], optimize.vparm[2][i],
					 optimize.vparm[3][i], optimize.vparm[4][i]);
			}

			else
			{
				fprintf(ioQQQ,"The number of variable options on this line makes no sense to me2\n");
				cdEXIT(EXIT_FAILURE);
			}

			fprintf( ioQQQ, "\n %s\n", chLine );
			if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
			{
				if( grid.paramValuesFromList[i].size() != 0U )
				{
					fprintf( ioQQQ, " %li parameters read from list.\n", grid.numParamValues[i] );
				}				
				else
				{
					fprintf( ioQQQ, " %s increment is %.3g, the limits are %.3g to %.3g\n",
						 grid.lgLinearIncrements[i] ? "Linear" : "Log",
						 grid.paramIncrements[i], grid.paramLimits[i][0], grid.paramLimits[i][1] );
				}
			}
			else
				fprintf( ioQQQ, " Initial increment is %.3g, the limits are %.3g to %.3g\n", 
					 optimize.vincr[i], optimize.varang[i][0], optimize.varang[i][1] );
		}
	}

	if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
	{
		if( called.lgTalk )
		{
			if( cpu.i().lgMPI() )
				fprintf( ioQQQ, "\n Running in MPI grid mode on %ld CPUs. ", cpu.i().nCPU() );
			else if( grid.lgParallel )
				fprintf( ioQQQ, "\n Running in parallel grid mode on %d CPUs. ", grid.useCPU );
			else
				fprintf( ioQQQ, "\n Running in single-CPU grid mode. " );
			fprintf( ioQQQ, "I will now start to write the input files.\n\n" );
		}

		for( j=0; j < optimize.nvary; j++ )
			ptem[j] = grid.paramLimits[j][0]; 
		for( j=optimize.nvary; j < LIMPAR; j++ )
		{
			ptem[j] = 0.f; 
			grid.paramIncrements[j] = 0.f;
			grid.lgLinearIncrements[j] = false;
		}

		gridXspec(ptem,optimize.nvary);

		if( called.lgTalk )
		{
			fprintf( ioQQQ, " **************************************************\n" );
			fprintf( ioQQQ, " **************************************************\n" );
			fprintf( ioQQQ, " **************************************************\n" );
			fprintf( ioQQQ, "\n Writing input files has been completed.\n\n\n" );
		}
	}
	else
	{
		called.lgTalk = false;
		/* this flag is needed to turn print on to have effect */
		called.lgTalkIsOK = false;

		lgAbort = lgOptimize_do();
	}

	return lgAbort;
}
