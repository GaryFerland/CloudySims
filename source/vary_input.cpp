/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*vary_input sets input lines to feed into cloudy in optimization runs */
#include "cddefines.h"
#include "input.h"
#include "save.h"
#include "grid.h"

void vary_input(bool *lgLimOK,
		int grid_index)
{
	long int i, 
	  np;

	DEBUG_ENTRY( "vary_input()" );

	// this would indicate int overflow, but is mainly there to keep the compiler from
	// whining about an unused variable...
	if( grid_index < -1 )
		TotalInsanity();

	/* set up chCardSav(n) array like Gary's input file, using input
	 * variable parameters p(i), and format information held in
	 * the common block /parmv/. Results written to common /kardsv/.
	 */

	/* will be set false if limit to a variable exceeded
	 * this is returned to calling code as problem indication*/
	*lgLimOK = true;

	if( cpu.i().lgMaster() || !grid.lgGrid )
		fprintf( ioQQQ, " **************************************************\n" );

	/* echo the variable input lines for this run */
	for( i=0; i < optimize.nvary; i++ )
	{
		bool lgLimitHit = false;

		np = optimize.nvfpnt[i];

		// check if the keyword _LOG is present; the optimizer may not work
		// correctly if it is not optimizing logarithmic quantities.
		//
		// exceptions are the commands ILLUMINATE and RATIO since they vary
		// quantities of order unity anyway, and the commands DLAW and FUDGE
		// since they are entirely defined by the user.
		//
		// it is ok not to convert to upper case first since the command line
		// image is completely under our own control.
		if( !optimize.lgOptimizeAsLinear[i] )
		{
			if( !nMatch( " LOG", optimize.chVarFmt[i] ) )
			{
				fprintf( ioQQQ, " vary_input: internal error - keyword _LOG not found!\n" );
				TotalInsanity();
			}
		}

		/* write formatted to the character string chCardSav(np),
		 * using the format held in chVarFmt(np) */

		/* >>chng 05 aug 09, by RP, both were == change to > and < */
		if( grid.paramIncrements[i] >= 0. &&
		    ( optimize.vparm[0][i] < optimize.varang[i][0] ||
		      optimize.vparm[0][i] > optimize.varang[i][1] ) )
		{
			*lgLimOK = false;
			lgLimitHit = true;
		}
		if( grid.paramIncrements[i] < 0. &&
		    ( optimize.vparm[0][i] > optimize.varang[i][0] ||
		      optimize.vparm[0][i] < optimize.varang[i][1] ) )
		{
			*lgLimOK = false;
			lgLimitHit = true;
		}

		/* now generate the actual command with parameter,
		 * there will be from 1 to 3 numbers on the line */
		if( optimize.nvarxt[i] == 1 )
		{
			/* case with 1 parameter */
			sprintf( input.chCardSav[np], optimize.chVarFmt[i], optimize.vparm[0][i] );
		}

		else if( optimize.nvarxt[i] == 2 )
		{
			/* case with 2 parameters */
			sprintf( input.chCardSav[np], optimize.chVarFmt[i], optimize.vparm[0][i],
				 optimize.vparm[1][i] );
		}

		else if( optimize.nvarxt[i] == 3 )
		{
			/* case with 3 parameters */
			sprintf( input.chCardSav[np], optimize.chVarFmt[i], optimize.vparm[0][i],
				 optimize.vparm[1][i], optimize.vparm[2][i] );
		}

		else if( optimize.nvarxt[i] == 4 )
		{
			/* case with 4 parameters */
			sprintf( input.chCardSav[np], optimize.chVarFmt[i], optimize.vparm[0][i],
				 optimize.vparm[1][i], optimize.vparm[2][i], optimize.vparm[3][i] );
		}

		else if( optimize.nvarxt[i] == 5 )
		{
			/* case with 5 parameters */
			sprintf( input.chCardSav[np], optimize.chVarFmt[i], 
				 optimize.vparm[0][i], optimize.vparm[1][i], optimize.vparm[2][i], 
				 optimize.vparm[3][i], optimize.vparm[4][i]);
		}

		else
		{
			fprintf(ioQQQ,"The number of variable options on this line makes no sense to me5\n");
			cdEXIT(EXIT_FAILURE);
		}

		strncpy( input.chCardStrip[np], input.chCardSav[np], INPUT_LINE_LENGTH );

		if( cpu.i().lgMaster() || !grid.lgGrid )
		{
			fprintf( ioQQQ, " %s\n", input.chCardSav[np] );
			if( lgLimitHit )
				fprintf( ioQQQ, " >>> Limit to variable exceeded.\n" );
		}
	}

	if( cpu.i().lgMaster() && grid.lgGrid )
	{
		// write the line images to an input script, one file for each grid point
		fstream io;
		string fnam = GridPointPrefix(grid_index) + save.chRedirectPrefix + ".in";
		open_data( io, fnam.c_str(), mode_w );
		// input.nSave has unusual definition, it is one less than the number of lines stored
		for( int i=0; i <= input.nSave; ++i )
			if( input.InclLevel[i] == 0 )
				io << input.chCardSav[i] << endl;
	}

	return;
}
