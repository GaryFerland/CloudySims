/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParsePlot parse the plot command */
/*ParsePlotRangeOption parse range from plot command */
#include "cddefines.h"
#include "plot.h"
#include "rfield.h"
#include "parser.h"

/*ParsePlotRangeOption parse range from plot command, located below */
STATIC void ParsePlotRangeOption(Parser &p);

/*ParsePlotRangeContin parse range option for continuum on map command, located below */
STATIC void ParsePlotRangeContin(Parser &p);

void ParsePlot(Parser &p )
{

	DEBUG_ENTRY( "ParsePlot()" );

	/* total number of plots so far */
	plotCom.nplot += 1;

	/* plots are turned on */
	plotCom.lgPlotON = true;

	/* make sure we have not hit the limit, the dimension of variables in plot.h */
	if( plotCom.nplot > NDPLOT )
	{
		fprintf( ioQQQ, 
			" Too many plots; the limit is%3ld This one ignored.\n", 
		  NDPLOT );
		plotCom.nplot = NDPLOT;
	}

	if( p.nMatch(" MAP") )
	{
		/*  option to make plot of heating - cooling map */
		strcpy( plotCom.chPType[plotCom.nplot-1], " MAP" );
	}

	else if( p.nMatch("CONT") )
	{
		/* option to make "raw" plot, in internal units */
		if( p.nMatch(" RAW") )
		{
			strcpy( plotCom.chPType[plotCom.nplot-1], "CRAW" );
			/* option to make diffuse continuum plot */
		}
		else if( p.nMatch("DIFF") )
		{
			strcpy( plotCom.chPType[plotCom.nplot-1], "DIFF" );
			/* this is emitted continuum */
		}
		else if( p.nMatch("EMIT") )
		{
			strcpy( plotCom.chPType[plotCom.nplot-1], "EMIT" );
			/* this is outward and attenuated continuum */
		}
		else if( p.nMatch("OUTW") )
		{
			strcpy( plotCom.chPType[plotCom.nplot-1], "OUTW" );
			/* this is reflected continuum */
		}
		else if( p.nMatch("REFL") )
		{
			strcpy( plotCom.chPType[plotCom.nplot-1], "REFL" );
			/* this is continuum in photons */
		}
		else if( p.nMatch("PHOT") )
		{
			strcpy( plotCom.chPType[plotCom.nplot-1], "CPHT" );
		}
		else
		{
			strcpy( plotCom.chPType[plotCom.nplot-1], "CONT" );
		}
	}

	else if( p.nMatch("OPAC") )
	{
		if( p.nMatch("ABSO") )
		{
			/* plot absorption opacity */
			strcpy( plotCom.chPType[plotCom.nplot-1], "OPAA" );
		}
		else if( p.nMatch("SCAT") )
		{
			/* plot scattering opacity */
			strcpy( plotCom.chPType[plotCom.nplot-1], "OPAS" );
		}
		else if( p.nMatch("TOTA") )
		{
			/* plot total opacity */
			strcpy( plotCom.chPType[plotCom.nplot-1], "OPAT" );
		}
		else
		{
			/* plot total opacity for default */
			strcpy( plotCom.chPType[plotCom.nplot-1], "OPAT" );
		}
	}
	else
	{
		fprintf( ioQQQ, " The second keyword on the PLOT command must be CONTINUUM, _MAP, or OPACITY.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* option to turn on trace of plot logic */
	if( p.nMatch("TRAC") )
	{
		plotCom.lgPltTrace[plotCom.nplot-1] = true;
	}
	else
	{
		plotCom.lgPltTrace[plotCom.nplot-1] = false;
	}

	/* option to set min and max x-axis */
	if( strcmp(plotCom.chPType[plotCom.nplot-1]," MAP") == 0 )
	{
		/* this will be map of cooling and heating vs temp */
		ParsePlotRangeOption(p);
	}
	else
	{
		/* continuum map */
		ParsePlotRangeContin(p);
	}

	return;
}

/*ParsePlotRangeOption parse range from plot command */
STATIC void ParsePlotRangeOption(Parser &p )
{
	bool  
	  lgLogOn;
	realnum a;

	DEBUG_ENTRY( "ParsePlotRangeOption()" );

	/* pltxmn is min for x axis of plot */
	plotCom.pltxmn[plotCom.nplot-1] = (realnum)p.FFmtRead();

	if( !p.nMatch("RANG") )
	{
		/* no lines were enterd, so use default temperature limits */
		plotCom.pltxmn[plotCom.nplot-1] = 1.;
		plotCom.pltxmx[plotCom.nplot-1] = 9.;
		lgLogOn = true;
	}
	else if( p.lgEOL() )
	{
		/* no lines were enterd, so use default temperature limits */
		plotCom.pltxmn[plotCom.nplot-1] = 1.;
		plotCom.pltxmx[plotCom.nplot-1] = 9.;
		lgLogOn = true;
	}
	else
	{
		/* number entered, now interprete it */
		if( plotCom.pltxmn[plotCom.nplot-1] <= 10. )
		{
			lgLogOn = true;
		}
		else
		{
			lgLogOn = false;
		}
		/* linear option for temperature */
		if( p.nMatch("LINE") )
			lgLogOn = false;
		/* lower temp was entered, now how about upper temp */
		plotCom.pltxmx[plotCom.nplot-1] = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			if( lgLogOn )
			{
				plotCom.pltxmx[plotCom.nplot-1] = 9.;
			}
			else
			{
				plotCom.pltxmx[plotCom.nplot-1] = 1e9;
			}
		}
		else
		{
			/*  second number was entered, check for sanity */
			if( plotCom.pltxmx[plotCom.nplot-1] <= plotCom.pltxmn[plotCom.nplot-1] )
			{
				fprintf( ioQQQ, 
					" The second (maximum) temperature for the map is%10.2e, but this is less than the first (minimum) temperature,%10.2e\n", 
				  plotCom.pltxmx[plotCom.nplot-1], plotCom.pltxmn[plotCom.nplot-1] );
				fprintf( ioQQQ, " HELP!  I am confused!!\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* now force to be log */
	if( !lgLogOn )
	{
		if( plotCom.pltxmx[plotCom.nplot-1] <= 0. || plotCom.pltxmn[plotCom.nplot-1] <= 
		  0. )
		{
			fprintf( ioQQQ, 
				" Limits for temperature are negative.  This is impossible.   Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		plotCom.pltxmx[plotCom.nplot-1] = (realnum)log10(plotCom.pltxmx[plotCom.nplot-1]);
		plotCom.pltxmn[plotCom.nplot-1] = (realnum)log10(plotCom.pltxmn[plotCom.nplot-1]);
	}

	/* check that min is less than max */
	if( plotCom.pltxmn[plotCom.nplot-1] == plotCom.pltxmx[plotCom.nplot-1] )
	{
		fprintf( ioQQQ, " Upper and lower plot boundaries are equal.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	else if( plotCom.pltxmn[plotCom.nplot-1] > plotCom.pltxmx[plotCom.nplot-1] )
	{
		a = plotCom.pltxmx[plotCom.nplot-1];
		plotCom.pltxmx[plotCom.nplot-1] = plotCom.pltxmn[plotCom.nplot-1];
		plotCom.pltxmn[plotCom.nplot-1] = a;
	}

	return;
}

/*ParsePlotRangeContin set range for map to parse range option on map command */
STATIC void ParsePlotRangeContin(Parser &p )
{
	bool 
	  lgLogOn;
	realnum a;

	DEBUG_ENTRY( "ParsePlotRangeContin()" );

	plotCom.pltxmn[plotCom.nplot-1] = (realnum)p.FFmtRead();

	if( !p.nMatch("RANG") )
	{
		/* no numbers on line, just set default limits to energy array */
		plotCom.pltxmn[plotCom.nplot-1] = (realnum)log10(rfield.emm());
		plotCom.pltxmx[plotCom.nplot-1] = (realnum)log10(rfield.egamry());
		lgLogOn = true;
	}
	else if( p.lgEOL() )
	{
		/* no numbers on line, just set default limits to energy array */
		plotCom.pltxmn[plotCom.nplot-1] = (realnum)log10(rfield.emm());
		plotCom.pltxmx[plotCom.nplot-1] = (realnum)log10(rfield.egamry());
		lgLogOn = true;
	}
	else
	{
		/* lower limit entered, now interprete it */
		if( plotCom.pltxmn[plotCom.nplot-1] < 0. )
		{
			lgLogOn = true;
		}
		else if( plotCom.pltxmn[plotCom.nplot-1] == 0.0 )
		{
			/* option for first number to be zero, lower edge of array */
			lgLogOn = false;
			plotCom.pltxmn[plotCom.nplot-1] = rfield.emm();
		}
		else
		{
			/* a positive number was entered */
			lgLogOn = false;
		}

		/* now look at upper limit */
		plotCom.pltxmx[plotCom.nplot-1] = (realnum)p.FFmtRead();
		if( p.lgEOL() || plotCom.pltxmx[plotCom.nplot-1] == 0.0 )
		{
			/* option for second number to be omitted or zero */
			if( lgLogOn )
			{
				plotCom.pltxmx[plotCom.nplot-1] = (realnum)log10(rfield.egamry());
			}
			else
			{
				plotCom.pltxmx[plotCom.nplot-1] = rfield.egamry();
			}
		}
	}

	/* now force to be log */
	if( !lgLogOn )
	{
		plotCom.pltxmx[plotCom.nplot-1] = (realnum)log10(plotCom.pltxmx[plotCom.nplot-1]);
		plotCom.pltxmn[plotCom.nplot-1] = (realnum)log10(plotCom.pltxmn[plotCom.nplot-1]);
	}

	/* check that min is less than max */
	if( plotCom.pltxmn[plotCom.nplot-1] == plotCom.pltxmx[plotCom.nplot-1] )
	{
		fprintf( ioQQQ, " Upper and lower plot boundaries are equal.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	else if( plotCom.pltxmn[plotCom.nplot-1] > plotCom.pltxmx[plotCom.nplot-1] )
	{
		a = plotCom.pltxmx[plotCom.nplot-1];
		plotCom.pltxmx[plotCom.nplot-1] = plotCom.pltxmn[plotCom.nplot-1];
		plotCom.pltxmn[plotCom.nplot-1] = a;
	}

	return;
}
