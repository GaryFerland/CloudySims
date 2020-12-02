/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtLineSum parse print line sum command to enter set of lines into sum  */
#include "cddefines.h"
#include "cddrive.h"
#include "radius.h"
#include "lines.h"
#include "parser.h"
/* this is the limit to the number of lines we can save */
#define	NRDSUM	300L
#include "prt.h"

static string strSMLab[NRDSUM];
static long int *ipLine;
static long nlsum;
static realnum *wavelength;

void ParsePrtLineSum(Parser &p)
{
	static bool lgFirst=true;

	/* remember whether we have been called before */

	DEBUG_ENTRY( "ParsePrtLineSum()" );

	/* >>chng 03 jan 23, if not first call, do not allocate space, 
	 * had aborted, which was bad in optized runs, or in a grid. 
	 * Bug caught by Melekh Bohdan */
	if( lgFirst )
	{
		/* do not malloc space again */
		lgFirst = false;
		wavelength = ((realnum *)MALLOC( sizeof(realnum )*NRDSUM ));
		
		/* create space for the array of array indices for lines*/
		ipLine = ((long int *)MALLOC(NRDSUM*sizeof(long)));
	}
	
	/* now read in lines */
	nlsum = 0;
	bool lgEND = false;
	while( !lgEND )
	{
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		
		if( ! p.hasCommand("END" ) )
		{
			if( nlsum >= NRDSUM )
			{
				fprintf( ioQQQ, 
							" Too many lines have been entered; the limit is %li.  Increase NRDSUM in PrtLineSum.\n", 
							NRDSUM );
				cdEXIT(EXIT_FAILURE);
			}
			
			LineID line = p.getLineID();
			if( !p.lgReachedEnd() )
			{
				fprintf( ioQQQ, "ParsePrtLineSum: found junk at end of input line:\n" );
				p.showLocation();
				cdEXIT(EXIT_FAILURE);
			}
			strSMLab[nlsum] = line.chLabel;
			wavelength[nlsum] = line.wave;
			++nlsum;
		}
		else
		{
			lgEND = true;
		}
	}
}
double PrtLineSum(void)
{
	long int i;

	/* remember whether we have been called before */

	double absint, 
	  relint ,
	  sum=-1.;

	DEBUG_ENTRY( "PrtLineSum()" );

	sum = 0.;
	/* this can be called during setup mode, in which case we do nothing */
	if( LineSave.ipass <= 0 )
	{ 
		return sum;
	}
	
	if( nzone == 1 )
	{
		bool	lgFail = false;
		for( i=0; i < nlsum; i++ )
		{
			/* save the array index for each line */
			if(( ipLine[i] = LineSave.findline( strSMLab[i].c_str(), wavelength[i]) ) <= 0 )
			{
				fprintf( ioQQQ, " PrtLineSum could not find line " );
				prt_line_err( ioQQQ, strSMLab[i].c_str(), wavelength[i] );
				lgFail = true;
			}
		}
		if( lgFail )
			cdEXIT(EXIT_FAILURE);
	}
	
	/* now sum the line */
	for( i=0; i < nlsum; i++ )
	{
		/* this version of chLine uses index, does not search*/
		cdLine_ip(ipLine[i],&relint,&absint);
		absint /= radius.Conv2PrtInten;
		sum += absint;
	}
	return sum;
}

