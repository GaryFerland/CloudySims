/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*save_line parse save lines command, or actually do the save lines output */
/*Save_Line_RT parse the save line rt command - read in a set of lines */
#include "cddefines.h"
#include "cddrive.h"
#include "radius.h"
#include "opacity.h"
#include "phycon.h"
#include "dense.h"
#include "lines.h"
#include "h2.h"
#include "prt.h"
#include "iso.h"
#include "parser.h"
#include "count_ptr.h"
#include "save.h"
/* this is the limit to the number of emission lines we can store */
#define	NPUNLM	200L

/* implement the save line xxx command.  cumulative, structure, and
 * emissivity all use same code base and variables, so only one can be used
 * at present */

class SaveLineList
{
public:
	string chPLab[NPUNLM];
	long int nLinesEntered;
	realnum wavelength[NPUNLM];
	long int ipLine[NPUNLM];
	bool lgRelativeIntensity, lgMustGetLines, lgMustPrintFirstTime;
	SaveLineList() : nLinesEntered(0), lgRelativeIntensity(false), lgMustGetLines(true),
						  lgMustPrintFirstTime(true) {}
};

static count_ptr<SaveLineList> linelist[LIMPUN];

void parse_save_line(Parser &p, 
		     /* true, return rel intensity, false, log of luminosity or intensity I */
		     bool lgLog3,
			  ostringstream& chHeader,
	        long ipPun)
{
	char chTemp[INPUT_LINE_LENGTH];

	// save return value of cdLine, 0 for success, -number of lines for fail
	long int i;

	DEBUG_ENTRY( "parse_save_line()" );

	/* very first time this routine is called, chDo is "READ" and we read
	 * in lines from the input stream.  The line labels and wavelengths
	 * are store locally, and output in later calls to this routine
	 * following is flag saying whether to do relative intensity or
	 * absolute emissivity */
	linelist[ipPun] = count_ptr<SaveLineList>(new SaveLineList);
	
	linelist[ipPun]->lgRelativeIntensity = lgLog3;
	
	/* number of lines we will save */
	linelist[ipPun]->nLinesEntered = 0;
	
	/* get the next line, and check for eof */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Hit EOF while reading line list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	while( !p.hasCommand("END") )
	{
		if( linelist[ipPun]->nLinesEntered >= NPUNLM )
		{
			fprintf( ioQQQ, 
						" Too many lines have been entered; the limit is %ld.  Increase variable NPUNLM in routine save_line.\n", 
						linelist[ipPun]->nLinesEntered );
			cdEXIT(EXIT_FAILURE);
		}

		LineID line = p.getLineID();
		if( !p.lgReachedEnd() )
		{
			fprintf( ioQQQ, "parse_save_line: found junk at end of input line:\n" );
			p.showLocation();
			cdEXIT(EXIT_FAILURE);
		}
		linelist[ipPun]->chPLab[linelist[ipPun]->nLinesEntered] = line.chLabel;
		linelist[ipPun]->wavelength[linelist[ipPun]->nLinesEntered] = line.wave;
		
		/* this is total number stored so far */
		++linelist[ipPun]->nLinesEntered;
		
		/* get next line and check for eof */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		
	}
	
	sncatf( chHeader, "#depth" );
	for( i=0; i < linelist[ipPun]->nLinesEntered; i++ )
	{
		sprt_wl( chTemp, linelist[ipPun]->wavelength[i] );
		sncatf( chHeader,
			"\t%s %s", linelist[ipPun]->chPLab[i].c_str(), chTemp ); 
	}
	sncatf( chHeader, "\n" );
}

void save_line(FILE * ioPUN, /* the file we will write to */
  const char *chDo, 
  // intrinsic or emergent line emission?
  bool lgEmergent,
  long ipPun
	)
{
	long int i;
	double a[NPUNLM], 
	  absint, 
	  relint;

	DEBUG_ENTRY( "save_line()" );

	/* it is possible that we will get here after an initial temperature
	 * too high abort, and the line arrays will not have been defined.
	 * do no lines in this case.  must still do save so that there
	 * is not a missing line in the grid save output */
	long nLinesNow = LineSave.nsum >0 ? linelist[ipPun]->nLinesEntered : 0;

	bool lgBadLine = false;
	if( nzone <= 1 && linelist[ipPun]->lgMustGetLines )
	{
		for( i=0; i < nLinesNow; i++ )
		{
			linelist[ipPun]->ipLine[i] = 
				LineSave.findline(linelist[ipPun]->chPLab[i].c_str(), linelist[ipPun]->wavelength[i]);
			if( linelist[ipPun]->ipLine[i] <= 0 )
			{
				// missed line - ignore if H2 line since large model may be off
				if( !h2.lgEnabled && strncmp( linelist[ipPun]->chPLab[i].c_str() , "H2  " , 4 )==0 )
				{
					if( linelist[ipPun]->lgMustPrintFirstTime )
					{
						/* it's an H2 line and H2 is not being done - ignore it */
						fprintf( ioQQQ,"\nPROBLEM Did not find an H2 line, the large model is not "
									"included, so I will ignore it.  Log intensity set to -30.\n" );
						fprintf( ioQQQ,"I will totally ignore any future missed H2 lines\n\n");
						linelist[ipPun]->lgMustPrintFirstTime = false;
					}
					/* flag saying to ignore this line */
					linelist[ipPun]->ipLine[i] = -2;
				}
				else
				{
					fprintf( ioQQQ, " save_line could not find line ");
					prt_line_err( ioQQQ, linelist[ipPun]->chPLab[i].c_str(), linelist[ipPun]->wavelength[i] );
					linelist[ipPun]->ipLine[i] = -1;
					lgBadLine = true;
				}
			}
		}
		linelist[ipPun]->lgMustGetLines = false;
		if( lgBadLine )
		{
			cdEXIT(EXIT_FAILURE);
		}
	}

	if( strcmp(chDo,"PUNS") == 0 )
	{
		/* save lines emissivity command */
		/* save lines structure command */

		for( i=0; i < nLinesNow; i++ )
		{
			/* this version of cdEmis uses index, does not search, do not call if line could not be found */
			/* test on case where we abort before first zone is done
			 * this happens in grid when temperature bounds of code
			 * are exceeded.  In this case return small float */
			if( lgAbort && nzone <=1 )
				a[i] = SMALLFLOAT;
			else if( linelist[ipPun]->ipLine[i]>0 )
				cdEmis_ip(linelist[ipPun]->ipLine[i],&a[i],lgEmergent);
			else
				a[i] = 1e-30f;
		}
	}

	else if( strcmp(chDo,"PUNC") == 0 )
	{
		/* save lines cumulative command */		
		for( i=0; i < nLinesNow; i++ )
		{
			if ( (lgAbort && nzone <=1) || linelist[ipPun]->ipLine[i]<=0 )
			{
				a[i] = 0.;
			}
			else
			{
				cdLine_ip(linelist[ipPun]->ipLine[i],&relint,&absint,lgEmergent);
				if( linelist[ipPun]->lgRelativeIntensity )
					/* relative intensity case */
					a[i] = relint;
				else
					/* emissivity or luminosity case */
					a[i] = absint;
			}
		}		
	}
	else if( strcmp(chDo,"PUNO") == 0 )
	{
		/* save lines optical depth some command */		
		for( i=0; i < nLinesNow; i++ )
		{
			if ( (lgAbort && nzone <=1) || linelist[ipPun]->ipLine[i]<=0 )
			{
				a[i] = 0.;
			}
			else
			{
				TransitionProxy tr = LineSave.lines[linelist[ipPun]->ipLine[i]].getTransition();
				if (tr.associated())
					a[i] = tr.Emis().TauIn()*SQRTPI;
				else
					a[i] = 0.;
			}
		}		
	}
	else
	{
		fprintf( ioQQQ, 
			" unrecognized key for save_line=%4.4s\n", 
		  chDo );
		cdEXIT(EXIT_FAILURE);
	}

	fprintf( ioPUN, "%.5e", radius.depth_mid_zone );
	
	for( i=0; i < nLinesNow; i++ )
	{
		fprintf( ioPUN, "\t%.4e", a[i] );
	}
	fprintf( ioPUN, "\n" );

	return;
}

#define LIMLINE 10
static long int line_RT_type[LIMLINE] = 
  {LONG_MIN , LONG_MIN ,LONG_MIN , LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN },
	line_RT_ipISO[LIMLINE] =  
  {LONG_MIN , LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN },
		line_RT_nelem[LIMLINE] =  
  {LONG_MIN , LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN },
			line_RT_ipHi[LIMLINE] =  
  {LONG_MIN , LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN },
				line_RT_ipLo[LIMLINE] = 
  {LONG_MIN , LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN };
static bool lgMustPrintHeader=true;
static long int nLine=-1;

/*Save_Line_RT parse the save line rt command - read in a set of lines */
void Parse_Save_Line_RT(Parser &p)
{
	/* save line rt parameters */
	DEBUG_ENTRY( "Parse_Save_Line_RT()" );

	/* very first time this routine is called, chDo is "READ" and we read
	 * in lines from the input stream.  The line labels and wavelengths
	 * are store locally, and output in later calls to this routine */
	
	/* say that we must print the header */
	lgMustPrintHeader = true;
	
	/* get the next line, and check for eof */
	nLine = 0;
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Hit EOF while reading line list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	do
	{
		if(nLine>=LIMLINE )
		{
			fprintf(ioQQQ," PUNCH RT has too many lines - increase LIMLINE in save_line.cpp\n");
			cdEXIT(EXIT_FAILURE);
		}
		
		/* right now it just does lines in the iso sequences */
		line_RT_type[nLine] = (long int)p.FFmtRead();
		line_RT_ipISO[nLine] = (long int)p.FFmtRead();
		line_RT_nelem[nLine] = (long int)p.FFmtRead();
		line_RT_ipHi[nLine] = (long int)p.FFmtRead();
		line_RT_ipLo[nLine] = (long int)p.FFmtRead();
		
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, 
						" there must be five numbers on this line\n");
			p.PrintLine(ioQQQ);
			cdEXIT(EXIT_FAILURE);
		}
		
		/* now increment number of lines */
		++nLine;
		
		/* now get next line until we hit eof or the word END */
		p.getline();
	} while( !p.m_lgEOF && !p.nMatch( "END") );
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Save_Line_RT hit end of file looking for END of RT lines\n");
		p.PrintLine(ioQQQ);
		cdEXIT(EXIT_FAILURE);
	}
}

void Save_Line_RT( 
	FILE * ioPUN )
{
	/* save line rt parameters */

	DEBUG_ENTRY( "Save_Line_RT()" );


	static char chLabel[LIMLINE][30];
	long int n;
	if( lgMustPrintHeader )
	{
		fprintf( ioPUN , "Line\tP(con,inc)\tAul\tgl\tgu\n");
		for( n=0; n<nLine; ++n )
		{
			TransitionProxy tr = iso_sp[line_RT_ipISO[n]][line_RT_nelem[n]].trans(line_RT_ipHi[n],line_RT_ipLo[n]);
			/* print info for header of file, line id and pump rate */
			sprintf( chLabel[n] , "%s ", 
					chLineLbl(tr).c_str() );
			fprintf( ioPUN , "%s ", chLabel[n] );
			fprintf( ioPUN , "%.4e ",
					tr.Emis().pump());
			fprintf( ioPUN , "%.4e ",
					tr.Emis().Aul());
			fprintf( ioPUN , "%.0f ",
					(*tr.Lo()).g());
			fprintf( ioPUN , "%.0f ",
					(*tr.Hi()).g());
			fprintf( ioPUN , "\n");
			
			if( line_RT_type[n]!=0. )
			{
				/* for now code only exists for H He like iso - assert this */
				fprintf( ioQQQ, 
							" Save_Line_RT only H, He like allowed for now\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		fprintf( ioPUN , "Line\tTauIn\tPopLo\tPopHi\tCul\tk(line)\tk(con,abs)\tk(con,scat)\n");
		lgMustPrintHeader = false;
	}
	
	fprintf(ioPUN, "RADIUS\t%e\tDEPTH\t%e\tTe\t%e\tNe\t%e\n",
			  radius.Radius_mid_zone ,
			  radius.depth_mid_zone ,
			  phycon.te ,
			  dense.eden );
	for( n=0; n<nLine; ++n )
	{
		TransitionProxy tr = iso_sp[line_RT_ipISO[n]][line_RT_nelem[n]].trans(line_RT_ipHi[n],line_RT_ipLo[n]);

		/* index for line within continuum array */
		long int ipCont = tr.ipCont();
		fprintf( ioPUN , "%s ", chLabel[n] );
		fprintf( ioPUN , "\t%e\t%e\t%e",
					tr.Emis().TauIn() ,
					(*tr.Lo()).Pop(),
					(*tr.Hi()).Pop()
			);
		fprintf( ioPUN , "\t%e",
					tr.Coll().ColUL( colliders ) / dense.EdenHCorr
			);
		
		fprintf( ioPUN , "\t%e\t%e\t%e\n",
					tr.Emis().PopOpc(),
					opac.opacity_abs[ipCont-1] ,
					opac.opacity_sct[ipCont-1]
			);
	}
}
 
#	undef LIMELM

