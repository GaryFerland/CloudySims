/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SetPrintLineCol set some parameters for the main line block & wl format */
/*prt_LineLabels save all labels and wavelengths for emission line array */
/*sprt_wl write wavelength to string - must be kept parallel with prt_wl */
/*prt_wl - print floating wavelength in Angstroms, in output format */
/* prt_line print the line label, followed by the line wl, and the wavelength of the closest match, if given */
/* prt_line_inlist print line suitable for output list, label not enclosed in quotation marks */
#include "cddefines.h"
#include "lines.h"
#include "prt.h"
#include "generic_state.h"

t_prt prt;
t_line_col prt_linecol;

void t_line_col::zero()
{
	DEBUG_ENTRY( "t_line_col::zero()" );
	/* set the format of the main output line list */
	absint_len = 8;
	relint_len = 9;
	col_gap_len = 6;
}


/*SetPrintLineCol set some parameters for the main line block & wl format */
void SetPrintLineCol ()
{
	LineSave.wl_length = (int) LineSave.sig_figs + 2;

	/* set the format of the main output line list */
	prt_linecol.column_len = (NCHLAB-1) + LineSave.wl_length + prt_linecol.absint_len + prt_linecol.relint_len + 3;

	prt_linecol.col_gap.assign( prt_linecol.col_gap_len, ' ' );
	prt_linecol.relint_outrange.assign( prt_linecol.relint_len , '*' );

	return;
}



/*prt_wl print floating wavelength in Angstroms, in output format */
void prt_wl( FILE *ioOUT , realnum wl )
{
	char chString[100];
	DEBUG_ENTRY( "prt_wl()" );

	sprt_wl( chString , wl );

	fprintf(ioOUT, "%.*s", LineSave.wl_length, chString );
	return;
}

/* write wavelength to string */
void sprt_wl( char *chString , realnum wl )
{
	char chUnits[10];

	DEBUG_ENTRY( "sprt_wl()" );

	/* print in A unless > 1e4, then use microns */
	if( wl > 1e8 )
	{
		/* centimeters */
		strcpy( chUnits , "c" );
		wl /= 1e8;
	}
	else if( wl > 1e4 )
	{
		/* microns */
		strcpy( chUnits , "m" );
		wl /= 1e4;
	}
	else if( wl == 0. )
	{
		strcpy( chUnits , " " );
	}
	else
	{
		/* Angstroms units */
		strcpy( chUnits , "A" );
	}

	/* want total of LineSave.sig_figs sig figs */
	if( wl==0. )
	{
		sprintf(chString, "%*i", LineSave.wl_length-1, 0 );
	}
	else
	{
		int n = LineSave.sig_figs - 1 - (int)log10(wl);
		if( n > 0 )
		{
			sprintf( chString, "%*.*f", (int)LineSave.sig_figs+1, n, wl );
		}
		else if (wl < INT_MAX)
		{
			sprintf(chString, "%*i", (int)LineSave.sig_figs+1, (int)wl );
		}
		else
		{
			sprintf(chString, "%*c", (int)LineSave.sig_figs+1, '*' );
		}
	}

	strcat( chString , chUnits );
	return;
}

/*prt_LineLabels save all labels and wavelengths for emission line array */
void prt_LineLabels(
	/* io file handle */
	FILE * ioOUT ,
	/* print all if true, if false then do not print parts of 
	 * transferred lines */
	bool lgPrintAll )
{
	long int i;

	DEBUG_ENTRY( "prt_LineLabels()" );

	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSave.lines[i].isSeparator() )
		{
			/*fprintf( ioOUT, "%s ", LineSave.lines[i].chALab );*/
			fprintf( ioOUT, "####\t%s",LineSave.chHoldComments[(int)LineSave.lines[i].wavelength()] ); 
		}
		else
		{
			if( !lgPrintAll &&
				 ( LineSave.lines[i].isInward() ||
					LineSave.lines[i].isCollisional() ||
					LineSave.lines[i].isPump() ||
					LineSave.lines[i].isHeat() )
				)
				/* option to do not print lots of redundant labels 
				 * lgPrintAll is false by default set true with LONG option
				 * on save line labels command */
				continue;
			/* this format chosen to be identical to that used by final */
			fprintf( ioOUT, "%li\t%s\t", 
						i,
						LineSave.lines[i].label().c_str() );
			/* skip over leading spaces - a formatting problem */
			long int j = 0;
			string comment = LineSave.lines[i].chComment();
			while( comment[j]!='\0' && comment[j]==' ')
				++j;
			/* comment entered when line intensity generated */
			fprintf( ioOUT , "\t# %s" , comment.substr(j).c_str() );
		}
		fprintf( ioOUT, "\n" );
	}
	return;
}

/* prt_line_err produce an error message containing the line label and wavelength,
 * 		followed, if given, by the wavelength of the closest line of the same label */
void prt_line_err ( FILE *ioOUT, const char *label, realnum wvlng )
{
	fprintf( ioOUT, "with label (between quotes) \"%s\" and wavelength ", label );
	prt_wl ( ioOUT, wvlng );
	fprintf( ioOUT, ".\n" );

	return;
}



/* prt_line_inlist print line suitable for output list, label not enclosed in quotation marks */
void prt_line_inlist ( FILE *ioOUT, const char *label, realnum wvlng )
{
	fprintf( ioOUT, "%-*s\t", NCHLAB-1, label );
	prt_wl ( ioOUT, wvlng );

	return;
}


void t_prt_matrix::zero()
{
	DEBUG_ENTRY( "t_prt_matrix::zero()" );

	for( int i = 0; i < CHARS_SPECIES; ++i )
	{
		species[ i ] = '\0';
		speciesLevels[ i ] = '\0';
	}
}

void t_prt_matrix::setSpecies( const string &sspec )
{
	DEBUG_ENTRY( "t_prt_matrix::setSpecies()" );

	strncpy( species, sspec.c_str(), CHARS_SPECIES-1 );
	species[CHARS_SPECIES-1] = '\0';
	char *p = strchr( species, '[' );
	char *p1 = p;
	// Erase levels specification, e.g., "[*:*]", if any
	if( p != NULL )
		*p = '\0';

	strncpy( speciesLevels, sspec.c_str(), CHARS_SPECIES-1 );
	speciesLevels[CHARS_SPECIES-1] = '\0';
	if( p1 == NULL )
	{
		if( strlen( speciesLevels ) + 3 < CHARS_SPECIES )
		{
			strcat( speciesLevels, "[:]" );
			fprintf( ioQQQ, "NOTE: Converted '%s' to '%s'\n",
					species, speciesLevels ); 
		}
		else
		{
			fprintf( ioQQQ, "PROBLEM: Cannot convert '%s' to '%s[:]'. Species label too wide.\n",
					species, speciesLevels ); 
			cdEXIT( EXIT_FAILURE );
		}
	}
}

void t_prt_matrix::resolveLevels()
{
	DEBUG_ENTRY( "t_prt_matrix::`resolveSpecies()" );

	if( strlen( speciesLevels ) == 0 )
		return;

	getLevelsGeneric( speciesLevels, true, speciesLevelList );
}

void t_prt_matrix::prtRates( const long nlevels_local, const multi_arr<double,2,C_TYPE> &a,
				valarray<double> &b )
{
	DEBUG_ENTRY( "t_prt_matrix::prtRates()" );

	if( speciesLevelList.size() == 0 )
		return;

	fprintf( ioQQQ, "'%s' lvl / creation /=>rates", species);
	for( vector<long>::iterator ipLo = speciesLevelList.begin();
		 ipLo != speciesLevelList.end(); ++ipLo )
	{
		if( ipLo == speciesLevelList.begin() )
			fprintf( ioQQQ, "\t%3ld", *ipLo+1 );
		else
			fprintf( ioQQQ, "\t%11ld", *ipLo+1 );
	}
	fprintf( ioQQQ, "\n" );

	for( vector<long>::iterator ipLo = speciesLevelList.begin();
		 ipLo != speciesLevelList.end(); ++ipLo )
	{
		if( *ipLo >= nlevels_local )
			continue;
		fprintf( ioQQQ, "%3ld\t %.4e", *ipLo+1, b[ *ipLo ] );
		for( vector<long>::iterator ipHi = speciesLevelList.begin();
			 ipHi != speciesLevelList.end(); ++ipHi )
		{
			fprintf( ioQQQ, "\t%11.4e", a[ *ipLo ][ *ipHi ] );
		}
		fprintf( ioQQQ, "\n" );
	}
}
