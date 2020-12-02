/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*wcnint initialize stack or warnings, cautions, notes */
/*warnin enter warnings at the end of the calculations into large stack */
/*notein enter a note about calculation into comment array */
/*bangin called by routine comment to enter surprise into comment stack */
/*caunin called by comment to enter caution into comment stack */
#include "cddefines.h"
#include "warnings.h"

t_warnings warnings;

void t_warnings::zero(void)
{

	DEBUG_ENTRY( "t_warnings::zero()" );

	/* this sub is called first, to initialize the variables */
	nwarn = 0;
	ncaun = 0;
	nnote = 0;
	nbang = 0;
	return;
}

/*warnin enter warnings at the end of the calculations into large stack */
void t_warnings::warnin(const char *chLine)
{

	DEBUG_ENTRY( "t_warnings::warnin()" );

	if( nwarn >= LIMWCN )
	{
		static bool lgFirst=true;
		if( lgFirst )
			fprintf( ioQQQ, 
				" Too many warnings have been entered; increase the value of LIMWCN everywhere in the code.\n" );
		lgFirst = false;
	}
	else
	{
		strcpy( chWarnln[nwarn], chLine  );
	}

	++nwarn;
	return;
}

/*notein enter a note about calculation into comment array */
void t_warnings::notein(const char *chLine)
{

	DEBUG_ENTRY( "notein()" );

	if( nnote >= LIMWCN )
	{
		static bool lgFirst=true;
		if( lgFirst )
			fprintf( ioQQQ, 
				" Too many notes have been entered; increase the value of LIMWCN everywhere in the code.\n" );
		lgFirst = false;
	}
	else
	{
		strcpy( chNoteln[nnote], chLine );
	}

	++nnote;
	return;
}

/*bangin called by routine comment to enter surprise into comment stack */
void t_warnings::bangin(const char *chLine)
{

	DEBUG_ENTRY( "t_warnings::bangin()" );

	if( nbang >= LIMWCN )
	{
		static bool lgFirst=true;
		if( lgFirst )
			fprintf( ioQQQ, 
				" Too many surprises have been entered; increase the value of LIMWCN everywhere in the code.\n" );
		lgFirst = false;
	}
	else
	{
		strcpy( chBangln[nbang], chLine );
	}

	++nbang;
	return;
}

/*caunin called by comment to enter caution into comment stack */
void t_warnings::caunin(const char *chLine)
{

	DEBUG_ENTRY( "t_warnings::caunin()" );

	if( ncaun >= LIMWCN )
	{
		static bool lgFirst=true;
		if( lgFirst )
			fprintf( ioQQQ, 
				" Too many cautions have been entered; increase the value of LIMWCN everywhere in the code.\n" );
		lgFirst = false;
	}
	else
	{
		strcpy( chCaunln[ncaun], chLine );
	}

	++ncaun;
	return;
}
