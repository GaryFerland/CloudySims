/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAge parse parameters off the age command */
#include "cddefines.h"
#include "timesc.h"
#include "parser.h"

void ParseAge( Parser &p )
{
	DEBUG_ENTRY( "ParseAge()" );

	/* set age for the cloud
	 * various timescales will be checked in AgeCheck, called in comment */

	/* key " off" turns age off */
	if( p.lgEOL() && (!p.nWord(" OFF")) )
	{
		fprintf( ioQQQ, " The age must be on this line.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	timesc.CloudAgeSet = parse_input_time( p );

	return;
}
