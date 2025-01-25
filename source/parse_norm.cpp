/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseNorm parse parameters on the normalize command */
#include "cddefines.h"
#include "lines.h"
#include "input.h"
#include "parser.h"
#include "service.h"

void ParseNorm(Parser &p)
{
	DEBUG_ENTRY( "ParseNorm()" );

	/* >>chng 01 aug 23, insist on a line label */
	/* 
	 * get possible label - must do first since it can contain a number.*/
	/* is there a double quote on the line?  if so then this is a line label */
	if( p.nMatch(  "\"" ) )
	{
		LineSave.NormLine = p.getLineID(false);
	}
	else
	{
		fprintf( ioQQQ, "The normalize command does not have a valid line.\n" );
		fprintf( ioQQQ, "A label must be specified, between double quotes, like \"H  1\" 4861.32\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	LineSave.ScaleNormLine = p.FFmtRead();
	if( p.lgEOL() )
		LineSave.ScaleNormLine = 1.;

	/* confirm that scale factor is positive */
	if( LineSave.ScaleNormLine <= 0. )
	{
		fprintf( ioQQQ, "The scale factor for relative intensities must be greater than zero.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
