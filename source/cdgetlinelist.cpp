/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdGetLineList routine to read in master list of emission line wavelengths and ids, for
 * generating loc grids */
#include "cddefines.h"
#include "cddrive.h"
#include "lines.h"
#include "parser.h"

/* return value is number of lines, -1 if file could not be opened */
long int cdGetLineList( 
	/* chFile is optional filename, if void then use BLRLineList,
	 * if not void then use file specified */
	const char chFile[],
	/* array of null term strings giving line labels */
	vector<string>& chLabels,
	/* a 1-d array of line wavelengths */
	vector<realnum>& wl)
{
	DEBUG_ENTRY( "cdGetLineList()" );

	/* first check that cdInit has been called, since we may have to write
	 * error output */
	if( !lgcdInitCalled )
	{
		fprintf(stderr," cdInit must be called before cdGetLineList.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* use default filename LineList_BLR.dat if void string, else use file specified */
	const char* chFilename = ( strlen(chFile) == 0 ) ? "LineList_BLR.dat" : chFile;

	/* we will check local space first, then on path if not present */
	FILE* ioData = open_data( chFilename, "r", AS_LOCAL_DATA_TRY );

	if( ioData == NULL )
	{
		/* did not find file, return -1 */
		return -1;
	}

	// make sure we are not leaking memory
	ASSERT( chLabels.size() == 0 && wl.size() == 0 );

	Parser p;
	char chLine[FILENAME_PATH_LENGTH_2];

	/* actually read and save the lines */
	while( read_whole_line( chLine, (int)sizeof(chLine), ioData ) != NULL )
	{
		// strip eol
		for( size_t i=0; ; ++i )
		{
			if( chLine[i] == '\n' || chLine[i] == '\r' )
				chLine[i] = '\0';
			if( chLine[i] == '\0' )
				break;
		}

		p.setline(chLine);
		/* this checks for in-file EOF markers */
		if( p.last() )
			break;

		if( p.isComment() )
			continue;

		LineID line = p.getLineID();
		if( !p.lgReachedEnd() )
		{
			fprintf( ioQQQ, "cdGetLineList: found junk at end of input line in file %s:\n", chFile );
			p.showLocation();
			cdEXIT(EXIT_FAILURE);
		}
		chLabels.push_back(line.chLabel);
		wl.push_back(line.wave);
	}

	fclose( ioData );

	/* return number of lines we found */
	return chLabels.size();
}
