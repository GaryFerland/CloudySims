/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* input_readarray read input commands from array where images are stored *
 * returns chCard, which will have <=80 characters before eol    *
 * line image is up and low case                                 */
/*input_init initial input_readarray array for storing line images at start of calculation */
/*lgInputComment - parse comment - check if argument is comment string */
#include "cddefines.h"
#include "trace.h"
#include "input.h"

t_input input;

void t_input::zero()
{
	DEBUG_ENTRY( "t_input::zero()" );
	/* some titles and line images */
	const int INPUT_LINE_LENGTH_EMPTY=75;
	ASSERT(INPUT_LINE_LENGTH_EMPTY <= INPUT_LINE_LENGTH);
	for( long i=0; i<INPUT_LINE_LENGTH_EMPTY-1; ++i)
	{
		chTitle[i] = ' ';
	}
	chTitle[INPUT_LINE_LENGTH_EMPTY-1] = '\0';
}

/* lgIsCommentSeq - is the string pointer s pointing to a comment sequence?
 * lgColumn0 indicates whether we are in column 0 or not; the set of valid
 *   comment characters depends on this
 * if lgReportVisible is true, visible comments will be reported, otherwise not */
bool lgIsCommentSeq( const char *s, bool lgColumn0, bool lgReportVisible )
{
	DEBUG_ENTRY( "lgIsCommentSeq()" );

	if( strncmp( s, "##", 2 ) == 0 )
		return true;
	else if( strncmp( s, "//", 2 ) == 0 || s[0] == '%' )
	{
		input.lgDeprecatedComment = true;
		return true;
	}
	else if( lgColumn0 && s[0] == '*' )
	{
		input.lgDeprecatedComment = true;
		return true;
	}
	else if( !lgColumn0 && s[0] == ';' )
	{
		input.lgDeprecatedComment = true;
		return true;
	}
	else if( lgReportVisible && s[0] == '#' )
		return true;
	else
		return false;
}

/* lgIsExpungedCommentSeq - does the string s start with an old-style comment? */
bool lgIsExpungedCommentSeq( const char *s )
{
	DEBUG_ENTRY( "lgIsExpungedCommentSeq()" );

	if( strncmp( s, "C ", 2 ) == 0 )
		return true;
	else
		return false;
}

/* lgInputComment - parse comment - check if argument is comment string, 
 * either upper or lower case -
 * returns true if line is a comment, false if not */
bool lgInputComment( const char *chLine )
{
	DEBUG_ENTRY( "lgInputComment()" );

	/* should not call this routine with null line */
	if( chLine[0] == 0 )
		TotalInsanity();

	return lgIsCommentSeq( chLine, true, true );
}

/* lgInputEOF - is this line an EOF marker? */
bool lgInputEOF( const char *chLine )
{
	if( chLine[0] == '\n' || chLine[0] == '\r' || chLine[0] == '\0' || chLine[0] == ' ' )
		return true;
	else if( strncmp( chLine, "***", 3 ) == 0 )
		return true;
	else
		return false;
}

/* StripComment- strips comment part off the command line s
 * if lgStripVisible is false, visible comments are retained
 * hidden comments are always stripped
 * this routine also removes underscores, brackets, and EOL characters */
void StripComment( string& s, bool lgStripVisible )
{
	DEBUG_ENTRY( "StripComment()" );

	for( size_t p=0; p < s.length(); )
	{
		bool lgColumn0 = ( p == 0 );
		if( s[p] == '\"' )
		{
			string buf;
			p = GetString( s, p, buf );
		}
		else if( lgIsCommentSeq(&s[p],lgColumn0,lgStripVisible) )
		{
			s.erase(p);
			break;
		}
		else if( lgIsCommentSeq(&s[p],lgColumn0,true) )
		{
			break;
		}
		else if( s[p] == '_' )
		{
			s[p++] = ' ';
			input.lgUnderscoreFound = true;
		}
		else if( s[p] == '[' || s[p] == ']' )
		{
			s[p++] = ' ';
			input.lgBracketFound = true;
		}
		else
		{
			++p;
		}
	}
	// now erase eol character. this should be a separate loop since it would
	// otherwise be skipped if GetString() doesn't find the second quote...
	for( size_t p=0; p < s.length(); ++p )
	{
		if( s[p] == '\n' || s[p] == '\r' )
		{
			s.erase(p);
			break;
		}
	}
}

/* GetString: retrieve a string between double quotes
 * s : string to be parsed
 * s[p] : place in string where to start parsing, should point to first set of double quotes
 * buf : buffer that will hold the string between double quotes
 * return value : pointer just beyond second set of double quotes, or string::npos in case of failure
 *                (second set of double quotes wasn't found) */
size_t GetString( const string& s, size_t p, string& buf )
{
	DEBUG_ENTRY( "GetString()" );

	ASSERT( s[p] == '\"' );

	buf.clear();
	for( ++p; p < s.length(); )
	{
		if( s[p] == '\\' )
			buf.push_back( GetEscape(s,p) );
		else if( s[p] == '\"' )
			return ++p;
		else
			buf.push_back( s[p++] );
	}
	// no second set of double quotes was found
	buf.clear();
	return string::npos;
}

char GetEscape( const string& s, size_t& p )
{
	DEBUG_ENTRY( "GetEscape()" );

	// This routine is the placeholder for treating character escape
	// sequences. For the moment we treat none. The routine assumes that
	// an esacpe sequence will always generate a single character. This is
	// true for all escape sequences except unicode sequences...
	// on exit, p will point one character beyond the escape sequence.
	return s[p++];
}

/*input_init initial input_readarray array for storing line images at start of calculation */
void t_input::init(void)
{
	DEBUG_ENTRY( "t_input::init()" );

	/* this sub must be called before calling READAR to get line images
	 * it simply sets the pointer to set up reading the images
	 * */
	/* this is usual case, read from the start of array, the commands */
	nRead = -1;

	return;
}

void t_input::echo( FILE *ipOUT )
{
	for( long i=0; i <= nSave; ++i )
		if( InclLevel[i] == 0 && lgVisible[i] )
			fprintf( ipOUT, "%s\n", chCardSav[i] );
}

/**input_readarray read input commands from array where images are stored *
 * returns the command in chCard */
void t_input::readarray(char *chCardStripped,
			char *chCardComment,
			bool *lgEOF)
{
	DEBUG_ENTRY( "t_input::readarray()" );

	/* usual case, reading commands from start of array
	 * nRead points to one plus the array element with the next line, it is
	 * one on the first call, which references line[0] */
	++nRead;

	/* nSave points to the last line array element that was saved,
	 * so it is one less than the number of lines read.  the last element
	 * containing a line image is [input.nSave].  There is a -1 for
	 * nRead to bring it onto the same c counting scale as nSave */
	if( nRead > nSave )
	{
		*lgEOF = true;
	}
	else
	{
		/* get the line image */
		strcpy( chCardStripped, chCardStrip[nRead] );
		strcpy( chCardComment, chCardSav[nRead] );
		*lgEOF = false;
	}

	/* if any "trace" appeared on a command line, then this flag was set
	 * so print the input command before it is parsed */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "t_input::readarray returns=%s=\n", chCardComment );
	}

	return;
}

/** input_readvector: read numbers from the file chFile and store them in a vector */
void input_readvector(const char* chFile, /**< file name to read from */
					  vector<double> &vec,/**< vector - the numbers that were read from the input line(s) */
					  bool* lgError)      /**< did an error occur reading the file? */
{
	DEBUG_ENTRY( "input_readvector()" );

	vec.clear();

	fstream ioDATA;
	open_data( ioDATA, chFile, mode_r, AS_LOCAL_ONLY );

	while( true )
	{
		double tmpValue;
		ioDATA >> tmpValue;
		if( !ioDATA.fail() )
			vec.push_back(tmpValue);
		else
			break;
	}
	*lgError = !ioDATA.eof();
}

/** input_readvector: read n numbers from the file chFile and store them in vector[] */
void input_readvector(const char* chFile, /**< file name to read from */
					  double vector[],    /**< vector[n] - the numbers that were read from the input line(s) */
					  long n,             /**< number of elements in vector[] that we need to read */
					  bool* lgEOF)        /**< was EOF reached before enough numbers were read? */
{
	DEBUG_ENTRY( "input_readvector()" );

	fstream ioDATA;
	open_data( ioDATA, chFile, mode_r, AS_LOCAL_ONLY );

	for( long i=0; i < n; ++i )
		ioDATA >> vector[i];

	*lgEOF = !ioDATA.good();
	return;
}
