/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef INPUT_H_
#define INPUT_H_

/* input.h */

#include "module.h"

/** limit to number of line images that can be read in */
const int NKRD = 4000;

/** lgIsCommentSeq - is the string pointer s pointing to a comment sequence?
 * lgColumn0 indicates whether we are in column 0 or not; the set of valid
 *   comment characters depends on this
 * if lgReportVisible is true, visible comments will be reported, otherwise not */
bool lgIsCommentSeq( const char *s, bool lgColumn0, bool lgReportVisible );

/** lgIsExpungedCommentSeq - does the string s start with an old-style comment? */
bool lgIsExpungedCommentSeq( const char *s );

/** lgInputComment - parse comment - check if argument is comment string, 
 * either upper or lower case -
 * returns true if line is a comment, false if not 
 \param *chLine the input line string
 */
bool lgInputComment( const char *chLine );

/** lgInputEOF - is this line an EOF marker? */
bool lgInputEOF( const char *chLine );

/** StripComment- strips comment part off the command line s
 * if lgStripVisible is false, visible comments are retained
 * hidden comments are always stripped
 * this routine also removes EOL characters */
void StripComment( string& s, bool lgStripVisible );

/** GetString: retrieve a string between double quotes
 * s : string to be parsed
 * s[p] : place in string where to start parsing, should point to first set of double quotes
 * buf : buffer that will hold the string between double quotes
 * return value : pointer just beyond second set of double quotes, or string::npos in case of failure
 *                (second set of double quotes wasn't found) */
size_t GetString( const string& s, size_t p, string& buf );

/** GetEscape:  This routine is the placeholder for treating character escape sequences.
 * For the moment we treat none. On exit, *p will point one character beyond the escape sequence. */
char GetEscape( const string& s, size_t& p );

/** input_readvector: read numbers from the file chFile and store them in a vector */
void input_readvector(const char* chFile, /**< file name to read from */
		      vector<double> &vec,/**< vector - the numbers that were read from the input line(s) */
		      bool* lgError);       /**< did an error occur reading the file? */

/** input_readvector: read n numbers from the file chFile and store them in vector[] */
void input_readvector(const char* chFile, /**< file name to read from */
		      double vector[],    /**< vector[n] - the numbers that were read from the input line(s) */
		      long n,             /**< number of elements in vector[] that we need to read */
		      bool* lgEOF);       /**< was EOF reached before enough numbers were read? */

struct t_input : public module {

	const char *chName() const
	{
		return "input";
	}

	void zero();

	void comment(t_warnings&) {}

	/** we will save the original (not caped) image of the line here */
	char chCardSav[NKRD][INPUT_LINE_LENGTH];

	/** same as above, but this has visible comments stripped, should only be used in parser */ 
	char chCardStrip[NKRD][INPUT_LINE_LENGTH];

	/** title entered with the title command */
	char chTitle[INPUT_LINE_LENGTH];

	/** delimiter character for file names, / for *nix, \\ for win */
	char chDelimiter[3];

	/** keep track of the include level of each input line
	 * level 0 - main input file
	 * level 1 - init file included from main input file
	 * level 2 - init file included from level 1 init file
	 *   etc...
	 * curInclLevel is the current value of the include level */
	int InclLevel[NKRD];
	int curInclLevel;

	/** is this line visible, or should it be hidden due to the
	 * HIDE keyword and/or the PRINT ON / OFF commands... */
	bool lgVisible[NKRD];

	/** current visibility status due to the PRINT ON / OFF commands */
	bool lgVisibilityStatus;

	/** one less than the total number of lines read in with cdRead */
	long int nSave;

	/** this points to the command we are now parsing, within the stack of commands */
	long int nRead;

	/** this is set true if an init file was used in the input deck */
	bool lgInitPresent;

	/** this is set true if underscore present in input stream, which was
	 * set to space */
	bool lgUnderscoreFound;

	/** this is set true if left or right bracket, [ or ], present in input stream, which was
	 * set to space */
	bool lgBracketFound;

	/** this is set true if a deprecated form of comment was used in the input script */
	bool lgDeprecatedComment;

	/** set true with no buffering command, used to print comment at end */
	bool lgSetNoBuffering;

	t_input()
	{
		// subdirectory delimiter character
#		ifdef _WIN32
		strcpy( chDelimiter, "\\" );
#		else
		strcpy( chDelimiter, "/" );
#		endif
		/* will be set true if no buffering command entered */
		lgSetNoBuffering = false;
	}


/** get the next input command off the command stack 
 * if more then copy into chCard and set lgEOF false,
 * if all command processed then set lgEOF true 
 \param *chCard the input line string
 \param *lgEOF true if hit end of file
 */
	// friend CodeSmell
	private:
	friend class Parser;
	void readarray(
		char *chCardStripped, 
		char *chCardComment, 
		bool *lgEOF);

	public:
	void echo( FILE *ipOUT);

	/** called when 'init' command hit, to reset counters for
	 * placing line images within the storage array */
	void init(void);


	};
extern t_input input;



#endif /* INPUT_H_ */
