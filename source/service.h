/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef SERVICE_H_
#define SERVICE_H_

#include <string>
#include <vector>

// remove whitespace from the end of a string
void trimTrailingWhiteSpace( std::string &str );
// remove whitespace form the end of a char array
void trimTrailingWhiteSpace( char *str );

/** split_mode defines how the routine Split generates substrings
 * SPM_RELAX: multiple adjacent separators will be coalesced into one
 *            this way you can never get an empty substring
 * SPM_KEEP_EMPTY: multiple adjacent separators will result in empty
 *                 substrings to be added to the list
 * SPM_STRICT: empty substrings are illegal */
typedef enum { SPM_RELAX, SPM_KEEP_EMPTY, SPM_STRICT } split_mode;

/** Split: split a string into substrings using "sep" as separator */
void Split(const std::string& str,   // input string
	   const std::string& sep,   // separator, may be multiple characters
	   std::vector<std::string>& lst, // the separated items will be appended here
	   split_mode mode);    // see above

inline bool FindAndReplace(string& str,
			   const string& substr,
			   const string& newstr)
{
	string::size_type ptr = str.find( substr );
	if( ptr != string::npos )
		str.replace( ptr, substr.length(), newstr );
	return ptr != string::npos;
}

inline bool FindAndErase(string& str,
			 const string& substr)
{
	return FindAndReplace( str, substr, "" );
}

void service(double tau, double a, double beta);

/** wr_block: write <len> bytes of data from buffer <*ptr> into open binary FILE* <fdes> */
inline void wr_block(const void *ptr,
		     size_t len,
		     FILE *fdes)
{
	if( fwrite(ptr,len,size_t(1),fdes) != 1 ) {
		printf( "wr_block: error writing to file\n" );
		fclose(fdes);
		cdEXIT(EXIT_FAILURE);
	}
}

/** wr_block: write <len> bytes of data from buffer <*ptr> into unformatted file <fnam> */
inline void wr_block(const void *ptr,
		     size_t len,
		     const char *fnam)
{
	FILE *fdes = open_data( fnam, "wb" );
	wr_block( ptr, len, fdes );
	fclose(fdes);
}

/** rd_block: read <len> bytes of data into buffer <*ptr> from open binary FILE* <fdes> */
inline void rd_block(void *ptr,
		     size_t len,
		     FILE *fdes)
{
	if( fread(ptr,len,size_t(1),fdes) != 1 ) {
		printf( "rd_block: error reading from file\n" );
		fclose(fdes);
		cdEXIT(EXIT_FAILURE);
	}
}

/** rd_block: read <len> bytes of data into buffer <*ptr> from unformatted file <fnam> */
inline void rd_block(void *ptr,
		     size_t len,
		     const char *fnam)
{
	FILE *fdes = open_data( fnam, "rb", AS_LOCAL_ONLY );
	rd_block( ptr, len, fdes );
	fclose(fdes);
}

#endif /* SERVICE_ */
