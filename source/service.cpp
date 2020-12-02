/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* 
* a set of  routines that are widely used across the code for various
* housekeeping chores.  These do not do any physics and are unlikely to
* change over time.  The prototypes are in cddefines.h and so are 
* automatically picked up by all routines 
*/
/*FFmtRead scan input line for free format number */
/*caps convert input command line (through eol) to ALL CAPS */
/*ShowMe produce request to send information to GJF after a crash */
/*AnuUnit produce continuum energy in arbitrary units */
/*cap4 convert first 4 char of input line chLab into chCAP all in caps, null termination */
/*insane set flag saying that insanity has occurred */
/*nMatch determine whether match to a keyword occurs on command line,
 * return value is 0 if no match, and position of match within string if hit */
/*fudge enter fudge factors, or some arbitrary number, with fudge command*/
/*qip compute pow(x,n) for positive integer n through repeated squares */
/*dsexp safe exponential function for doubles */
/*sexp safe exponential function */
/*TestCode set flag saying that test code is in place */
/*CodeReview - placed next to code that needs to be checked */
/*fixit - say that code needs to be fixed */
/*broken set flag saying that the code is broken, */
/*dbg_printf is a debug print routine that was provided by Peter Teuben,
 * as a component from his NEMO package.  It offers run-time specification
 * of the level of debugging */
/*qg32 32 point Gaussian quadrature, original Fortran given to Gary F by Jim Lattimer */
/*TotalInsanity general error handler for something that cannot happen */
/*BadRead general error handler for trying to read data, but failing */
/*MyMalloc wrapper for malloc().  Returns a good pointer or dies. */
/*spsort netlib routine to sort array returning sorted indices */
/*chLineLbl use information in line transfer arrays to generate a line label *
 * this label is null terminated */
/*chIonLbl use information in line array to generate a null terminated ion label in "Fe 2" */
/*csphot returns photoionization cross section from opacity stage using std pointers */
/*MyAssert a version of assert that fails gracefully */
/*RandGauss normal random variate generator */
/*MyGaussRand a wrapper for RandGauss, see below */

#include "cdstd.h"
#include <cstdarg>	/* ANSI variable arg macros */
#include "cddefines.h"
#include "service.h"
#include "cddrive.h"
#include "called.h"
#include "opacity.h"
#include "rfield.h"
#include "hextra.h"
#include "struc.h"
#include "fudgec.h"
#include "broke.h"
#include "trace.h"
#include "input.h"
#include "save.h"
#include "version.h"
#include "warnings.h"
#include "conv.h"
#include "atmdat.h"
#include "mole.h"
#include "prt.h"
#include "integrate.h"

#ifdef __CYGWIN__
extern "C" { int vsnprintf(char*, size_t, const char*, va_list); }
#endif

/*read_whole_line - safe version of fgets - read a line, 
 * return null if cannot read line or if input line is too long */
char *read_whole_line( char *chLine , int nChar , FILE *ioIN )
{
	char *chRet;

	DEBUG_ENTRY( "read_whole_line()" );

	/* wipe the buffer to prevent the code from accidentally picking up on previous input */
	memset( chLine, 0, (size_t)nChar );

	/* this always writes a '\0' character, even if line does not fit in buffer
	 * the terminating newline is copied only if the line does fit in the buffer */
	if( (chRet = fgets( chLine, nChar, ioIN )) == NULL )
		return NULL;

	long lineLength = strlen( chRet );
	//fprintf(ioQQQ , "DEBUG reading:%s\n" , chLine);
	//fprintf(ioQQQ , "DEBUG length is %li nChar is %i \n", lineLength , nChar);

	/* return null if input string is longer than nChar-1 (including terminating newline),
	 * the longest we can read. Print and return null but chLine still has as much of
	 * the line as could be placed in the buffer */
	if( lineLength>=nChar-1 )
	{
		if( called.lgTalk )
			fprintf( ioQQQ, "DISASTER PROBLEM read_whole_line found input"
			" with a line too long to be read, limit is %i char.  "
			"Start of line follows:\n%s\n",
			nChar , chLine );

		lgAbort = true;
		return NULL;
	}
	return chRet;
}

/** Split: split a string into substrings using "sep" as separator */
void Split(const string& str,   // input string
	   const string& sep,   // separator, may be multiple characters
	   vector<string>& lst, // the separated items will be appended here
	   split_mode mode)     // SPM_RELAX, SPM_KEEP_EMPTY, or SPM_STRICT; see cddefines.h
{
	DEBUG_ENTRY( "Split()" );

	bool lgStrict = ( mode == SPM_STRICT );
	bool lgKeep = ( mode == SPM_KEEP_EMPTY );
	bool lgFail = false;
	string::size_type ptr1 = 0;
	string::size_type ptr2 = str.find( sep );
	string sstr = str.substr( ptr1, ptr2-ptr1 );
	if( sstr.length() > 0 )
		lst.push_back( sstr );
	else {
		if( lgStrict ) lgFail = true;
		if( lgKeep ) lst.push_back( sstr );
	}
	while( ptr2 != string::npos ) {
		// the separator is skipped
		ptr1 = ptr2 + sep.length();
		if( ptr1 < str.length() ) {
			ptr2 = str.find( sep, ptr1 );
			sstr = str.substr( ptr1, ptr2-ptr1 );
			if( sstr.length() > 0 )
				lst.push_back( sstr );
			else {
				if( lgStrict ) lgFail = true;
				if( lgKeep ) lst.push_back( sstr );
			}
		}
		else {
			ptr2 = string::npos;
			if( lgStrict ) lgFail = true;
			if( lgKeep ) lst.push_back( "" );
		}
	}
	if( lgFail )
	{
		fprintf( ioQQQ, " A syntax error occurred while splitting the string: \"%s\"\n", str.c_str() );
		fprintf( ioQQQ, " The separator is \"%s\". Empty substrings are not allowed.\n", sep.c_str() );
		cdEXIT(EXIT_FAILURE);
	}
}

// remove whitespace from the end of a string
void trimTrailingWhiteSpace( string &str )
{
	size_t pos = str.find_last_not_of(" \t");
	// If this fails, everything is whitespace.
	if ( pos != string::npos )
		str.erase( pos+1 );
	else
		str.clear();
	return;
}

// remove whitespace from the end of a string
void trimTrailingWhiteSpace( char *str )
{
	int pos = strlen( str );
	while( pos > 0 && (str[pos-1]==' ' || str[pos-1]=='\t' ))
		--pos;
	str[pos] = '\0';
	// If this fails, everything is whitespace.
	// ASSERT( pos != 0 ); 
	return;
}

/* a version of assert that fails gracefully */
void MyAssert(const char *file, int line, const char *comment)
{
	DEBUG_ENTRY( "MyAssert()" );

	fprintf(ioQQQ,"\n\n\n PROBLEM DISASTER\n An assert has been thrown, this is bad.\n");
	fprintf(ioQQQ," %s\n",comment);
	fprintf(ioQQQ," It happened in the file %s at line number %i\n", file, line );
	fprintf(ioQQQ," This is iteration %li, nzone %li, fzone %.2f, lgSearch=%c.\n", 
		iteration , 
		nzone ,
		fnzone ,
		TorF(conv.lgSearch) );

	ShowMe();
#	ifdef OLD_ASSERT
	cdEXIT(EXIT_FAILURE);
#	endif
}

/*AnuUnit produce continuum energy in arbitrary units, as determined by ChkUnits() */
double AnuUnit(realnum energy_ryd)
{
	DEBUG_ENTRY( "AnuUnit()" );

	return Energy((double)energy_ryd).get(save.chConSavEnr[save.ipConPun]);
}

/*ShowMe produce request to send information to GJF after a crash */
void ShowMe(void)
{
	DEBUG_ENTRY( "ShowMe()" );

	/* print info if output unit is defined */
	if( ioQQQ != NULL )
	{
		/* >>chng 06 mar 02 - check if molecular but cosmic rays are ignored */
		molezone* h2 = findspecieslocal("H2");
		// molecular species may not be set up yet, so check for NULL pointer...
		if( (hextra.cryden == 0.) && h2 != NULL && h2->xFracLim > 0.1 )
		{
			fprintf( ioQQQ, " >>> \n >>> \n >>> Cosmic rays are not included and the gas is molecular.  "
				"THIS IS KNOWN TO BE UNSTABLE.  Add cosmic rays and try again.\n >>> \n >>>\n\n");
		}
		else
		{
			fprintf( ioQQQ, "\n\n\n" );
			fprintf( ioQQQ, "           vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv \n" );
			fprintf( ioQQQ, "          > PROBLEM DISASTER PROBLEM DISASTER.      <\n" );
			fprintf( ioQQQ, "          > Sorry, something bad has happened.      <\n" );
			fprintf( ioQQQ, "          > Please post this on the Cloudy web site <\n" );
			fprintf( ioQQQ, "          > discussion board at www.nublado.org     <\n" );
			fprintf( ioQQQ, "          > Please send all following information:  <\n" );
			fprintf( ioQQQ, "           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n" );
			fprintf( ioQQQ, "\n\n" );


			fprintf( ioQQQ, " Cloudy version number is %s\n", 
				t_version::Inst().chVersion );
			fprintf( ioQQQ, " %s\n\n", t_version::Inst().chInfo );

			fprintf( ioQQQ, "%5ld warnings,%3ld cautions,%3ld temperature failures.  Messages follow.\n", 
			  warnings.nwarn, warnings.ncaun, conv.nTeFail );

			/* print the warnings first */
			cdWarnings(ioQQQ);

			/* now print the cautions */
			cdCautions(ioQQQ);

			/* now output the commands */
			cdPrintCommands(ioQQQ);

			if( input.lgInitPresent )
			{
				fprintf(ioQQQ," This input stream included one or more init files.\n");
				fprintf(ioQQQ," If these files are not part of the standard Cloudy distribution\n"); 
				fprintf(ioQQQ," then I will need a copy of them too.\n");
			}
		}
	}
	return;
}

/*cap4 convert first 4 char of input line chLab into chCAP all in caps, null termination */
void cap4(char *chCAP ,	/* output string, cap'd first 4 char of chLab, with null terminating */
	  const char *chLab)	/* input string ending with null byte */
{
	DEBUG_ENTRY( "cap4()" );

	/* convert 4 character string in chLab to ALL CAPS in chCAP */
	for( long i=0; i < 4; i++ )
	{
		/* toupper is function in ctype that converts to upper case */
		chCAP[i] = toupper( chLab[i] );
		if( chLab[i] == '\0' )
			break;
	}

	/* now end string with null byte */
	chCAP[4] = '\0';
	return;
}

/*uncaps convert input command line (through eol) to all lowercase */
void uncaps(char *chCard )
{
	DEBUG_ENTRY( "uncaps()" );

	long i = 0;
	while( chCard[i] != '\0' )
	{
		chCard[i] = tolower( chCard[i] );
		++i;
	}
	return;
}

/*caps convert input command line (through eol) to all uppercase */
void caps(char *chCard )
{
	DEBUG_ENTRY( "caps()" );

	long i = 0;
	while( chCard[i] != '\0' )
	{
		chCard[i] = toupper( chCard[i] );
		++i;
	}
	return;
}

static const double pos_pow10[] = {
	1.e+000, 1.e+001, 1.e+002, 1.e+003, 1.e+004, 1.e+005, 1.e+006, 1.e+007, 1.e+008, 1.e+009,
	1.e+010, 1.e+011, 1.e+012, 1.e+013, 1.e+014, 1.e+015, 1.e+016, 1.e+017, 1.e+018, 1.e+019,
	1.e+020, 1.e+021, 1.e+022, 1.e+023, 1.e+024, 1.e+025, 1.e+026, 1.e+027, 1.e+028, 1.e+029,
	1.e+030, 1.e+031, 1.e+032, 1.e+033, 1.e+034, 1.e+035, 1.e+036, 1.e+037, 1.e+038, 1.e+039,
	1.e+040, 1.e+041, 1.e+042, 1.e+043, 1.e+044, 1.e+045, 1.e+046, 1.e+047, 1.e+048, 1.e+049,
	1.e+050, 1.e+051, 1.e+052, 1.e+053, 1.e+054, 1.e+055, 1.e+056, 1.e+057, 1.e+058, 1.e+059,
	1.e+060, 1.e+061, 1.e+062, 1.e+063, 1.e+064, 1.e+065, 1.e+066, 1.e+067, 1.e+068, 1.e+069,
	1.e+070, 1.e+071, 1.e+072, 1.e+073, 1.e+074, 1.e+075, 1.e+076, 1.e+077, 1.e+078, 1.e+079,
	1.e+080, 1.e+081, 1.e+082, 1.e+083, 1.e+084, 1.e+085, 1.e+086, 1.e+087, 1.e+088, 1.e+089,
	1.e+090, 1.e+091, 1.e+092, 1.e+093, 1.e+094, 1.e+095, 1.e+096, 1.e+097, 1.e+098, 1.e+099,
	1.e+100, 1.e+101, 1.e+102, 1.e+103, 1.e+104, 1.e+105, 1.e+106, 1.e+107, 1.e+108, 1.e+109,
	1.e+110, 1.e+111, 1.e+112, 1.e+113, 1.e+114, 1.e+115, 1.e+116, 1.e+117, 1.e+118, 1.e+119,
	1.e+120, 1.e+121, 1.e+122, 1.e+123, 1.e+124, 1.e+125, 1.e+126, 1.e+127, 1.e+128, 1.e+129,
	1.e+130, 1.e+131, 1.e+132, 1.e+133, 1.e+134, 1.e+135, 1.e+136, 1.e+137, 1.e+138, 1.e+139,
	1.e+140, 1.e+141, 1.e+142, 1.e+143, 1.e+144, 1.e+145, 1.e+146, 1.e+147, 1.e+148, 1.e+149,
	1.e+150, 1.e+151, 1.e+152, 1.e+153, 1.e+154, 1.e+155, 1.e+156, 1.e+157, 1.e+158, 1.e+159,
	1.e+160, 1.e+161, 1.e+162, 1.e+163, 1.e+164, 1.e+165, 1.e+166, 1.e+167, 1.e+168, 1.e+169,
	1.e+170, 1.e+171, 1.e+172, 1.e+173, 1.e+174, 1.e+175, 1.e+176, 1.e+177, 1.e+178, 1.e+179,
	1.e+180, 1.e+181, 1.e+182, 1.e+183, 1.e+184, 1.e+185, 1.e+186, 1.e+187, 1.e+188, 1.e+189,
	1.e+190, 1.e+191, 1.e+192, 1.e+193, 1.e+194, 1.e+195, 1.e+196, 1.e+197, 1.e+198, 1.e+199,
	1.e+200, 1.e+201, 1.e+202, 1.e+203, 1.e+204, 1.e+205, 1.e+206, 1.e+207, 1.e+208, 1.e+209,
	1.e+210, 1.e+211, 1.e+212, 1.e+213, 1.e+214, 1.e+215, 1.e+216, 1.e+217, 1.e+218, 1.e+219,
	1.e+220, 1.e+221, 1.e+222, 1.e+223, 1.e+224, 1.e+225, 1.e+226, 1.e+227, 1.e+228, 1.e+229,
	1.e+230, 1.e+231, 1.e+232, 1.e+233, 1.e+234, 1.e+235, 1.e+236, 1.e+237, 1.e+238, 1.e+239,
	1.e+240, 1.e+241, 1.e+242, 1.e+243, 1.e+244, 1.e+245, 1.e+246, 1.e+247, 1.e+248, 1.e+249,
	1.e+250, 1.e+251, 1.e+252, 1.e+253, 1.e+254, 1.e+255, 1.e+256, 1.e+257, 1.e+258, 1.e+259,
	1.e+260, 1.e+261, 1.e+262, 1.e+263, 1.e+264, 1.e+265, 1.e+266, 1.e+267, 1.e+268, 1.e+269,
	1.e+270, 1.e+271, 1.e+272, 1.e+273, 1.e+274, 1.e+275, 1.e+276, 1.e+277, 1.e+278, 1.e+279,
	1.e+280, 1.e+281, 1.e+282, 1.e+283, 1.e+284, 1.e+285, 1.e+286, 1.e+287, 1.e+288, 1.e+289,
	1.e+290, 1.e+291, 1.e+292, 1.e+293, 1.e+294, 1.e+295, 1.e+296, 1.e+297, 1.e+298, 1.e+299,
	1.e+300, 1.e+301, 1.e+302, 1.e+303, 1.e+304, 1.e+305, 1.e+306, 1.e+307, 1.e+308
};

static const int max_pow10 = int(sizeof(pos_pow10)/sizeof(pos_pow10[0]) - 1);

static const double neg_pow10[] = {
	1.e-000, 1.e-001, 1.e-002, 1.e-003, 1.e-004, 1.e-005, 1.e-006, 1.e-007, 1.e-008, 1.e-009,
	1.e-010, 1.e-011, 1.e-012, 1.e-013, 1.e-014, 1.e-015, 1.e-016, 1.e-017, 1.e-018, 1.e-019,
	1.e-020, 1.e-021, 1.e-022, 1.e-023, 1.e-024, 1.e-025, 1.e-026, 1.e-027, 1.e-028, 1.e-029,
	1.e-030, 1.e-031, 1.e-032, 1.e-033, 1.e-034, 1.e-035, 1.e-036, 1.e-037, 1.e-038, 1.e-039,
	1.e-040, 1.e-041, 1.e-042, 1.e-043, 1.e-044, 1.e-045, 1.e-046, 1.e-047, 1.e-048, 1.e-049,
	1.e-050, 1.e-051, 1.e-052, 1.e-053, 1.e-054, 1.e-055, 1.e-056, 1.e-057, 1.e-058, 1.e-059,
	1.e-060, 1.e-061, 1.e-062, 1.e-063, 1.e-064, 1.e-065, 1.e-066, 1.e-067, 1.e-068, 1.e-069,
	1.e-070, 1.e-071, 1.e-072, 1.e-073, 1.e-074, 1.e-075, 1.e-076, 1.e-077, 1.e-078, 1.e-079,
	1.e-080, 1.e-081, 1.e-082, 1.e-083, 1.e-084, 1.e-085, 1.e-086, 1.e-087, 1.e-088, 1.e-089,
	1.e-090, 1.e-091, 1.e-092, 1.e-093, 1.e-094, 1.e-095, 1.e-096, 1.e-097, 1.e-098, 1.e-099,
	1.e-100, 1.e-101, 1.e-102, 1.e-103, 1.e-104, 1.e-105, 1.e-106, 1.e-107, 1.e-108, 1.e-109,
	1.e-110, 1.e-111, 1.e-112, 1.e-113, 1.e-114, 1.e-115, 1.e-116, 1.e-117, 1.e-118, 1.e-119,
	1.e-120, 1.e-121, 1.e-122, 1.e-123, 1.e-124, 1.e-125, 1.e-126, 1.e-127, 1.e-128, 1.e-129,
	1.e-130, 1.e-131, 1.e-132, 1.e-133, 1.e-134, 1.e-135, 1.e-136, 1.e-137, 1.e-138, 1.e-139,
	1.e-140, 1.e-141, 1.e-142, 1.e-143, 1.e-144, 1.e-145, 1.e-146, 1.e-147, 1.e-148, 1.e-149,
	1.e-150, 1.e-151, 1.e-152, 1.e-153, 1.e-154, 1.e-155, 1.e-156, 1.e-157, 1.e-158, 1.e-159,
	1.e-160, 1.e-161, 1.e-162, 1.e-163, 1.e-164, 1.e-165, 1.e-166, 1.e-167, 1.e-168, 1.e-169,
	1.e-170, 1.e-171, 1.e-172, 1.e-173, 1.e-174, 1.e-175, 1.e-176, 1.e-177, 1.e-178, 1.e-179,
	1.e-180, 1.e-181, 1.e-182, 1.e-183, 1.e-184, 1.e-185, 1.e-186, 1.e-187, 1.e-188, 1.e-189,
	1.e-190, 1.e-191, 1.e-192, 1.e-193, 1.e-194, 1.e-195, 1.e-196, 1.e-197, 1.e-198, 1.e-199,
	1.e-200, 1.e-201, 1.e-202, 1.e-203, 1.e-204, 1.e-205, 1.e-206, 1.e-207, 1.e-208, 1.e-209,
	1.e-210, 1.e-211, 1.e-212, 1.e-213, 1.e-214, 1.e-215, 1.e-216, 1.e-217, 1.e-218, 1.e-219,
	1.e-220, 1.e-221, 1.e-222, 1.e-223, 1.e-224, 1.e-225, 1.e-226, 1.e-227, 1.e-228, 1.e-229,
	1.e-230, 1.e-231, 1.e-232, 1.e-233, 1.e-234, 1.e-235, 1.e-236, 1.e-237, 1.e-238, 1.e-239,
	1.e-240, 1.e-241, 1.e-242, 1.e-243, 1.e-244, 1.e-245, 1.e-246, 1.e-247, 1.e-248, 1.e-249,
	1.e-250, 1.e-251, 1.e-252, 1.e-253, 1.e-254, 1.e-255, 1.e-256, 1.e-257, 1.e-258, 1.e-259,
	1.e-260, 1.e-261, 1.e-262, 1.e-263, 1.e-264, 1.e-265, 1.e-266, 1.e-267, 1.e-268, 1.e-269,
	1.e-270, 1.e-271, 1.e-272, 1.e-273, 1.e-274, 1.e-275, 1.e-276, 1.e-277, 1.e-278, 1.e-279,
	1.e-280, 1.e-281, 1.e-282, 1.e-283, 1.e-284, 1.e-285, 1.e-286, 1.e-287, 1.e-288, 1.e-289,
	1.e-290, 1.e-291, 1.e-292, 1.e-293, 1.e-294, 1.e-295, 1.e-296, 1.e-297, 1.e-298, 1.e-299,
	1.e-300, 1.e-301, 1.e-302, 1.e-303, 1.e-304, 1.e-305, 1.e-306, 1.e-307
};

static const int min_pow10 = -int(sizeof(neg_pow10)/sizeof(neg_pow10[0]) - 1);

/*FFmtRead scan input line for free format number */
double FFmtRead(const char *chCard, 
		long int *ipnt, 
		long int length,
		bool *lgEOL)
{
	DEBUG_ENTRY( "FFmtRead()" );

	char chr = '\0';
	const char *eol_ptr = &chCard[length]; // eol_ptr points one beyond end of buffer
	const char *ptr = min(&chCard[*ipnt-1],eol_ptr); // ipnt is on fortran scale!

	if( *ipnt <= 0 )
	{
		fprintf(ioQQQ, "PROBLEM FFmtRead called with index <= 0, ipnt is %li\n",*ipnt);
		fprintf(ioQQQ, "Line image: %s\n", chCard);
		TotalInsanity();
	}
	else if( *ipnt > length )
	{
		fprintf(ioQQQ, "PROBLEM FFmtRead called with index > length, ipnt is %li length is %li\n",
				*ipnt, length);
		fprintf(ioQQQ, "Line image: %s\n", chCard);
		TotalInsanity();
	}

	while( true )
	{
		if( ptr >= eol_ptr || ( chr = *ptr++ ) == '\0' )
		{
			*ipnt = length+1;
			*lgEOL = true;
			return 0.;
		}
		const char *lptr = ptr;
		char lchr = chr;
		if( ( lchr == '-' || lchr == '+' ) && lptr < eol_ptr )
			lchr = *lptr++;
		if( lchr == '.'  && lptr < eol_ptr )
			lchr = *lptr;
		if( isdigit(lchr) )
			break;
	}

	//double atofval = atof(ptr-1);

	double number = 0.0;
	int exponent=0, sign=1, expsign=1, scale=0;
	bool lgCommaFound = false, lgLastComma = false, foundpoint = false, foundexp = false;
	do
	{
		lgCommaFound = lgLastComma;
		if( chr == ',' )
		{
			/* don't complain about comma if it appears after number,
			   as determined by exiting loop before this sets lgCommaFound */
			lgLastComma = true;
		}
		else if (isdigit(chr))
		{
			int digit = (chr - '0');
			if (foundexp)
			{
				exponent = 10*exponent+digit;
			}
			else 
			{
				number = 10.0*number+digit;
				if (foundpoint)
					++scale;
			}
		}
		else if (chr == '-')
		{
			if (foundexp)
			{
				if (exponent != 0 || expsign != 1)
					break;
				expsign = -1;
			}
			else
			{
				if (number != 0 || sign != 1)
					break;
				sign = -1;
			}
		}
		else if (tolower(chr) == 'e')
		{
			if (foundexp)
				break;
			foundexp = true;
		}
		else if (chr == '.')
		{
			if (foundpoint)
				break;
			foundpoint = true;
		}
		if( ptr == eol_ptr )
			break;
		chr = *ptr++;
	}
	while( isdigit(chr) || chr == '.' || chr == '-' || chr == '+' || chr == ',' || chr == 'e' || chr == 'E' );

	if( lgCommaFound )
	{
		fprintf( ioQQQ, " PROBLEM - a comma was found embedded in a number, this is deprecated.\n" );
		fprintf(ioQQQ, "== %-80s ==\n",chCard);
	}

	int expo = expsign*exponent-scale;
	double value = sign*number;
	// numbers produced by FFmtRead() should ideally be accurate to 1 ULP, but certainly
	// better than 3 ULP (the rest of the code relies on this).
	// To achieve this we use a lookup table of powers of 10, which is also fast...
	if( expo >= 0 )
	{
		while( expo > max_pow10 )
		{
			value *= pos_pow10[max_pow10];
			expo -= max_pow10;
		}
		value *= pos_pow10[expo];
	}
	else
	{
		while( expo < min_pow10 )
		{
			value *= neg_pow10[-min_pow10];
			expo -= min_pow10;
		}
		value *= neg_pow10[-expo];
	}

	//ASSERT(fp_equal(value,atofval,2));
	//fprintf(ioQQQ,"%g %g %g\n",value == 0 ? atofval : atofval/value-1.,value,atofval);
	
	*ipnt = (long)(ptr - chCard); // ptr already points 1 beyond where next read should start
	*lgEOL = false;
	return value;
}

/*nMatch determine whether match to a keyword occurs on command line,
 * return value is 0 if no match, and position of match within string if hit */
long nMatch(const char *chKey, 
	    const char *chCard)
{
	const char *ptr;
	long Match_v;

	DEBUG_ENTRY( "nMatch()" );

	ASSERT( strlen(chKey) > 0 );

	if( ( ptr = strstr_s( chCard, chKey ) ) == NULL )
	{
		/* did not find match, return 0 */
		Match_v = 0L;
	}
	else
	{
		/* return position within chCard (fortran scale) */
		Match_v = (long)(ptr-chCard+1);
	}
	return Match_v;
}

/* fudge enter fudge factors, or some arbitrary number, with fudge command
 * other sections of the code access these numbers by calling fudge
 * fudge(0) returns the first number that was entered
 * prototype for this function is in cddefines.h so that it can be used without
 * declarations 
 * fudge(-1) queries the routine for the number of fudge parameters that were entered,
 * zero returned if none */
double fudge(long int ipnt)
{
	double fudge_v;

	DEBUG_ENTRY( "fudge()" );

	if( ipnt < 0 )
	{
		/* this is special case, return number of arguments */
		fudge_v = fudgec.nfudge;
		fudgec.lgFudgeUsed = true;
	}
	else if( ipnt >= fudgec.nfudge )
	{
		fprintf( ioQQQ, " FUDGE factor not entered for array number %3ld\n", 
		  ipnt );
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		fudge_v = fudgec.fudgea[ipnt];
		fudgec.lgFudgeUsed = true;
	}
	return fudge_v;
}

/* want to define this only if no native os support exists */
#ifndef HAVE_POWI

/* powi.c - calc x^n, where n is an integer! */

/* Very slightly modified version of power() from Computer Language, Sept. 86,
	pg 87, by Jon Snader (who extracted the binary algorithm from Donald Knuth,
	"The Art of Computer Programming", vol 2, 1969).
	powi() will only be called when an exponentiation with an integer
	exponent is performed, thus tests & code for fractional exponents were 
	removed.
 */

double powi( double x, long int n )	/* returns:  x^n */
/* x;	 base */
/* n;	 exponent */
{
	double p;	/* holds partial product */

	DEBUG_ENTRY( "powi()" );

	if( x == 0 )
		return 0.;

	/* test for negative exponent */
	if( n < 0 )
	{	
		n = -n;
		x = 1/x;
	}

	p = is_odd(n) ? x : 1;	/* test & set zero power */

	/*lint -e704 shift right of signed quantity */
	/*lint -e720 Boolean test of assignment */
	while( n >>= 1 )
	{	/* now do the other powers */
		x *= x;			/* sq previous power of x */
		if( is_odd(n) )	/* if low order bit set */
			p *= x;		/*	then, multiply partial product by latest power of x */
	}
	/*lint +e704 shift right of signed quantity */
	/*lint +e720 Boolean test of assignment */
	return p;
}

#endif /* HAVE_POWI */

/* efficiently calculate x^(p/q) */
double powpq(double x, int p, int q)
{
	DEBUG_ENTRY( "powpq()" );

	if( q < 0 )
	{
		p = -p;
		q = -q;
	}

	if( q == 1 )
		return powi(x,p);
	else if( q == 2 )
		return powi(sqrt(x),p);
	else if( q == 3 )
		return powi(cbrt(x),p);
	else if( q == 4 )
		return powi(sqrt(sqrt(x)),p);
	else if( q == 6 )
		return powi(sqrt(cbrt(x)),p);
	else if( q == 8 )
		return powi(sqrt(sqrt(sqrt(x))),p);
	else if( q == 9 )
		return powi(cbrt(cbrt(x)),p);
	else
		return pow(x,double(p)/double(q));
}

long ipow( long m, long n )	/* returns:  m^n */
/* m;		 base */
/* n;		 exponent */
{
	long p;	/* holds partial product */

	DEBUG_ENTRY( "ipow()" );

	if( m == 0 || (n < 0 && m > 1) )
		return 0L;
	/* NOTE: negative exponent always results in 0 for integers!
	 * (except for the case when m==1 or -1) */

	if( n < 0 )
	{	/* test for negative exponent */
		n = -n;
		m = 1/m;
	}

	p = is_odd(n) ? m : 1;	/* test & set zero power */

	/*lint -e704 shift right of signed quantity */
	/*lint -e720 Boolean test of assignment */
	while( n >>= 1 )
	{	/* now do the other powers */
		m *= m;			/* sq previous power of m */
		if( is_odd(n) )	/* if low order bit set */
			p *= m;		/*	then, multiply partial product by latest power of m */
	}
	/*lint +e704 shift right of signed quantity */
	/*lint +e720 Boolean test of assignment */
	return p;
}

#ifndef HAVE_STRNLEN

// this routine cannot clash with library version since it has C++ linkage
size_t strnlen(const char *s, size_t maxlen)
{
	for( size_t i=0; i < maxlen; ++i )
		if( s[i] == '\0' )
			return i;
	return maxlen;
}

#endif

//
// sncatf() is fully equivalent to snprintf() apart from the fact that it
// concatenates the output to an existing string rather than replacing it.
//
// this routine was taken from
// http://stackoverflow.com/questions/2674312/how-to-append-strings-using-sprintf
//
// the return value is the length the string in buf[] would have had _if_ buf[]
// would have been large enough to hold it. Hence this value can be used to test
// for truncation: if ret_val >= bufSize then the string in buf[] was truncated.
//
size_t sncatf( char* buf, size_t bufSize, const char* fmt, ... )
{
	size_t result;
	va_list args;
	size_t len = strnlen( buf, bufSize );

	va_start( args, fmt );
	result = vsnprintf( buf + len, bufSize - len, fmt, args );
	va_end( args );

	return result + len;
}

size_t sncatf( ostringstream& str, const char* fmt, ... )
{
	DEBUG_ENTRY( "sncatf()" );
	size_t result;
	size_t len = str.tellp();
	va_list args;
	char tmp[64];
	char *tmp1=tmp;

	va_start( args, fmt );
	result = vsnprintf( tmp, sizeof(tmp), fmt, args );
	va_end( args );

	if (result >= sizeof(tmp))
	{
		tmp1 = new char[result+1];
		va_start( args, fmt );
		result = vsnprintf( tmp1, result+1, fmt, args );
		va_end( args );
	}
	
	str << tmp1;

	if (tmp1 != tmp)
		delete [] tmp1;

	return result+len;
}

/*PrintE82 - series of routines to mimic 1p, e8.2 fortran output */
/***********************************************************
 * contains the following sets of routines to get around   *
 * the MS C++ compilers unusual exponential output.        *
 * PrintEfmt <= much faster, no overhead with unix         *
 * PrintE93                                                *
 * PrintE82                                                *
 * PrintE71                                                *
 **********************************************************/

#ifdef _MSC_VER
/**********************************************************/
/*
 * Instead of printf'ing with %e or %.5e or whatever, call
 * efmt("%e", val) and print the result with %s.  This lets
 * us work around bugs in VS C 6.0.
 */
char *PrintEfmt(const char *fmt, double val /* , char *buf */) 
{
	static char buf[30]; /* or pass it in */

	DEBUG_ENTRY( "PrintEfmt()" );

	/* create the string */
	sprintf(buf, fmt, val);

	/* code to fix incorrect ms v e format.  works only for case where there is
	 * a leading space in the format - for formats like 7.1, 8.2, 9.3, 10.4, etc
	 * result will have 1 too many characters */
	char *ep , buf2[30];

	/* msvc behaves badly in different ways for positive vs negative sign vals,
	 * if val is positive must create a leading space */
	if( val >= 0.)
	{
		strcpy(buf2 , " " );
		strcat(buf2 , buf);
		strcpy( buf , buf2);
	}

	/* allow for both e and E formats */
	if((ep = strchr_s(buf, 'e')) == NULL)
	{
		ep = strchr_s(buf, 'E');
	}

	/* ep can still be NULL if val is Inf or NaN */
	if(ep != NULL) 
	{
		/* move pointer two char past the e, to pick up the e and sign */
		ep += 2;

		/* terminate buf where the e is, *ep points to this location */
		*ep = '\0';

		/* skip next char, */
		++ep;

		/* copy resulting string to return string */
		strcat( buf, ep );
	}
	return buf;
}
#endif

/**********************************************************/
void PrintE82( FILE* ioOUT, double value )
{
	double frac , xlog , xfloor , tvalue;
	int iExp;

	DEBUG_ENTRY( "PrintE82()" );

	if( value < 0 )
	{
		fprintf(ioOUT,"********");
	}
	else if( value <= DBL_MIN )
	{
		fprintf(ioOUT,"0.00E+00");
	}
	else
	{
		/* round number off for 8.2 format (not needed since can't be negative) */
		tvalue = value;
		xlog = log10( tvalue );
		xfloor = floor(xlog);
		/* this is now the fractional part */
		if (xfloor < 0.)
			frac = tvalue*exp10(-xfloor);
		else
			frac = (10.*tvalue)*exp10(-(xfloor+1.));
		/*this is the possibly signed exponential part */
		iExp = (int)xfloor;
		if( frac>9.9945 )
		{
			frac /= 10.;
			iExp += 1;
		}
		/* print the fractional part*/
		fprintf(ioOUT,"%.2f",frac);
		/* E for exponent */
		fprintf(ioOUT,"E");
		/* if positive throw a + sign*/
		if(iExp>=0 )
		{
			fprintf(ioOUT,"+");
		}
		fprintf(ioOUT,"%.2d",iExp);
	}
	return;
}
/*
 *==============================================================================
 */
void PrintE71( FILE* ioOUT, double value )
{
	double frac , xlog , xfloor , tvalue;
	int iExp;

	DEBUG_ENTRY( "PrintE71()" );

	if( value < 0 )
	{
		fprintf(ioOUT,"*******");
	}
	else if( value <= DBL_MIN )
	{
		fprintf(ioOUT,"0.0E+00");
	}
	else
	{
		/* round number off for 8.2 format (not needed since can't be negative) */
		tvalue = value;
		xlog = log10( tvalue );
		xfloor = floor(xlog);
		/* this is now the fractional part */
		if (xfloor < 0.)
			frac = tvalue*exp10(-xfloor);
		else
			frac = (10.*tvalue)*exp10(-(xfloor+1.));
		/*this is the possibly signed exponential part */
		iExp = (int)xfloor;
		if( frac>9.9945 )
		{
			frac /= 10.;
			iExp += 1;
		}
		/* print the fractional part*/
		fprintf(ioOUT,"%.1f",frac);
		/* E for exponent */
		fprintf(ioOUT,"E");
		/* if positive throw a + sign*/
		if(iExp>=0 )
		{
			fprintf(ioOUT,"+");
		}
		fprintf(ioOUT,"%.2d",iExp);
	}
	return;
}

/*
 *==============================================================================
 */
void PrintE93( FILE* ioOUT, double value )
{
	double frac , xlog , xfloor, tvalue;
	int iExp;

	DEBUG_ENTRY( "PrintE93()" );

	if( value < 0 )
	{
		fprintf(ioOUT,"*********");
	}
	else if( value <= DBL_MIN )
	{
		fprintf(ioOUT,"0.000E+00");
	}
	else
	{
		/* round number off for 9.3 format, neg numb not possible */
		tvalue = value;
		xlog = log10( tvalue );
		xfloor = floor(xlog);
		/* this is now the fractional part */
		if (xfloor < 0.)
			frac = tvalue*exp10(-xfloor);
		else
			frac = (10.*tvalue)*exp10(-(xfloor+1.));
		/*this is the possibly signed exponential part */
		iExp = (int)xfloor;
		if( frac>9.99949 )
		{
			frac /= 10.;
			iExp += 1;
		}
		/* print the fractional part*/
		fprintf(ioOUT,"%5.3f",frac);
		/* E for exponent */
		fprintf(ioOUT,"E");
		/* if positive throw a + sign*/
		if(iExp>=0 )
		{
			fprintf(ioOUT,"+");
		}
		fprintf(ioOUT,"%.2d",iExp);
	}
	return;
}

/*TotalInsanity general error handler for something that cannot happen */
NORETURN void TotalInsanity(void)
{
	DEBUG_ENTRY( "TotalInsanity()" );

	/* something that cannot happen, happened,
	 * if this message is triggered, simply place a breakpoint here
	 * and debug the error */
	fprintf( ioQQQ, " Something that cannot happen, has happened.\n" );
	fprintf( ioQQQ, " This is TotalInsanity, I live in %s.\n", __FILE__ );
	ShowMe();

	cdEXIT(EXIT_FAILURE);
}

/*BadRead general error handler for trying to read data, but failing */
NORETURN void BadRead(void)
{
	DEBUG_ENTRY( "BadRead()" );

	/* read failed */
	fprintf( ioQQQ, " A read of internal input data has failed.\n" );
	fprintf( ioQQQ, " This is BadRead, I live in %s.\n", __FILE__ );
	ShowMe();

	cdEXIT(EXIT_FAILURE);
}

/*sexp safe exponential function */
sys_float sexp(sys_float x)
{
	sys_float sexp_v;

	DEBUG_ENTRY( "sexp()" );

	/* SEXP_LIMIT is 84 in cddefines.h */
	if( x < SEXP_LIMIT )
	{
		sexp_v = exp(-x);
	}
	else
	{
		sexp_v = 0.f;
	}
	return sexp_v;
}

/*sexp safe exponential function */
double sexp(double x)
{
	double sexp_v;

	DEBUG_ENTRY( "sexp()" );

	/* SEXP_LIMIT is 84 in cddefines.h */
	if( x < SEXP_LIMIT )
	{
		sexp_v = exp(-x);
	}
	else
	{
		sexp_v = 0.;
	}
	return sexp_v;
}


/*dsexp safe exponential function for doubles */
double dsexp(double x)
{
	double dsexp_v;

	DEBUG_ENTRY( "dsexp()" );

	if( x < DSEXP_LIMIT )
	{
		dsexp_v = exp(-x);
	}
	else
	{
		dsexp_v = 0.;
	}
	return dsexp_v;
}

/*TestCode set flag saying that test code is in place 
 * prototype in cddefines.h */
void TestCode(void)
{
	DEBUG_ENTRY( "TestCode()" );

	/* called if test code is in place */
	lgTestCodeCalled = true;
	return;
}

/*broken set flag saying that the code is broken, */
void broken(void)
{
	DEBUG_ENTRY( "broken()" );

	broke.lgBroke = true;
	return;
}

/*fixit say that code needs to be fixed */
void fixit_base(const char* func,
					 const char* file, 
					 int         line, 
					 const char* reason
	)
{
	DEBUG_ENTRY( "fixit_base()" );

	broke.lgFixit = true;
	ostringstream oss;
	oss << "-- At " << file << ":"<< line << " in " << func <<"()\n";
	oss << reason;
	FixitList::Inst().list.push_back(oss.str());
	return;
}

/*CodeReview placed next to code that needs to be checked */
void CodeReview(void)
{
	DEBUG_ENTRY( "CodeReview()" );

	broke.lgCheckit = true;
	return;
}

/** dprintf -- version of fprintf which prepends DEBUG */
int dprintf(FILE *fp, const char *format, ...)
{
	va_list ap;
	int i1, i2;

	DEBUG_ENTRY( "dprintf()" );
	va_start(ap,format);
	i1 = fprintf(fp,"DEBUG ");
	if (i1 >= 0)
		i2 = vfprintf(fp,format,ap);
	else
		i2 = 0;
	if (i2 < 0)
		i1 = 0;
	va_end(ap);

	return i1+i2;
}

int fprintf (const Output& stream, const char *format, ...)
{
	va_list ap;
	va_start(ap,format);
	int i = vfprintf(stream.fptr(), format, ap);
	va_end(ap);
	return i;
}

int dprintf(const Output& stream, const char *format, ...)
{
	DEBUG_ENTRY( "dprintf()" );
	va_list ap;
	va_start(ap,format);
	int i1 = fprintf(stream.fptr(),"DEBUG ");
	int i2;
	if (i1 >= 0)
		i2 = vfprintf(stream.fptr(),format,ap);
	else
		i2 = 0;
	if (i2 < 0)
		i1 = 0;
	va_end(ap);

	return i1+i2;
}



/* dbg_printf is a debug print routine that was provided by Peter Teuben,
 * as a component from his NEMO package.  It offers run-time specification
 * of the level of debugging */
int dbg_printf(int debug, const char *fmt, ...)
{
	va_list ap;
	int i=0;

	DEBUG_ENTRY( "dbg_printf()" );

	/* print this debug message? (debug_level not currently used)*/
	if(debug <= trace.debug_level) 
	{		
		va_start(ap, fmt);	

		i = vfprintf(ioQQQ, fmt, ap);
		/* drain ioQQQ */
		fflush(ioQQQ);
		va_end(ap);
	}
	return i;
}


/*qg32 32 point Gaussian quadrature, originally given to Gary F by Jim Lattimer */
double qg32(
	double xl, /*lower limit to integration range*/
	double xu, /*upper limit to integration range*/
	/*following is the pointer to the function that will be evaluated*/
	double (*fct)(double) )
{
	double a = 0.5*(xu+xl), 
	  b = xu-xl, 
	  y = 0.;

	DEBUG_ENTRY( "qg32()" );

	/********************************************************************************
	 *                                                                              *
	 *  32-point Gaussian quadrature                                                *
	 *  xl  : the lower limit of integration                                        *
	 *  xu  : the upper limit                                                       *
	 *  fct : the (external) function                                               *
	 *  returns the value of the integral                                           *
	 *                                                                              *
	 * simple call to integrate sine from 0 to pi                                   *
	 * double agn = qg32( 0., 3.141592654 ,  sin );                                 *
	 *                                                                              *
	 *******************************************************************************/

	double weights[16] = {
		.35093050047350483e-2, .81371973654528350e-2, .12696032654631030e-1, .17136931456510717e-1,
		.21417949011113340e-1, .25499029631188088e-1, .29342046739267774e-1, .32911111388180923e-1,
		.36172897054424253e-1, .39096947893535153e-1, .41655962113473378e-1, .43826046502201906e-1,
		.45586939347881942e-1, .46922199540402283e-1, .47819360039637430e-1, .48270044257363900e-1};

	double c[16] = {
		.498631930924740780, .49280575577263417, .4823811277937532200, .46745303796886984000,
		.448160577883026060, .42468380686628499, .3972418979839712000, .36609105937014484000,
		.331522133465107600, .29385787862038116, .2534499544661147000, .21067563806531767000,
		.165934301141063820, .11964368112606854, .7223598079139825e-1, .24153832843869158e-1};

	for( int i=0; i<16; i++)
	{
		y += b * weights[i] * ((*fct)(a+b*c[i]) + (*fct)(a-b*c[i]));
	}

	/* the answer */
	return y;
}

/*spsort netlib routine to sort array returning sorted indices */
void spsort(
	  /* input array to be sorted */
	  realnum x[], 
	  /* number of values in x */
	  long int n, 
	  /* permutation output array */
	  long int iperm[], 
	  /* flag saying what to do - 1 sorts into increasing order, not changing
	   * the original vector, -1 sorts into decreasing order. 2, -2 change vector */
	  int kflag, 
	  /* error condition, should be 0 */
	  int *ier)
{
	/*
	 ****BEGIN PROLOGUE  SPSORT
	 ****PURPOSE  Return the permutation vector generated by sorting a given
	 *            array and, optionally, rearrange the elements of the array.
	 *            The array may be sorted in increasing or decreasing order.
	 *            A slightly modified quicksort algorithm is used.
	 ****LIBRARY   SLATEC
	 ****CATEGORY  N6A1B, N6A2B
	 ****TYPE      SINGLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
	 ****KEY WORDS NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
	 ****AUTHOR  Jones, R. E., (SNLA)
	 *           Rhoads, G. S., (NBS)
	 *           Wisniewski, J. A., (SNLA)
	 ****DESCRIPTION
	 *
	 *   SPSORT returns the permutation vector IPERM generated by sorting
	 *   the array X and, optionally, rearranges the values in X.  X may
	 *   be sorted in increasing or decreasing order.  A slightly modified
	 *   quicksort algorithm is used.
	 *
	 *   IPERM is such that X(IPERM(I)) is the Ith value in the rearrangement
	 *   of X.  IPERM may be applied to another array by calling IPPERM,
	 *   SPPERM, DPPERM or HPPERM.
	 *
	 *   The main difference between SPSORT and its active sorting equivalent
	 *   SSORT is that the data are referenced indirectly rather than
	 *   directly.  Therefore, SPSORT should require approximately twice as
	 *   long to execute as SSORT.  However, SPSORT is more general.
	 *
	 *   Description of Parameters
	 *      X - input/output -- real array of values to be sorted.
	 *          If ABS(KFLAG) = 2, then the values in X will be
	 *          rearranged on output; otherwise, they are unchanged.
	 *      N - input -- number of values in array X to be sorted.
	 *      IPERM - output -- permutation array such that IPERM(I) is the
	 *              index of the value in the original order of the
	 *              X array that is in the Ith location in the sorted
	 *              order.
	 *      KFLAG - input -- control parameter:
	 *            =  2  means return the permutation vector resulting from
	 *                  sorting X in increasing order and sort X also.
	 *            =  1  means return the permutation vector resulting from
	 *                  sorting X in increasing order and do not sort X.
	 *            = -1  means return the permutation vector resulting from
	 *                  sorting X in decreasing order and do not sort X.
	 *            = -2  means return the permutation vector resulting from
	 *                  sorting X in decreasing order and sort X also.
	 *      IER - output -- error indicator:
	 *          =  0  if no error,
	 *          =  1  if N is zero or negative,
	 *          =  2  if KFLAG is not 2, 1, -1, or -2.
	 ****REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
	 *                 for sorting with minimal storage, Communications of
	 *                 the ACM, 12, 3 (1969), pp. 185-187.
	 ****ROUTINES CALLED  XERMSG
	 ****REVISION HISTORY  (YYMMDD)
	 *   761101  DATE WRITTEN
	 *   761118  Modified by John A. Wisniewski to use the Singleton
	 *           quicksort algorithm.
	 *   870423  Modified by Gregory S. Rhoads for passive sorting with the
	 *           option for the rearrangement of the original data.
	 *   890620  Algorithm for rearranging the data vector corrected by R.
	 *           Boisvert.
	 *   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
	 *   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
	 *   920507  Modified by M. McClain to revise prologue text.
	 *   920818  Declarations section rebuilt and code restructured to use
	 *           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
	 ****END PROLOGUE  SPSORT
	 *     .. Scalar Arguments ..
	 */
	long int i, 
	  ij, 
	  il[21], 
	  indx, 
	  indx0, 
	  istrt, 
	  istrt_, 
	  iu[21], 
	  j, 
	  k, 
	  kk, 
	  l, 
	  lm, 
	  lmt, 
	  m, 
	  nn;
	realnum r, 
	  ttemp;

	DEBUG_ENTRY( "spsort()" );

	/*     .. Array Arguments .. */
	/*     .. Local Scalars .. */
	/*     .. Local Arrays .. */
	/*     .. External Subroutines .. */
	/*     .. Intrinsic Functions .. */
	/****FIRST EXECUTABLE STATEMENT  SPSORT */
	*ier = 0;
	nn = n;
	if( nn < 1 )
	{
		*ier = 1;
		return;
	}
	else
	{
		kk = labs(kflag);
		if( kk != 1 && kk != 2 )
		{
			*ier = 2;
			return;
		}
		else
		{

			/* Initialize permutation vector to index on f scale
			 * */
			for( i=0; i < nn; i++ )
			{
				iperm[i] = i+1;
			}

			/* Return if only one value is to be sorted */
			if( nn == 1 )
			{ 
				--iperm[0];
				return;
			}

			/* Alter array X to get decreasing order if needed */
			if( kflag <= -1 )
			{
				for( i=0; i < nn; i++ )
				{
					x[i] = -x[i];
				}
			}

			/* Sort X only */
			m = 1;
			i = 1;
			j = nn;
			r = .375e0;
		}
	}

	while( true )
	{
		if( i == j )
			goto L_80;
		if( r <= 0.5898437e0 )
		{
			r += 3.90625e-2;
		}
		else
		{
			r -= 0.21875e0;
		}

L_40:
		k = i;

		/*     Select a central element of the array and save it in location L
		 * */
		ij = i + (long)((j-i)*r);
		lm = iperm[ij-1];

		/*     If first element of array is greater than LM, interchange with LM
		 * */
		if( x[iperm[i-1]-1] > x[lm-1] )
		{
			iperm[ij-1] = iperm[i-1];
			iperm[i-1] = lm;
			lm = iperm[ij-1];
		}
		l = j;

		/*     If last element of array is less than LM, interchange with LM
		 * */
		if( x[iperm[j-1]-1] < x[lm-1] )
		{
			iperm[ij-1] = iperm[j-1];
			iperm[j-1] = lm;
			lm = iperm[ij-1];

			/*        If first element of array is greater than LM, interchange
			 *        with LM
			 * */
			if( x[iperm[i-1]-1] > x[lm-1] )
			{
				iperm[ij-1] = iperm[i-1];
				iperm[i-1] = lm;
				lm = iperm[ij-1];
			}
		}

		/* Find an element in the second half of the array which is smaller
		 * than LM */
		while( true )
		{
			l -= 1;
			if( x[iperm[l-1]-1] <= x[lm-1] )
			{

				/* Find an element in the first half of the array which is greater
				 * than LM */
				while( true )
				{
					k += 1;
					if( x[iperm[k-1]-1] >= x[lm-1] )
						break;
				}

				/* Interchange these elements */
				if( k > l )
					break;
				lmt = iperm[l-1];
				iperm[l-1] = iperm[k-1];
				iperm[k-1] = lmt;
			}
		}

		/* Save upper and lower subscripts of the array yet to be sorted */
		if( l - i > j - k )
		{
			il[m-1] = i;
			iu[m-1] = l;
			i = k;
			m += 1;
		}
		else
		{
			il[m-1] = k;
			iu[m-1] = j;
			j = l;
			m += 1;
		}

L_90:
		if( j - i >= 1 )
			goto L_40;
		if( i == 1 )
			continue;
		i -= 1;

		while( true )
		{
			i += 1;
			if( i == j )
				break;
			lm = iperm[i];
			if( x[iperm[i-1]-1] > x[lm-1] )
			{
				k = i;

				while( true )
				{
					iperm[k] = iperm[k-1];
					k -= 1;

					if( x[lm-1] >= x[iperm[k-1]-1] )
						break;
				}
				iperm[k] = lm;
			}
		}

		/* Begin again on another portion of the unsorted array */
L_80:
		m -= 1;
		if( m == 0 )
			break;
		/*lint -e644 not explicitly initialized */
		i = il[m-1];
		j = iu[m-1];
		/*lint +e644 not explicitly initialized */
		goto L_90;
	}

	/* Clean up */
	if( kflag <= -1 )
	{
		for( i=0; i < nn; i++ )
		{
			x[i] = -x[i];
		}
	}

	/* Rearrange the values of X if desired */
	if( kk == 2 )
	{

		/* Use the IPERM vector as a flag.
		 * If IPERM(I) < 0, then the I-th value is in correct location */
		for( istrt=1; istrt <= nn; istrt++ )
		{
			istrt_ = istrt - 1;
			if( iperm[istrt_] >= 0 )
			{
				indx = istrt;
				indx0 = indx;
				ttemp = x[istrt_];
				while( iperm[indx-1] > 0 )
				{
					x[indx-1] = x[iperm[indx-1]-1];
					indx0 = indx;
					iperm[indx-1] = -iperm[indx-1];
					indx = labs(iperm[indx-1]);
				}
				x[indx0-1] = ttemp;
			}
		}

		/* Revert the signs of the IPERM values */
		for( i=0; i < nn; i++ )
		{
			iperm[i] = -iperm[i];
		}
	}

	for( i=0; i < nn; i++ )
	{
		--iperm[i];
	}
	return;
}

/*MyMalloc wrapper for malloc().  Returns a good pointer or dies. 
 * memory is filled with NaN
 * >>chng 05 dec 14, do not set to NaN since tricks debugger 
 * routines within code do not call this or malloc, but rather MALLOC
 * which is resolved into MyMalloc or malloc depending on whether 
 * NDEBUG is set by the compiler to indicate "not debugging",
 * in typical negative C style */
STATIC void MyMalloc_base(void* ptr,
			  size_t size, /*use same type as library function MALLOC*/ 
			  const char *chFile, 
			  int line)
{
	DEBUG_ENTRY( "MyMalloc_base()" );

	/* debug branch for printing malloc args */
	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC)
		{
			static long int kount=0, nTot=0;
			nTot += (long)size;
			fprintf(ioQQQ,"%li\t%li\t%li\n", 
				kount , 
				(long)size , 
				nTot );
			++kount;
		}
	}

	if( ptr == NULL )
	{
		fprintf(ioQQQ,"DISASTER MyMalloc could not allocate %lu bytes.  Exit in MyMalloc.",
			(unsigned long)size );
		fprintf(ioQQQ,"MyMalloc called from file %s at line %i.\n",
			chFile , line );

		if( struc.nzlim>2000 )
			fprintf(ioQQQ,"This may have been caused by the large number of zones."
			" %li zones were requested.  Is this many zones really necessary?\n",
			struc.nzlim );

		cdEXIT(EXIT_FAILURE);
	}

	/* flag -DNOINIT will turn off this initialization which can fool valgrind/purify */
#	if !defined(NDEBUG) && !defined(NOINIT)

	size_t nFloat = size/4;
	size_t nDouble = size/8;
	sys_float *fptr = static_cast<sys_float*>(ptr);
	double *dptr = static_cast<double*>(ptr);

	/* >>chng 04 feb 03, fill memory with invalid numbers, PvH */
	/* on IA32/AMD64 processors this will generate NaN's for both float and double;
	 * on most other (modern) architectures it is likely to do the same... */
	/* >>chng 05 dec 14, change code to generate signaling NaN's for most cases (but not all!) */
	if( size == nDouble*8 )
	{
		/* this could be an array of doubles as well as floats -> we will hedge our bets
		 * we will fill the array with a pattern that is interpreted as all signaling
		 * NaN's for doubles, and alternating signaling and quiet NaN's for floats:
		 * byte offset:  0      4      8     12     16
		 * double        | SNaN        | SNan        |
		 * float         | SNaN | QNaN | SNan | QNaN |  (little-endian, e.g. Intel, AMD, alpha)
		 * float         | QNaN | SNaN | QNan | SNaN |  (big-endian, e.g. Sparc, PowerPC, MIPS) */
		set_NaN( dptr, (long)nDouble );
	}
	else if( size == nFloat*4 )
	{
		/* this could be an arrays of floats, but not doubles -> init to all float SNaN */
		set_NaN( fptr, (long)nFloat );
	}
	else
	{
		memset( ptr, 0xff, size );
	}

#	endif /* !defined(NDEBUG) && !defined(NOINIT) */
}

void *MyMalloc(size_t size,
	       const char *chFile, 
	       int line)
{
	ASSERT( size > 0 );

	void *ptr = malloc(size);
	MyMalloc_base(ptr, size, chFile, line);
	return ptr;
}

/* function to facilitate addressing opacity array */
double csphot(
	/* INU is array index pointing to frequency where opacity is to be evaluated
	 * on f not c scale */
	long int inu, 
	/* ITHR is pointer to threshold*/
	long int ithr, 
	/* IOFSET is offset as defined in opac0*/
	long int iofset)
{
	double csphot_v;

	DEBUG_ENTRY( "csphot()" );

	csphot_v = opac.OpacStack[inu-ithr+iofset-1];
	return csphot_v;
}

/*RandGauss normal Gaussian random number generator
 * the user must call srand to set the seed before using this routine.

 * the random numbers will have a mean of xMean and a standard deviation of s

 * The convention is for srand to be called when the command setting 
 * the noise is parsed 

 * for very small dispersion there are no issues, but when the dispersion becomes
 * large the routine will find negative values - so most often used in this case
 * to find dispersion in log space 
 * this routine will return a normal Gaussian - must be careful in how this is
 * used when adding noise to physical quantity */
/*
NB - following from Ryan Porter:
I discovered that I unintentionally created an antisymmetric skew in my 
Monte Carlo.  RandGauss is symmetric in log space, which means it is not 
symmetric in linear space.  But to get the right standard deviation you 
have to take 10^x, where x is the return from RandGauss.  The problem is 
10^x will happen less frequently than 10^-x, so without realizing it, the 
average "tweak" to every piece of atomic data in my Monte Carlo run was 
not 1.0 but something greater than 1.0, causing every single line to have 
an average Monte Carlo emissivity greater than the regular value.  Any place
that takes 10^RandGauss() needs to be adjusted if what is intended is +/- x. */
double RandGauss(
	/* mean value */
	double xMean, 
	/*standard deviation s */
	double s )
{
	double  x1, x2, w, yy1;
	static double yy2=-BIGDOUBLE;
	static int use_last = false;

	DEBUG_ENTRY( "RandGauss()" );

	if( use_last )
	{
		yy1 = yy2;
		use_last = false;
	}
	else
	{
		do {
			x1 = 2.*genrand_real3() - 1.;
			x2 = 2.*genrand_real3() - 1.;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt((-2.0*log(w))/w);
		yy1 = x1 * w;
		yy2 = x2 * w;
		use_last = true;
	}
	return xMean + yy1 * s;
}

/* MyGaussRand takes as input a percent uncertainty less than 50%
 * (expressed as 0.5). The routine then assumes this input variable represents one
 * standard deviation about a mean of unity, and returns a random number within
 * that range.  A hard cutoff is imposed at two standard deviations, which 
 * eliminates roughly 5% of the normal distribution.  In other words, the routine
 * returns a number in a normal distribution with standard deviation equal to
 * the input.  The number will be between 1-3*stdev and 1+3*stdev. */
double MyGaussRand( double PctUncertainty )
{
	double StdDev;
	double result;

	DEBUG_ENTRY( "MyGaussRand()" );

	ASSERT( PctUncertainty < 0.5 );
	/* We want this "percent uncertainty" to represent one standard deviation */
	StdDev = PctUncertainty;

	do
	{
		/*result = exp10(  RandGauss( 0., logStdDev ) );*/
		result = 1.+RandGauss( 0., StdDev );
	}
	/* only allow values that are within 3 standard deviations */
	while( (result < 1.-3.*PctUncertainty) || (result > 1.+3.*PctUncertainty) );

	ASSERT( result>0. && result<2. );
	return result;
}

/*plankf evaluate Planck function for any cell at current electron temperature */
double plankf(long int ip)
{
	double plankf_v;

	DEBUG_ENTRY( "plankf()" );

	/* evaluate Planck function
	 * argument is pointer to cell energy in ANU
	 * return photon flux for cell IP */
	if( rfield.ContBoltz[ip] <= 0. )
	{
		plankf_v = 1e-35;
	}
	else
	{
		plankf_v = 6.991e-21*POW2(FR1RYD*rfield.anu(ip))/
			(1./rfield.ContBoltz[ip] - 1.)*FR1RYD*4.;
	}
	return plankf_v;
}

void CloudyPrintReference()
{
	fstream io;
	string line;
	open_data( io, "citation_cloudy.txt", mode_r );
	while( SafeGetline( io, line ) )
	{
		if( line[0] == '#' )
			continue;
		// replace XXXX with actual version number
		size_t p = line.find( "XXXX" );
		if( p != string::npos )
			line.replace( p, 4, t_version::Inst().chVersion );
		fprintf( ioQQQ, "%s\n", line.c_str() );
	}
}

void DatabasePrintReference()
{
	fstream io;
	string line;
	open_data( io, "citation_data.txt", mode_r );
	while( SafeGetline( io, line ) )
	{
		if( line[0] == '#' )
			continue;
		// replace XXXX with actual version number
		size_t p = line.find( "XXXX" );
		if( p != string::npos )
			line.replace( p, 4, atmdat.chVersion );
		fprintf( ioQQQ, "%s\n", line.c_str() );
	}
}

// this routine was taken from
// http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
// it is copyrighted by a creative commons license
// http://creativecommons.org/licenses/by-sa/3.0/
//
// safe version of getline() that correctly handles all types of EOL lf, crlf and cr...
// it has been modified such that it does not produce a spurious empty line at the end of a file
// this way it is compatible with the standard getline() (at least with g++/linux).
istream& SafeGetline(istream& is, string& t)
{
	t.clear();

	// The characters in the stream are read one-by-one using a streambuf.
	// That is faster than reading them one-by-one using the istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	istream::sentry se(is, true);
	streambuf* sb = is.rdbuf();

	while( true )
	{
		int c = sb->sbumpc();
		switch (c)
		{
		case '\n':
			if( sb->sgetc() == EOF )
				is.setstate(ios::eofbit);
			return is;
		case '\r':
			if( sb->sgetc() == '\n' )
				sb->sbumpc();
			if( sb->sgetc() == EOF )
				is.setstate(ios::eofbit);
			return is;
		case EOF:
			// Also handle the case when the last line has no line ending
			is.setstate(ios::eofbit);
			return is;
		default:
			t += (char)c;
		}
	}
}
