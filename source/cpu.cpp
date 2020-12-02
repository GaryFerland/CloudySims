/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/**\file cpu.cpp implement hardware dependent definitions */

#include "cdstd.h"
#include <locale.h>

#if defined(__ia64) && defined(__INTEL_COMPILER)
extern "C" unsigned long fpgetmask();
extern "C" void fpsetmask(unsigned long);
#endif	

#if defined(__sun) || defined(__sgi)
#include <ieeefp.h>
#if defined(HAVE_SUNMATH) || defined(FLUSH_DENORM_TO_ZERO)
#include <sunmath.h>
#endif
#endif

#if defined(__alpha) && defined(__linux__) && defined(__GNUC__)
#define __USE_GNU
#include <fenv.h>
#endif

#if defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x5140
#include <fenv.h>
#endif

#if defined(__unix) || defined(__APPLE__)
#include <unistd.h>
#endif

#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

/* the redefinition of float in cddefines.h can cause problems in system headers
 * hence these includes MUST come after the system header includes above */
#include "cddefines.h"
#include "path.h"
#include "trace.h"
#include "service.h"
#include "thirdparty.h"
#include "version.h"
#include "prt.h"

STATIC NORETURN void AbortErrorMessage( const char* fname, const vector<string>& PathList, access_scheme scheme );
STATIC string check_mult_path( const char* fname, const vector<string>& PathList, access_scheme scheme, bool lgRead );

// Use Schwartz/nifty counter to ensure that global policy class
// is set up before other globals/statics, and deleted last.
t_cpu_i *t_cpu::m_i;
static int cpu_count = 0;
t_cpu::t_cpu()
{
	if (0 == cpu_count++)
	{
		m_i = new t_cpu_i;
	}
}
t_cpu::~t_cpu()
{
	if (0 == --cpu_count)
	{
		delete m_i;
	}
}

/* NB NB - this constructor needs to be called before any of the user code is executed !! */
t_cpu_i::t_cpu_i()
{
	DEBUG_ENTRY( "t_cpu_i()" );

	// set up signal handlers so that we can control what happens...
	set_signal_handlers();

	p_exit_status.resize( ES_TOP, "--undefined--" );
	p_exit_status[ES_SUCCESS]             = "ok";
	p_exit_status[ES_FAILURE]             = "early termination";
	p_exit_status[ES_WARNINGS]            = "warnings";
	p_exit_status[ES_BOTCHES]             = "botched monitors";
	p_exit_status[ES_CLOUDY_ABORT]        = "cloudy abort";
	p_exit_status[ES_BAD_ASSERT]          = "failed assert";
	p_exit_status[ES_BAD_ALLOC]           = "failed memory alloc";
	p_exit_status[ES_OUT_OF_RANGE]        = "array bound exceeded";
	p_exit_status[ES_DOMAIN_ERROR]        = "math domain error";
	p_exit_status[ES_USER_INTERRUPT]      = "user interrupt";
	p_exit_status[ES_TERMINATION_REQUEST] = "process killed";
	p_exit_status[ES_ILLEGAL_INSTRUCTION] = "illegal instruction";
	p_exit_status[ES_FP_EXCEPTION]        = "fp exception";
	p_exit_status[ES_SEGFAULT]            = "segmentation fault";
	p_exit_status[ES_BUS_ERROR]           = "bus error";
	p_exit_status[ES_UNKNOWN_SIGNAL]      = "unknown signal";
	p_exit_status[ES_UNKNOWN_EXCEPTION]   = "unknown exception";

	/* >>chng 05 dec 14, add test of endianness of the CPU, PvH */
	endian.c[0] = 0x12;
	endian.c[1] = 0x34;
	endian.c[2] = 0x56;
	endian.c[3] = 0x78;

	/* >>chng 05 dec 15, add signaling NaN for float and double to cpu struct, PvH */
	/* in C++ this should be replaced by numeric_limits<TYPE>::signaling_NaN() */
	if( sizeof(sys_float) == 4 )
	{
#		ifdef __mips
		/* definition of signaling and quiet NaN is reversed on MIPS */
		Float_SNaN_Value = 0xffffffff;
#		else
		if( big_endian() || little_endian() )
		{
			/* this should work on most modern CPU's */
			Float_SNaN_Value = 0xffbfffff;
		}
		else
		{
			/* this is an unusual CPU -> bit pattern for SNaN is unknown */
			Float_SNaN_Value = -1;
		}
#		endif
	}
	else
	{
		/* this is an unusual CPU -> bit pattern for SNaN is unknown */
		Float_SNaN_Value = -1;
	}

#	ifdef HAVE_INT64

	if( sizeof(double) == 8 )
	{
#		ifdef __mips
		/* definition of signaling and quiet NaN is reversed on MIPS */
		Double_SNaN_Value = INT64_LIT(0xffffffffffffffff);
#		else
		/* this should work on most modern CPU's */
		Double_SNaN_Value = INT64_LIT(0xfff7ffffffbfffff);
#		endif
	}
	else
	{
		/* this is an unusual CPU -> bit pattern for SNaN is unknown */
		Double_SNaN_Value = -1;
	}

#	else

	if( sizeof(double) == 8 )
	{
#		ifdef __mips
		/* definition of signaling and quiet NaN is reversed on MIPS */
		Double_SNaN_Value[0] = 0xffffffff;
		Double_SNaN_Value[1] = 0xffffffff;
#		else
		if( big_endian() )
		{
			/* this should work on most modern CPU's */
			Double_SNaN_Value[0] = 0xfff7ffff;
			Double_SNaN_Value[1] = 0xffbfffff;
		}
		else if( little_endian() )
		{
			/* this should work on most modern CPU's */
			Double_SNaN_Value[0] = 0xffbfffff;
			Double_SNaN_Value[1] = 0xfff7ffff;
		}
		else
		{
			/* this is an unusual CPU -> bit pattern for SNaN is unknown */
			Double_SNaN_Value[0] = -1;
			Double_SNaN_Value[1] = -1;
		}
#		endif
	}
	else
	{
		/* this is an unusual CPU -> bit pattern for SNaN is unknown */
		Double_SNaN_Value[0] = -1;
		Double_SNaN_Value[1] = -1;
	}

#	endif

	/* set FP environment to trap FP exceptions */
	enable_traps();

	ioStdin = stdin;
	ioQQQ = stdout;
	ioPrnErr = stderr;
	lgPrnErr = false;

	test_float = FLT_MIN;
	test_double = DBL_MIN;

	/* default is for failed asserts not to abort */
	p_lgAssertAbort = false;

	const char *str;

	/* determine the no. of CPUs on this machine; used by PHYMIR, grid command, .... */
#	if defined(_SC_NPROCESSORS_ONLN)  /* Linux, Sun Sparc, DEC Alpha, MacOS (OS releases >= 10.4) */
#if defined(__APPLE__) /* MacOS only use physical cores*/
	size_t sizeOfInt = sizeof(int);
	int physicalCores;
	sysctlbyname("hw.physicalcpu", &physicalCores, &sizeOfInt, NULL, 0);
	n_avail_CPU = int(physicalCores);
#else
	n_avail_CPU = sysconf(_SC_NPROCESSORS_ONLN);
#endif
#	elif defined(_SC_NPROC_ONLN)      /* SGI Iris */
	n_avail_CPU = sysconf(_SC_NPROC_ONLN);
#	elif defined(_SC_CRAY_NCPU)       /* Cray */
	n_avail_CPU = sysconf(_SC_CRAY_NCPU);
#	elif defined(_WIN32)              /* Microsoft Windows */
	str = getenv( "NUMBER_OF_PROCESSORS" );
	if( str != NULL )
	{
		int found = sscanf( str, "%ld", &n_avail_CPU );
		if( found != 1 )
			n_avail_CPU = 1;
	}
	else
	{
		n_avail_CPU = 1;
	}
#	elif defined(HW_AVAILCPU)         /* MacOS, BSD variants */
#if defined(__APPLE__) /* MacOS only use physical cores*/
	size_t sizeOfInt = sizeof(int);
	int physicalCores;
	sysctlbyname("hw.physicalcpu", &physicalCores, &sizeOfInt, NULL, 0);
	n_avail_CPU = int(physicalCores);
#else
	int mib[2];
	size_t len = sizeof(n_avail_CPU);
	mib[0] = CTL_HW;
	mib[1] = HW_AVAILCPU;  // alternatively, try HW_NCPU;
	sysctl(mib, 2, &n_avail_CPU, &len, NULL, 0);
	if( n_avail_CPU < 1 ) 
	{
		mib[1] = HW_NCPU;
		sysctl(mib, 2, &n_avail_CPU, &len, NULL, 0);
		if( n_avail_CPU < 1 )
			n_avail_CPU = 1;
	}
#endif
#	else                              
	/* Other systems, supply no. of CPUs on OPTIMIZE PHYMIR command line */
	n_avail_CPU = 1;
#	endif
	/* the constructor is run before MPI starts, so the rank is not available yet */
#	ifdef MPI_ENABLED
	p_lgMPI = true;
#	else
	p_lgMPI = false;
#	endif
	/* the default is for all ranks to cooperate on the same sim */
	p_lgMPISingleRankMode = false;
	n_rank = 0;

#	ifdef _WIN32	
	str = getenv( "COMPUTERNAME" );
#	else
	str = getenv( "HOSTNAME" );
#	endif

	if( str != NULL )
		strncpy( HostName, str, STDLEN-1 );
	else
		strncpy( HostName, "unknown", STDLEN-1 );
	HostName[STDLEN-1] = '\0';

	/* pick up the path from the environment, if set by user */
	const char *path = getenv( "CLOUDY_DATA_PATH" );

#	ifdef _WIN32
	string separator( ";" );
	p_chDirSeparator = '\\';
#	else
	string separator( ":" );
	p_chDirSeparator = '/';
#	endif

	// if the environment variable was not set, use the preprocessor symbol CLOUDY_DATA_PATH
	// this is normally defined in the Makefile, but can also be set in path.h (deprecated)
	string chSearchPathRaw = ( path != NULL ) ? string( path ) : string( CLOUDY_DATA_PATH );
	// the current working directory should be first
	chSearchPathRaw = "." + separator + chSearchPathRaw + separator;

#ifdef CLOUDY_ROOT
	// CLOUDY_ROOT is only defined in the Makefile
	string chCloudyRoot = string( CLOUDY_ROOT );

	// expand "+" to the default search path for this installation
	while( FindAndReplace( chSearchPathRaw, separator + "+" + separator,
			       separator + chCloudyRoot + p_chDirSeparator + "data" + separator ) ) {}
#endif

	Split( chSearchPathRaw, separator, chSearchPath, SPM_RELAX );

	for( vector<string>::iterator p=chSearchPath.begin(); p != chSearchPath.end(); ++p )
	{
#ifdef HAVE_REALPATH
		// clean up path
		char* ptr = realpath( p->c_str(), NULL );
		if( ptr != NULL )
		{
			*p = ptr;
			free( ptr );
		}
#endif

		/* get last valid char */
		char chEnd = *p->rbegin();

		/* make sure path ends with directory separator */
		if( chEnd != p_chDirSeparator )
			*p += p_chDirSeparator;

#if 0
		// enable this code block once the search path is implemented in its final form
		// test it by running a sim inside the data directory that e.g. uses grains
		// it should not produce any warnings about multiple grain files along the path

		// check if this path component is already present
		if( find( chSearchPath.begin(), p, *p ) != p )
			chSearchPath.erase(p);
#endif
	}

	nFileDone = 0;
	nMD5Mismatch = 0;

	// avoid vagueries of different locales on different machines...
	setlocale(LC_ALL, "POSIX");

	getMD5sums( "md5datafiles.dat" );
}

void t_cpu_i::enable_traps() const
{
	/* >>chng 01 aug 07, added code to circumvent math library bug with g++ on
	 * alpha-linux machines, see bug report 51072 on http://bugzilla.redhat.com, PvH */
	/* >>chng 01 apr 17, added code for Solaris and SGI operating systems, PvH */
	/* this routine contains no code for alphas or crays, they do not need
	 * special code to enable FP exceptions since they are enabled by default */

	/* there is no command line option on MS Visual Studio to force crash */
#	if defined(_MSC_VER)
	volatile unsigned int NewMask;

	/* | is a bitwise inclusive or, turns on bits
	 * 0|0 = 0
	 * 0|1 = 1|0 = 1|1 = 1 */
	NewMask = _EM_ZERODIVIDE | _EM_OVERFLOW | _EM_INVALID;
	/* ~ is the unary bitwise complement - all bits flip */ 
	NewMask = ~NewMask;
	_controlfp( NewMask , _MCW_EM );

	/* this is the code for Linux PC (but not Linux alpha) to force crash */
	/* >>chng 04 apr 26, added support for AMD64, enable FPE traps for SSE/SSE2, PvH */
	/* >>chng 06 aug 12, added support for Apple MacOSX, and hopefully also Solaris x86, PvH */
	// do not enable trapping FPEs for Clang compiler since it is not supported in v3.4 and later
#	elif defined(__GNUC__) && ( defined(__i386) || defined(__amd64) ) && !defined(__clang__)
	volatile unsigned int Old_Mask, New_Mask;
#	if defined(__SSE__) || defined(__SSE2__)
	volatile unsigned int SSE_Mask;
#	endif

#	define _FPU_MASK_IM  0x01  /* Invalid          */
#	define _FPU_MASK_DM  0x02  /* Denormalized     */
#	define _FPU_MASK_ZM  0x04  /* Division-by-zero */
#	define _FPU_MASK_OM  0x08  /* Overflow         */
#	define _FPU_MASK_UM  0x10  /* Underflow        */
#	define _FPU_MASK_PM  0x20  /* Inexact          */

	/*  | is a bitwise inclusive or, turns on bits     */
	/* 0|0 = 0                                         */
	/* 0|1 = 1|0 = 1|1 = 1                             */

	/*  ~ is the unary bitwise complement - all bits flip */

	/* this enables FPE traps for regular i387 FP instructions */

	volatile unsigned int UnMask = ~((unsigned int)( _FPU_MASK_ZM | _FPU_MASK_IM  | _FPU_MASK_OM ));

	__asm__ volatile("fnstcw %0" : "=m" (*&Old_Mask));

	New_Mask = Old_Mask & UnMask;

	__asm__ volatile("fldcw %0" : : "m" (*&New_Mask));

#	if defined(__SSE__) || defined(__SSE2__)

#	ifndef USE_DENORM
	/* using this causes denormalized numbers to be flushed to zero,
	 * which will speed up the code on Pentium 4 processors */
	SSE_Mask = 0x9900;
#	else
	/* this version allows denormalized numbers to be retained */
	SSE_Mask = 0x1900;
#	endif

	/* this enables FPE traps for SSE/SSE2 instructions */

	__asm__ volatile( "ldmxcsr %0" : : "m" (*&SSE_Mask) );

#	endif

	/* this is for IA64 systems running g++ or icc (e.g. SGI, HP, ...) */
#	elif defined(__ia64)

#	define FPSR_TRAP_VD     (1 << 0)        /* invalid op trap disabled */
#	define FPSR_TRAP_DD     (1 << 1)        /* denormal trap disabled */
#	define FPSR_TRAP_ZD     (1 << 2)        /* zero-divide trap disabled */
#	define FPSR_TRAP_OD     (1 << 3)        /* overflow trap disabled */
#	define FPSR_TRAP_UD     (1 << 4)        /* underflow trap disabled */
#	define FPSR_TRAP_ID     (1 << 5)        /* inexact trap disabled */

#	define FPSR_SF0_FTZ     (1 << 6)        /* flush denormalized numbers to zero */

#	if defined(__GNUC_EXCL__)
	/* __asm__ instructions are not supported by icc as of v9.0 */
#	define _IA64_REG_AR_FPSR    40

#	define ia64_getreg( regnum )       __asm__ volatile( "mov %0=ar%1" : "=r" (fpsr) : "i"(regnum) )
#	define ia64_setreg( regnum, val )  __asm__ volatile( "mov ar%0=%1" :: "i" (regnum), "r"(val): "memory" )
#	define ia64_serialize              __asm__ volatile( "srlz.i" );

	volatile unsigned long fpsr, flags = FPSR_TRAP_VD | FPSR_TRAP_ZD | FPSR_TRAP_OD;

	ia64_getreg( _IA64_REG_AR_FPSR );
	fpsr &= ~flags;
#	if defined(FLUSH_DENORM_TO_ZERO)
	fpsr |= FPSR_SF0_FTZ;
#	endif
	ia64_setreg( _IA64_REG_AR_FPSR, fpsr );
	/* this prevents RAW and WAW dependency violations in case this ever gets inlined... */
	ia64_serialize;

#	elif defined(__INTEL_COMPILER)
	/* this is for icc on IA64 SGI machines */
	unsigned long fpsr = fpgetmask();
	fpsr |= FPSR_TRAP_VD | FPSR_TRAP_ZD | FPSR_TRAP_OD;
	fpsetmask( fpsr );
#	endif /* defined(__GNUC_EXCL__) */

	/* this is for Solaris and SGI to force crash */
#	elif defined(__sun) || defined(__sgi)

	fp_except mask;

	/* >>chng 05 dec 30, accept FLUSH_DENORM_TO_ZERO as a synonym for HAVE_SUNMATH, PvH */
#	if defined(HAVE_SUNMATH) || defined(FLUSH_DENORM_TO_ZERO)

	/* >>chng 01 oct 09, disable gradual underflow on ultrasparc whith g++
	 * (only needed for versions < 3.1 or >= 4.3.0, see Note 1).
	 *
	 * compile this module with:
	 *     g++ [ other options... ] -I<include-dir> -DHAVE_SUNMATH -c cpu.cpp
	 * link the program with:
	 *     g++ -L<library-dir> -o cloudy.exe *.o -lsunmath
	 *
	 * you probably need to use -I<include-dir> and -L<library-dir> to point the
	 * compiler/linker to the location of the sunmath.h header file and libsunmath.so
	 * library (e.g., -I/opt/SUNWspro/prod/include/cc -L/opt/SUNWspro/lib; note that
	 * the actual location may vary from one installation to another).
	 * See also bug report 4487 on http://gcc.gnu.org/bugzilla/
	 *
	 * Note 1: Starting with g++ 3.1, bug 4487 has been solved: -funsafe-math-optimizations
	 * will automatically disable gradual underflow. Hence using nonstandard_arithmetic()
	 * is no longer necessary. The option -funsafe-math-optimizations should be included
	 * both when compiling and linking:
	 *
	 * g++ [ other options... ] -funsafe-math-optimizations -c *.c
	 * g++ [ other options... ] -funsafe-math-optimizations -o cloudy.exe *.o
	 *
	 * Starting with g++ 4.3.0 the -funsafe-math-optimizations option can no longer be
	 * used as it implicitly enables -fno-trapping-math, which is unsafe for Cloudy
	 * because we do trap floating point exceptions.
	 *
	 * Note 2: Don't use nonstandard_arithmetic() with CC (the SunWorks/Forte compiler);
	 * use the -fast commandline option instead to disable gradual underflow (or use
	 * -fnonstd if you don't want all the other options enabled by -fast). The option
	 * -fast (or -fnonstd) should be included both when compiling and linking:
	 *
	 * CC [ other options... ] -fast -c *.c
	 * CC -fast -o cloudy.exe *.o
	 *
	 * PvH */
	nonstandard_arithmetic();
#	endif

	/* enable floating point exceptions on sun and sgi */
	mask = fpgetmask();
	mask = mask | FP_X_INV | FP_X_OFL | FP_X_DZ;
	fpsetmask(mask);

#	elif defined(__alpha) && defined(__linux__) && defined(__GNUC__)

	/* the following is not supported on all hardware platforms, but certainly for EV56
	 * and later. earlier versions may work as well, but that has not been tested.
	 * for details see https://bugzilla.redhat.com/bugzilla/show_bug.cgi?id=51072 */
#	ifdef FE_NONIEEE_ENV
	/* this prevents the infamous math library bug when compiling with gcc on alpha-linux
	 * machines. if this doesn't work on your system, the only alternative is to link
	 * against the Compaq math library: gcc *.o -lcpml -lm, or use ccc itself, PvH */
	fesetenv(FE_NONIEEE_ENV);
#	endif

#	elif defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x5140

	// from Oracle Developer Studio 12.5 onwards, FP traps are not enebled by default
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#	endif
}

void t_cpu_i::set_signal_handlers()
{
	DEBUG_ENTRY( "set_signal_handlers()" );

#ifdef CATCH_SIGNAL
#	ifdef __unix
	p_action.sa_handler = &signal_handler;
	sigemptyset( &p_action.sa_mask );
	p_action.sa_flags = SA_NODEFER;

	p_default.sa_handler = SIG_DFL;
	sigemptyset( &p_default.sa_mask );
	p_default.sa_flags = SA_NODEFER;

	for( int sig=1; sig <= 31; sig++ )
	{
		// is the signal valid?
		if( sigaction( sig, NULL, NULL ) == 0 )
			// these two are for suspending and resuming a job
			if( sig != SIGSTOP && sig != SIGCONT )
				sigaction( sig, action(), NULL );
	}
#	endif

#	ifdef _MSC_VER
	signal( SIGABRT, &signal_handler );
	signal( SIGFPE, &signal_handler );
	signal( SIGILL, &signal_handler );
	signal( SIGINT, &signal_handler );
	signal( SIGSEGV, &signal_handler );
	signal( SIGTERM, &signal_handler );
#	endif
#endif
}

void t_cpu_i::signal_handler(int sig)
{
	// when an FPE is caught, the mask is reset...
	cpu.i().enable_traps();
#	ifdef _MSC_VER
	// at this point the signal handler has reverted to the default handler
	signal( sig, &signal_handler );
#	endif
	throw bad_signal( sig );
}


void t_cpu_i::printDataPath() const
{
	fprintf(ioQQQ, "The path is:\n");
	for( vector<string>::size_type i=1; i < chSearchPath.size(); ++i )
		fprintf( ioQQQ, "   ==%s==\n", chSearchPath[i].c_str() );
}

// this routine generates a list of all full paths to the locations where we should look for the file
void t_cpu_i::getPathList( const char* fname, vector<string>& PathList, access_scheme scheme, bool lgRead ) const
{
	DEBUG_ENTRY( "getPathList()" );

	vector<string>::size_type begin, end;

	switch( scheme )
	{
	case AS_DATA_ONLY:
	case AS_DATA_ONLY_TRY:
	case AS_DATA_OPTIONAL:
		begin = 1;
		end = chSearchPath.size();
		break;
	case AS_LOCAL_DATA:
	case AS_LOCAL_DATA_TRY:
		begin = 0;
		end = chSearchPath.size();
		break;
	case AS_LOCAL_ONLY:
	case AS_LOCAL_ONLY_TRY:
	case AS_SILENT_TRY:
		begin = 0;
		end = 1;
		break;
	case AS_DEFAULT:
		ASSERT( !lgRead );
		// these values are not used, but need to be there
		// to avoid warnings about using uninitialized vars
		begin = 0;
		end = 1;
		break;
	default:
		TotalInsanity();
	}

	PathList.clear();
	string FileName( fname );
	if( lgRead )
	{
		for( vector<string>::size_type i=begin; i < end; ++i )
			PathList.push_back( chSearchPath[i] + FileName );
	}
	else
	{
		PathList.push_back( FileName );
	}
}

void t_cpu_i::getMD5sums( const char* fname )
{
	DEBUG_ENTRY( "getMD5sums()" );

	// this routine reads a file with expected md5sum values for all Cloudy data files
	// they will be stored in the map md5sum_expct[] and can be used to compare to the
	// actual md5sum of the data files that were read in.

	md5sum_expct.clear();
	// we cannot use open_data() here, likely because the t_cpu_i c'tor hasn't completed yet
	// getPathList() would get a NULL this-pointer in that case and crash on a segfault...
	// so we have to copy a skeleton version of the open_data() code here...
	vector<string> PathList;
	getPathList( fname, PathList, AS_DATA_ONLY_TRY, true );
	fstream io;
	vector<string>::const_iterator ptr;
	for( ptr=PathList.begin(); ptr != PathList.end() && !io.is_open(); ++ptr )
		io.open( ptr->c_str(), mode_r );
	// if the file is not found, we return quietly
	if( !io.is_open() )
		return;
	++nFileDone;
	// the file has been opened, now parse the contents
	string line;
	char dirSep = '#';
	while( getline( io, line ) )
	{
		// determine what directory separator is used in the file
		if( dirSep == '#' )
		{
			if( line.find( '/' ) != string::npos )
				dirSep = '/';
			if( line.find( '\\' ) != string::npos )
				dirSep = '\\';
		}
		// replace the directory separator if the local OS uses a different one
		if( dirSep != '#' && dirSep != p_chDirSeparator )
			while( FindAndReplace( line, string(1,dirSep), string(1,p_chDirSeparator) ) ) {}
		string md5sum, path;
		istringstream iss( line );
		iss >> md5sum >> path;
		md5sum_expct[path] = md5sum;
	}
}

STATIC NORETURN void AbortErrorMessage( const char* fname, const vector<string>& PathList, access_scheme scheme )
{
	DEBUG_ENTRY( "AbortErrorMessage()" );

	if( scheme == AS_DATA_OPTIONAL )
		// presence is optional -> make warning less scary...
		fprintf( ioQQQ, "\nI could not open the data file %s\n\n", fname );
	else
		fprintf( ioQQQ, "\nPROBLEM DISASTER I could not open the data file %s\n\n", fname );
	if( cpu.i().firstOpen() || scheme == AS_DATA_ONLY )
	{
		// failed on very first open -> most likely path is not correct
		// failed on AS_DATA_ONLY -> CLOUDY_DATA_PATH may point to obsolete data dir
		fprintf( ioQQQ, "Although there may be other reasons you have received this error,\n");
		fprintf( ioQQQ, "the most likely are that the path has not been properly set or\n");
		fprintf( ioQQQ, "that the path points to an old version of the data. It should\n");
		fprintf( ioQQQ, "point to the data directory you downloaded from the web site.\n\n");
		fprintf( ioQQQ, "Please use \"make\" to compile the code. This will automatically\n");
		fprintf( ioQQQ, "set the path correctly. Alternatively you can set the environment\n");
		fprintf( ioQQQ, "variable CLOUDY_DATA_PATH to point to the data directory using\n");
		fprintf( ioQQQ, "the shell command \nexport CLOUDY_DATA_PATH=\"/path/to/data\"\n");
		fprintf( ioQQQ, "from a bash command prompt.\n\n");
		cpu.i().printDataPath();
	}
	else
	{
		// failed on search including local directory -> most likely the file name
		// was mistyped on a compile command, or Cloudy is run in the wrong directory
		// if scheme == AS_DATA_OPTIONAL, this most likely is a stellar grid that is not installed.
		fprintf( ioQQQ, "These are all the paths I tried:\n" );
		for( vector<string>::const_iterator ptr=PathList.begin(); ptr != PathList.end(); ++ptr )
			fprintf( ioQQQ, "   ==%s==\n", ptr->c_str() );
		// AS_DATA_OPTIONAL files should provide their own message (currently only stellar grids)
		if( scheme != AS_DATA_OPTIONAL )
		{
			fprintf( ioQQQ, "\nAlthough there may be other reasons you have received this error,\n");
			fprintf( ioQQQ, "the most likely are that you mistyped the file name, or that you\n");
			fprintf( ioQQQ, "are running Cloudy in the wrong directory. If you are running a\n");
			fprintf( ioQQQ, "COMPILE command, this needs to be done in the data directory.\n\n");
			fprintf( ioQQQ, "It is is also possible that the path has not been properly set. It\n");
			fprintf( ioQQQ, "should point to the data directory you downloaded from the web site.\n");
			fprintf( ioQQQ, "Please use \"make\" to compile the code. This will automatically\n");
			fprintf( ioQQQ, "set the path correctly. Alternatively you can set the environment\n");
			fprintf( ioQQQ, "variable CLOUDY_DATA_PATH to point to the data directory using\n");
			fprintf( ioQQQ, "the shell command \nexport CLOUDY_DATA_PATH=\"/path/to/data\"\n");
			fprintf( ioQQQ, "from a bash command prompt.\n\n");
		}
	}
	fprintf(ioQQQ, "Sorry.\n\n\n");
	cdEXIT(EXIT_FAILURE);
}

STATIC string check_mult_path( const char* fname, const vector<string>& PathList, access_scheme scheme, bool lgRead )
{
	DEBUG_ENTRY( "check_mult_path()" );

	if( !lgRead )
	{
		ASSERT( PathList.size() == 1 );
		if( trace.lgTrace && scheme != AS_SILENT_TRY )
			fprintf( ioQQQ, " open_data writing %s\n", PathList[0].c_str() );
		return PathList[0];
	}

	vector<string>::const_iterator ptr;
	vector<string> PathSuccess;
	for( ptr=PathList.begin(); ptr != PathList.end(); ++ptr )
	{
		FILE* handle = sys_fopen( ptr->c_str(), "r" );
		if( trace.lgTrace && scheme != AS_SILENT_TRY )
		{
			fprintf( ioQQQ, " open_data trying to read %s found %c", ptr->c_str(), TorF(handle != NULL) );
			if( handle != NULL )
				fprintf( ioQQQ, " used %c", TorF(PathSuccess.size() == 0) );
			fprintf( ioQQQ, "\n" );
		}
		if( handle != NULL )
		{
			// don't store path if it is identical to a previous one
			if( find( PathSuccess.begin(), PathSuccess.end(), *ptr ) == PathSuccess.end() )
				PathSuccess.push_back( *ptr );
			fclose( handle );
		}
	}

	if( PathSuccess.size() > 1 && scheme != AS_SILENT_TRY )
	{
		fprintf( ioQQQ, "CAUTION: multiple matches for file %s found:\n", fname );
		for( size_t i=0; i < PathSuccess.size(); ++i )
			fprintf( ioQQQ, "   ==%s==\n", PathSuccess[i].c_str() );
		fprintf( ioQQQ, "Using the first match.\n" );
	}

	// return the first successful match, or an empty string if no match was found
	return ( PathSuccess.size() > 0 ) ? PathSuccess[0] : "";
}

FILE* open_data( const char* fname, const char* mode, access_scheme scheme )
{
	DEBUG_ENTRY( "open_data()" );

	// for mode "r" and "rb" the default is AS_DATA_ONLY, and for all
	// other modes AS_LOCAL_ONLY since the latter can overwrite the file
	string m = mode;
	bool lgRead = ( m == "r" || m == "rb" );
	if( lgRead && scheme == AS_DEFAULT )
		scheme = AS_DATA_ONLY;

	bool lgAbort = ( scheme == AS_DATA_ONLY || scheme == AS_DATA_OPTIONAL ||
			 scheme == AS_LOCAL_DATA || scheme == AS_LOCAL_ONLY ||
			 scheme == AS_DEFAULT );

	vector<string> PathList;
	cpu.i().getPathList( fname, PathList, scheme, lgRead );

	FILE* handle = NULL;
	string path = check_mult_path( fname, PathList, scheme, lgRead );
	if( path != "" )
	{
		if( lgRead )
			check_data( path.c_str(), fname );
		handle = sys_fopen( path.c_str(), mode );
	}

	if( handle == NULL )
	{
		if( lgAbort )
			AbortErrorMessage( fname, PathList, scheme );
	}
	else
		++cpu.i().nFileDone;

	return handle;
}

void open_data( fstream& stream, const char* fname, ios_base::openmode mode, access_scheme scheme )
{
	DEBUG_ENTRY( "open_data()" );

	// for mode_r and mode_rb the default is AS_DATA_ONLY, and for all
	// other modes AS_LOCAL_ONLY since the latter can overwrite the file
	bool lgRead = ( (mode&ios_base::out) == 0 );
	if( lgRead && scheme == AS_DEFAULT )
		scheme = AS_DATA_ONLY;

	bool lgAbort = ( scheme == AS_DATA_ONLY || scheme == AS_DATA_OPTIONAL ||
			 scheme == AS_LOCAL_DATA || scheme == AS_LOCAL_ONLY ||
			 scheme == AS_DEFAULT );

	vector<string> PathList;
	cpu.i().getPathList( fname, PathList, scheme, lgRead );

	ASSERT( !stream.is_open() );
	string path = check_mult_path( fname, PathList, scheme, lgRead );
	if( path != "" )
	{
		if( lgRead )
			check_data( path.c_str(), fname );
		stream.open( path.c_str(), mode );
	}

	if( !stream.is_open() )
	{
		if( lgAbort )
			AbortErrorMessage( fname, PathList, scheme );
	}
	else
		++cpu.i().nFileDone;
}

MPI_File open_data( const char* fname, int mode, access_scheme scheme )
{
	DEBUG_ENTRY( "open_data()" );

	// for mpi_mode_r the default is AS_DATA_ONLY, and for all
	// other modes AS_LOCAL_ONLY since the latter can overwrite the file
	bool lgRead = ( mode == mpi_mode_r );
	if( lgRead && scheme == AS_DEFAULT )
		scheme = AS_DATA_ONLY;

	bool lgAbort = ( scheme == AS_DATA_ONLY || scheme == AS_DATA_OPTIONAL ||
			 scheme == AS_LOCAL_DATA || scheme == AS_LOCAL_ONLY ||
			 scheme == AS_DEFAULT );

	vector<string> PathList;
	cpu.i().getPathList( fname, PathList, scheme, lgRead );

	int err = MPI_ERR_INTERN;
	MPI_File fh = MPI_FILE_NULL;
	string path = check_mult_path( fname, PathList, scheme, lgRead );
	if( path != "" )
	{
		if( lgRead )
			check_data( path.c_str(), fname );
		err = MPI_File_open( MPI_COMM_WORLD, const_cast<char*>(path.c_str()), mode, MPI_INFO_NULL, &fh );
	}

	if( err != MPI_SUCCESS )
	{
		if( lgAbort )
			AbortErrorMessage( fname, PathList, scheme );
		// just to be safe, the man page is not clear on this...
		fh = MPI_FILE_NULL;
	}
	else
		++cpu.i().nFileDone;

	return fh;
}

void check_data( const char* fpath, const char* fname )
{
	DEBUG_ENTRY( "check_data()" );

	if( !( t_version::Inst().lgRelease || t_version::Inst().lgReleaseBranch ) || !prt.lgPrintTime )
		return;

	map<string,string>::const_iterator ptr = cpu.i().md5sum_expct.find( fname );
	if( ptr != cpu.i().md5sum_expct.end() )
	{
		fstream ioFile( fpath, mode_r );
		if( ioFile.is_open() )
		{
			string md5sum = MD5datastream( ioFile );
			if( md5sum != ptr->second )
			{
				fprintf( ioQQQ, "NOTE: using modified data in %s.\n", fname );
				++cpu.i().nMD5Mismatch;
			}
		}
	}
}

/** define routines for setting single and double precision signaling NaN
 * The bit pattern for an SNaN is implementation defined, but this should
 * work on most modern CPU's. The system definition is preferred, so in
 * C++ this should be replaced by numeric_limits<TYPE>::signaling_NaN() */

void set_NaN(sys_float &x)
{
	if( sizeof(sys_float) == 4 )
		*reinterpret_cast<int32*>(&x) = cpu.i().Float_SNaN_Value;
	else
		x = -FLT_MAX;
}

void set_NaN(sys_float x[], /* x[n] */
	     long n)
{
	long i;

	if( sizeof(sys_float) == 4 )
	{
		int32 *y = reinterpret_cast<int32*>(x);
		for( i=0; i < n; i++ )
			*y++ = cpu.i().Float_SNaN_Value;
	}
	else
	{
		for( i=0; i < n; i++ )
			x[i] = -FLT_MAX;
	}
}

void set_NaN(double &x)
{
	if( sizeof(double) == 8 )
	{
#		ifdef HAVE_INT64
		*reinterpret_cast<int64*>(&x) = cpu.i().Double_SNaN_Value;
#		else
		int32 *y = reinterpret_cast<int32*>(&x);
		*y++ = cpu.i().Double_SNaN_Value[0];
		*y = cpu.i().Double_SNaN_Value[1];
#		endif
	}
	else
		x = -DBL_MAX;
}

/* set_NaN - set NaN */
void set_NaN(double x[], /* x[n] */
	     long n)
{
	long i;

	if( sizeof(double) == 8 )
	{
#		ifdef HAVE_INT64
		int64 *y = reinterpret_cast<int64*>(x);
		for( i=0; i < n; i++ )
			*y++ = cpu.i().Double_SNaN_Value;
#		else
		int32 *y = reinterpret_cast<int32*>(x);
		for( i=0; i < n; i++ )
		{
			*y++ = cpu.i().Double_SNaN_Value[0];
			*y++ = cpu.i().Double_SNaN_Value[1];
		}
#		endif
	}
	else
	{
		for( i=0; i < n; i++ )
			x[i] = -DBL_MAX;
	}
}

/** detect quiet and signaling NaNs in single precision FP */
bool MyIsnan(const sys_float &x)
{
	if( sizeof(sys_float) == 4 && FLT_MAX_EXP-FLT_MIN_EXP+3 == 256 )
	{
		const int32 *p = reinterpret_cast<const int32*>(&x);
		int32 r = *p & 0x7f800000; r ^= 0x7f800000;
		int32 s = *p & 0x007fffff;
		return ( r == 0 && s != 0 );
	}
	else
		/* we don't understand this CPU */
		return false;
}

/** detect quiet and signaling NaNs in double precision FP */
bool MyIsnan(const double &x)
{
	if( sizeof(double) == 8 && DBL_MAX_EXP-DBL_MIN_EXP+3 == 2048 )
	{
#		ifdef HAVE_INT64
		const int64 *p = reinterpret_cast<const int64*>(&x);
		int64 r = *p & INT64_LIT(0x7ff0000000000000); r ^= INT64_LIT(0x7ff0000000000000);
		int64 s = *p & INT64_LIT(0x000fffffffffffff);
		return ( r == 0 && s != 0 );
#		else
		const int32 *p = reinterpret_cast<const int32*>(&x);
		if( cpu.i().little_endian() )
		{
			int32 r = p[1] & 0x7ff00000; r ^= 0x7ff00000;
			int32 s = p[1] & 0x000fffff; s |= p[0];
			return ( r == 0 && s != 0 );
		}
		else if( cpu.i().big_endian() )
		{
			int32 r = p[0] & 0x7ff00000; r ^= 0x7ff00000;
			int32 s = p[0] & 0x000fffff; s |= p[1];
			return ( r == 0 && s != 0 );
		}
		else
			/* we don't understand this CPU */
			return false;
#		endif
	}
	else
		/* we don't understand this CPU */
		return false;
}
