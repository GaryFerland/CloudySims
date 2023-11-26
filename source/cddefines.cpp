/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* out-of-line constructor for assert -- put breakpoint in this
   routine to trap assert throws for IDEs without built-in facility. */
#include "cddefines.h"
#include "prt.h"
#include "lines_service.h"

FILE *ioQQQ;
FILE *ioStdin;
FILE* ioPrnErr;
bool lgTestCodeCalled; 
bool lgTestCodeEnabled;
bool lgPrnErr;
long int nzone;
double fnzone;
long int iteration;

bad_signal::bad_signal(int sig, void* ptr) : p_sig(sig)
{
	cpu.i().GenerateBacktrace(ptr);
}

bad_assert::bad_assert(const char* file, long line, const char* comment) :
	p_file(file), p_line(line), p_comment(comment)
{
	cpu.i().GenerateBacktrace(NULL);
}

cloudy_abort::cloudy_abort(const char* comment) : p_comment(comment)
{
	cpu.i().GenerateBacktrace(NULL);
}

realnum t_wavl::p_convertWvl() const
{
	if( p_type == WL_AIR || ( p_type == WL_NATIVE && prt.lgPrintLineAirWavelengths ) )
		// wavelength may be negative
		return ( p_wavl < 0_r ) ? -wlAirVac(-p_wavl) : wlAirVac(p_wavl);
	else
		return p_wavl;
}

string t_wavl::p_convertWvlString() const
{
	string chBuf;
	realnum wvac = p_convertWvl();
	sprt_wl(chBuf, wvac);
	return chBuf;
}
