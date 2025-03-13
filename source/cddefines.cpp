/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* out-of-line constructor for assert -- put breakpoint in this
   routine to trap assert throws for IDEs without built-in facility. */
#include "cddefines.h"
#include "lines.h"
#include "prt.h"

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
		return p_wlAirVac();
	else
		return p_wavl;
}

/* compute wavelength in vacuum given air wavelength */
realnum t_wavl::p_wlAirVac() const
{
	DEBUG_ENTRY( "p_wlAirVac()" );

	// wavelength may be negative
	double wlAir = fabs(p_wavl);
	double wlVacuum = wlAir;

	if( wlAir > 2000. )
	{
		double RefIndex_v = 1.;
		// iterate since wavenumber depends on wlVac not wlAir, but difference should be small
		for( int i=0; i < 2; ++i )
		{
			RefIndex_v = p_RefIndex(1e8 / wlVacuum);
			wlVacuum = wlAir * RefIndex_v;
		}
	}

	return sign(realnum(wlVacuum), p_wavl);
}

/* calculate the index of refraction for STP air using the line energy in wavenumbers */
double t_wavl::p_RefIndex(double EnergyWN) const
{
	DEBUG_ENTRY( "p_RefIndex()" );

	ASSERT( EnergyWN > 0. );

	double RefIndex_v = 1.0;

	/* only do index of refraction if longward of 2000A */
	if( EnergyWN < 5e4 )
	{
		/* xl is wavenumber in micron^-1, squared */
		double xl = pow2(EnergyWN * 1e-4);
		/* use a formula from 
		 *>>refer	air	index refraction	Peck & Reeder 1972, JOSA, 62, 8, 958 */
		RefIndex_v += 1e-8 * (8060.51 + 2480990.0 / (132.274 - xl) + 17455.7 / (39.32957 - xl));
	}

	ASSERT( RefIndex_v >= 1. );
	return RefIndex_v;
}

/* write wavelength to string */
string t_wavl::sprt_wl(const char* format) const
{
	DEBUG_ENTRY( "sprt_wl()" );

	realnum wl = wavlVac();

	if( prt.lgPrintLineAirWavelengths && wl > 0_r )
		wl /= p_RefIndex(1.e8/wl);

	if( format != NULL )
	{
		auto buflen = snprintf(NULL, 0, format, wl); // dry run to get buffer length
		vector<char> buf(buflen+1); // buflen does not include the terminating 0-byte
		snprintf(buf.data(), buflen+1, format, wl);
		return string(buf.data());
	}

	/* print in angstrom unless > 1e4, then use micron */
	string chUnits;
	if( wl > 1e8 )
	{
		/* centimeter */
		chUnits = "c";
		wl /= 1e8;
	}
	else if( wl > 1e4 )
	{
		/* micron */
		chUnits = "m";
		wl /= 1e4;
	}
	else if( wl == 0. )
	{
		chUnits = " ";
	}
	else
	{
		/* angstrom units */
		chUnits = "A";
	}

	ostringstream oss;
	/* want total of LineSave.sig_figs sig figs */
	if( wl==0. )
	{
		oss << setw(LineSave.wl_length-1) << 0;
	}
	else
	{
		int n = LineSave.sig_figs - 1 - (int)log10(wl);
		if( n > 0 )
		{
			oss << fixed << setw(LineSave.sig_figs+1) << setprecision(n) << wl;
		}
		else if( wl < (realnum)INT_MAX )
		{
			oss << setw(LineSave.sig_figs+1) << (int)wl;
		}
		else
		{
			oss << setw(LineSave.sig_figs+1) << "*";
		}
	}
	return oss.str() + chUnits;
}

/* write wavelength to output stream */
void t_wavl::prt_wl(FILE *ioOUT, const char* format) const
{
	DEBUG_ENTRY( "prt_wl()" );

	if( format != NULL )
		fprintf(ioOUT, "%s", sprt_wl(format).c_str() );
	else if( LineSave.wl_length > 0 )
		fprintf(ioOUT, "%.*s", LineSave.wl_length, sprt_wl().c_str() );
	else
		fprintf(ioOUT, "%s", sprt_wl().c_str() );
	return;
}
