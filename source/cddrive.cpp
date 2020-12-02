/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdDrive main routine to call cloudy under all circumstances) */
/*cdReasonGeo write why the model stopped and type of geometry on io file */
/*cdWarnings write all warnings entered into comment stack */
/*cdEms obtain the local emissivity for a line, for the last computed zone */
/*cdColm get the column density for a constituent  */
/*cdLine get the predicted line intensity, also index for line in stack */
/*cdLine_ip get the predicted line intensity, using index for line in stack */
/*cdCautions print out all cautions after calculation, on arbitrary io unit */
/*cdTemp_last routine to query results and return temperature of last zone */
/*cdDepth_depth get depth structure from previous iteration */
/*cdTimescales returns thermal, recombination, and H2 formation timescales */
/*cdSurprises print out all surprises on arbitrary unit number */
/*cdNotes print stack of notes about current calculation */
/*cdPressure_last routine to query results and return pressure of last zone */
/*cdTalk tells the code whether to print results or be silent */
/*cdOutput redirect output to arbitrary file */
/*cdInput redirect input from arbitrary file */
/*cdRead routine to read in command lines when cloudy used as subroutine */
/*cdErrors produce summary of all warnings, cautions, etc, on arbitrary io unit */
/*cdIonFrac get ionization fractions for a constituent */
/*cdTemp get mean electron temperature for any element */
/*cdCooling_last routine to query results and return cooling of last zone */
/*cdHeating_last routine to query results and return heating of last zone */
/*cdEDEN_last return electron density of last zone */
/*cdNoExec call this routine to tell code not to actually execute */
/*cdDate - puts date of code into string */
/*cdVersion produces string that gives version number of the code */
/*cdExecTime any routine can call this, find the time [s] since cdInit was called */
/*cdPrintCommands( FILE *) prints all input commands into file */
/*cdDrive main routine to call cloudy under all circumstances) */
/*cdNwcns get the number of cautions and warnings, to tell if calculation is ok */
/*debugLine provides a debugging hook into the main line array  */
/*cdEms_ip obtain the local emissivity for a line with known index */
/*cdnZone gets number of zones */
/*cdClosePunchFiles closes all the save files that have been used */
/*cdB21cm - returns B as measured by 21 cm */
/*cdPrtWL print line wavelengths in Angstroms in the standard format */
#include "cddefines.h"
#include "cddrive.h"

#include "trace.h"
#include "conv.h"
#include "init.h"
#include "lines.h"
#include "pressure.h"
#include "dense.h"
#include "radius.h"
#include "struc.h"
#include "mole.h"
#include "elementnames.h"
#include "mean.h"
#include "phycon.h"
#include "called.h"
#include "parse.h"
#include "input.h"
#include "taulines.h"
#include "version.h"
#include "thermal.h"
#include "grid.h"
#include "timesc.h"
#include "cloudy.h"
#include "warnings.h"
#include "broke.h"
#include "iso.h"
#include "save.h"
#include "rfield.h"
#include "lines_service.h"
#include "service.h"
#include "parser.h"
#include "generic_state.h"
#include "prt.h"

/*************************************************************************
 *
 * cdDrive - main routine to call cloudy - returns 0 if all ok, 1 if problems
 *
 ************************************************************************/

int cdDrive()
{
	bool lgBAD;

	DEBUG_ENTRY( "cdDrive()" );
	/*********************************
	 * main routine to call cloudy   *
	 *********************************/

	/* this is set false when code loaded, set true when cdInit called,
	 * this is check that cdInit was called first */
	if( !lgcdInitCalled )
	{
		printf(" cdInit was not called first - this must be the first call.\n");
		cdEXIT(EXIT_FAILURE);
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			"cdDrive: lgOptimr=%1i lgVaryOn=%1i lgNoVary=%1i input.nSave:%li\n",
			optimize.lgOptimr , optimize.lgVaryOn , optimize.lgNoVary, input.nSave );
	}

	/* should we call cloudy, or the optimization driver? 
	 * possible to have VARY on line without OPTIMIZE being set 
	 * lgNoVary set with "no optimize" command */
	if( optimize.lgOptimr && optimize.lgVaryOn && !optimize.lgNoVary && !grid.lgInsideGrid )
		optimize.lgVaryOn = true;
	else
		optimize.lgVaryOn = false;

	/* set up the frequency mesh, the SET CONTINUUM RESOLUTION factor was already parsed in cdRead() */
	/* the bool 'false' inidicates that the cell for the unit test is the highest entry in the regular mesh */
	rfield.InitMesh(false);
	rfield.nflux_with_check = rfield.ncells();
	/* there must be a cell above nflux for us to pass unity through the vol integrator */
	/*>>chng 04 oct 10, from nflux = nupper to nupper-1 since vectors must malloc to nupper, but
	 * will address [nflux] for unit continuum test */
	rfield.nflux = rfield.nflux_with_check-1;
	rfield.nPositive = rfield.nflux;

	/* one time initialization of core load - returns if already called 
	 * called here rather than in cdInit since at this point we know if
	 * single sim or grid */
	InitCoreload();

	if( optimize.lgVaryOn )
	{
		/* this branch called if optimizing or grid calculation */
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: calling grid_do\n");
		/* option to drive optimizer set if OPTIMIZE was in input stream */
		lgBAD = grid_do();
	}
	else
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: calling cloudy\n");

		/* optimize did not occur, only compute one model, call cloudy */
		lgBAD = cloudy();
	}

	/* reset flag saying that cdInit has not been called */
	lgcdInitCalled = false;

	if( lgAbort || lgBAD )
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: returning failure during call. \n");
		/* lgAbort set true if something wrong, so return lgBAD false. */
		return 1;
	}
	else
	{
		/* everything is ok, return 0 */
		return 0;
	}
}


/*************************************************************************
 *
 * cdReasonGeo write why the model stopped and type of geometry on io file 
 *
 ************************************************************************/

/*cdReasonGeo write why the model stopped and type of geometry on io file */
void cdReasonGeo(FILE * ioOUT)
{

	DEBUG_ENTRY( "cdReasonGeo()" );

	/*this is the reason the calculation stopped*/
	fprintf( ioOUT, "%s", warnings.chRgcln[0] );
	fprintf( ioOUT , "\n" );
	/* this is the geometry */
	fprintf( ioOUT, "%s", warnings.chRgcln[1] );
	fprintf( ioOUT , "\n" );
	return;
}


/*************************************************************************
 *
 * cdWarnings write all warnings entered into comment stack 
 *
 ************************************************************************/

/*cdWarnings write all warnings entered into comment stack */
void cdWarnings(FILE *ioPNT )
{
	long int i;

	DEBUG_ENTRY( "cdWarnings()" );

	if( warnings.nwarn > 0 )
		{
			for( i=0; i < warnings.nwarn; i++ )
			{
				/* these are all warnings that were entered in comment */
				fprintf( ioPNT, "%s", warnings.chWarnln[i] );
				fprintf( ioPNT, "\n" );
			}
		}

	return;
}


/*************************************************************************
 *
 * cdCautions print out all cautions after calculation, on arbitrary io unit 
 *
 ************************************************************************/

/*cdCautions print out all cautions after calculation, on arbitrary io unit */
void cdCautions(FILE * ioOUT)
{
	long int i;

	DEBUG_ENTRY( "cdCautions()" );

	if( warnings.ncaun > 0 )
	{
		for( i=0; i < warnings.ncaun; i++ )
		{
			fprintf( ioOUT, "%s", warnings.chCaunln[i] );
			fprintf( ioOUT, "\n" );
		}
	}
	return;
}

/*************************************************************************
 *
 * cdTimescales returns thermal, recombination, and H2 formation timescales 
 *
 ************************************************************************/

void cdTimescales(
	/* the thermal cooling timescale */
	double *TTherm , 
	/* the hydrogen recombination timescale */
	double *THRecom , 
	/* the H2 formation timescale */
	double *TH2 )
{

	DEBUG_ENTRY( "cdTimescales()" );

	/* these were all evaluated in AgeCheck, which was called by PrtComment */

	/* thermal or cooling timescale */
	*TTherm = timesc.time_therm_long;

	/* the hydrogen recombination timescale */
	*THRecom = timesc.time_Hrecom_long;

	/* longer of the the H2 formation and destruction timescales */
	*TH2 = MAX2( timesc.time_H2_Dest_longest , timesc.time_H2_Form_longest );
	return;
}


/*************************************************************************
 *
 * cdSurprises print out all surprises on arbitrary unit number 
 *
 ************************************************************************/

/*cdSurprises print out all surprises on arbitrary unit number */
void cdSurprises(FILE * ioOUT)
{
	long int i;

	DEBUG_ENTRY( "cdSurprises()" );

	if( warnings.nbang > 0 )
	{
		for( i=0; i < warnings.nbang; i++ )
		{
			fprintf( ioOUT, "%s\n", warnings.chBangln[i] );
		}
	}

	if  ( broke.lgPrintFixits )
	{
		fprintf (ioOUT, "\nActive fixits for this run\n");
		for (vector<string>::iterator it = FixitList::Inst().list.begin();
			  it != FixitList::Inst().list.end(); ++it)
		{
			fprintf(ioOUT,"%s\n",it->c_str());
		}
		fprintf (ioOUT, "\n");
	}

	return;
}


/*************************************************************************
 *
 * cdNotes print stack of notes about current calculation 
 *
 ************************************************************************/

/*cdNotes print stack of notes about current calculation */
void cdNotes(FILE * ioOUT)
{
	long int i;

	DEBUG_ENTRY( "cdNotes()" );

	if( warnings.nnote > 0 )
	{
		for( i=0; i < warnings.nnote; i++ )
		{
			fprintf( ioOUT, "%s", warnings.chNoteln[i] );
			fprintf( ioOUT, "\n" );
		}
	}
	return;
}

/*************************************************************************
 *
 * cdCooling_last routine to query results and return cooling of last zone 
 *
 ************************************************************************/

/*cdCooling_last routine to query results and return cooling of last zone */
double cdCooling_last() /* return cooling for last zone */
{
	return thermal.ctot;
}

/*************************************************************************
 *
 * cdVersion - puts version number of code into string 
 * incoming string must have sufficient length and will become null
 * terminated string
 *
 ************************************************************************/

void cdVersion(char chString[])
{
	strcpy( chString , t_version::Inst().chVersion );
	return;
}

/*************************************************************************
 *
 * cdDate - puts date of code into string 
 * incoming string must have at least 8 char and will become null
 * terminated string
 *
 ************************************************************************/

/* cdDate - puts date of code into string  */
void cdDate(char chString[])
{
	strcpy( chString , t_version::Inst().chDate );
	return;
}


/*************************************************************************
 *
 * cdHeating_last routine to query results and return heating of last zone
 *
 ************************************************************************/

/*cdHeating_last routine to query results and return heating of last zone */
double cdHeating_last() /* return heating for last zone */
{
	return thermal.htot;
}


/*************************************************************************
 *
 * cdEDEN_last return electron density of last zone
 *
 ************************************************************************/

/*cdEDEN_last return electron density of last zone */
double cdEDEN_last() /* return electron density for last zone */
{
	return dense.eden;
}

/*************************************************************************
 *
 * cdNoExec call this routine to tell code not to actually execute
 *
 ************************************************************************/

/*cdNoExec call this routine to tell code not to actually execute */
#include "noexec.h"

void cdNoExec()
{

	DEBUG_ENTRY( "cdNoExec()" );

	/* option to read in all input quantiles and NOT execute the actual model
	 * only check on input parameters - set by calling cdNoExec */
	noexec.lgNoExec = true;

	return;
}


/*************************************************************************
 *
 * cdSetExecTime routine to initialize variables keeping track of time at start of calculation
 *
 ************************************************************************/

/* set to false initially, then to true when cdSetExecTime() is called to
 * initialize the clock */
static bool lgCalled=false;

#if defined(_MSC_VER) 
/* _MSC_VER branches assume getrusage isn't implemented by MS */
struct timeval {
	long    tv_sec;         /* seconds */
	long    tv_usec;        /* microseconds */
};
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

/* will be used to save initial time */
static struct timeval before;

/* cdClock stores time since arbitrary datum in clock_dat           */
STATIC void cdClock(struct timeval *clock_dat)
{
	DEBUG_ENTRY( "cdClock()" );

/* >>chng 06 sep 2 rjrw: use long recurrence, fine grain UNIX clock *
 * -- maintain system dependences in a single routine               */
#if defined(_MSC_VER)
	clock_t clock_val;
	/* >>chng 05 dec 21, from above to below due to negative exec times */
	/*return (double)(clock() - before) / (double)CLOCKS_PER_SEC;*/
	clock_val = clock();
	clock_dat->tv_sec = clock_val/CLOCKS_PER_SEC;
	clock_dat->tv_usec = 1000000*(clock_val-(clock_dat->tv_sec*CLOCKS_PER_SEC))/CLOCKS_PER_SEC;
	/*>>chng 06 oct 05, this produced very large number, time typically 50% too long 
	clock_dat->tv_usec = 0;*/
#else
	struct rusage rusage;
	if(getrusage(RUSAGE_SELF,&rusage) != 0)
	{ 
		fprintf( ioQQQ, "DISASTER cdClock called getrusage with invalid arguments.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	clock_dat->tv_sec = rusage.ru_utime.tv_sec;
	clock_dat->tv_usec = rusage.ru_utime.tv_usec;
#endif

}
/* cdSetExecTime is called by cdInit when everything is initialized,
 * so that every time cdExecTime is called the elapsed time is returned */
void cdSetExecTime()
{
	cdClock(&before);
	lgCalled = true;
}
/* cdExecTime returns the elapsed time cpu time (sec) that has elapsed 
 * since cdInit called cdSetExecTime.*/
double cdExecTime()
{
	DEBUG_ENTRY( "cdExecTime()" );

	struct timeval clock_dat;
	/* check that we were properly initialized */
	if( lgCalled)
	{
		cdClock(&clock_dat);
		/*fprintf(ioQQQ,"\n DEBUG sec %.2e usec %.2e\n",
			(double)(clock_dat.tv_sec-before.tv_sec),
			1e-6*(double)(clock_dat.tv_usec-before.tv_usec));*/
		return (double)(clock_dat.tv_sec-before.tv_sec)+1e-6*(double)(clock_dat.tv_usec-before.tv_usec);
	}
	else
	{
		/* this is a big problem, we were asked for the elapsed time but
		 * the timer was not initialized by calling SetExecTime */
		fprintf( ioQQQ, "DISASTER cdExecTime was called before cdSetExecTime, impossible.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
}

// Maximum memory used, in kB
long cdMemory()
{
#if defined(_MSC_VER)
	return 0L;
#else
	struct rusage usage;
	if(getrusage(RUSAGE_SELF,&usage) != 0)
	{ 
		fprintf( ioQQQ, "DISASTER cdMemory called getrusage with invalid arguments.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return usage.ru_maxrss;
#endif
}

/*************************************************************************
 *
 * cdPrintCommands prints all input commands into file
 *
 ************************************************************************/

/* cdPrintCommands( FILE *)
 * prints all input commands into file */
void cdPrintCommands( FILE * ioOUT )
{
	long int i;
	fprintf( ioOUT, " Input commands follow:\n" );
	fprintf( ioOUT, "c ======================\n" );

	for( i=0; i <= input.nSave; i++ )
	{
		// exclude lines from init files
		if( input.InclLevel[i] == 0 )
			fprintf( ioOUT, "%s\n", input.chCardSav[i] );
	}
	fprintf( ioOUT, "c ======================\n" );
}


/*************************************************************************
 *
 * cdEmis obtain the local emissivity for a line, for the last computed zone
 *
 ************************************************************************/

/** cdEms obtain the local emissivity for a line, for the last computed zone */


void cdEmis(
	const char *chLabel,
	realnum wavelength, 
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// intrinsic or emergent
	bool lgEmergent )
{
	DEBUG_ENTRY( "cdEmis()" );

	long ipobs = LineSave.findline(chLabel, wavelength);
	if (ipobs)
		cdEmis(&LineSave.lines[ipobs],emiss,lgEmergent);
	else
		*emiss = 0.0;
	return;
}

void cdEmis_ip(
	/* index of the line in the stack */
	long int ipLine, 
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// intrinsic or emergent
	bool lgEmergent )
{
	DEBUG_ENTRY( "cdEmis_ip()" );

	/* returns the emissivity in a line - implements save lines structure 
	 * this version uses previously stored line index to save speed */
	ASSERT( ipLine >= 0 && ipLine < LineSave.nsum );
	cdEmis(&LineSave.lines[ipLine],emiss,lgEmergent);
}

/*************************************************************************
 *
 * cdColm get the column density for a constituent - 0 return if ok, 1 if problems 
 *
 ************************************************************************/

int cdColm(
	/* return value is zero if all ok, 1 if errors happened */
	/* 4-char + eol string that is first 
	 * 4 char of element name as spelled by Cloudy, upper or lower case */
	const char *chLabel,	

	/* integer stage of ionization, 1 for atom, 2 for A+, etc, 
	 * 0 is special flag for CO, H2, OH, or excited state */
	long int ion,

	/* the column density derived by the code [cm-2] */
	double *theocl )	
{
	long int nelem;
	/* use following to store local image of character strings */
	char chLABEL_CAPS[NCHLAB];

	DEBUG_ENTRY( "cdColm()" );
	*theocl = 0.;

	if( strlen(chLabel) > NCHLAB-1 )
	{
		fprintf( ioQQQ, " PROBLEM cdColm called with insane chLabel (between quotes) \"%s\", must be no more than %d characters long.\n",
					chLabel, NCHLAB-1 );
		return 1;
	}

	strcpy( chLABEL_CAPS, chLabel );

	/* convert element label to all caps */
	caps(chLABEL_CAPS);
	trimTrailingWhiteSpace( chLABEL_CAPS );

	genericState gs;

	/* zero ionization stage has special meaning.  The quantities recognized are
	 * the molecules, "H2  ", "OH  ", "CO  ", etc 
	 * "CII*" excited state C+ */
	if( ion < 0 )
	{
		fprintf( ioQQQ, " PROBLEM cdColm called with insane ion, =%li\n",
		  ion );
		return 1;
	}

	else if( ion == 0 )
	{
		/* this case molecular column density */
		/* want the *total* molecular hydrogen H2 column density */
		if( strcmp( chLABEL_CAPS , "H2" )==0 )
		{
			*theocl = findspecieslocal("H2")->column + findspecieslocal("H2*")->column;
		}
		/* H2g - ground H2 column density */
		else if( strcmp( chLABEL_CAPS , "H2G" )==0 )
		{
			gs.sp = findspecieslocal("H2");
			*theocl = column(gs);
		}
 		/* H2* - excited H2 - column density */
		else if( strcmp( chLABEL_CAPS , "H2*" )==0 )
		{
			gs.sp = findspecieslocal("H2*");
			*theocl = column(gs);
		}
		/* special option, "H2vJ" */
		else if( strncmp(chLABEL_CAPS , "H2" , 2 ) == 0 &&
					isdigit(chLABEL_CAPS[2]) && isdigit(chLABEL_CAPS[3]) )
		{
			long int iVib = chLABEL_CAPS[2] - '0';
			long int iRot = chLABEL_CAPS[3] - '0';
			if( iVib<0 || iRot < 0 )
			{
				fprintf( ioQQQ, " PROBLEM cdColm called with insane v,J for H2=\"%s\" caps=\"%s\"\n",
				  chLabel , chLABEL_CAPS );
				return 1;
			}
			*theocl = cdH2_colden( iVib , iRot );
		}

		/* ===========================================================*/
		/* end special case molecular column densities, start special cases
		 * excited state column densities */
		else if( strcmp( chLABEL_CAPS , "HE1*" )==0 )
		{
			if (dense.lgElmtOn[ipHELIUM])
			{
				gs.q = iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2s3S];
			}
			*theocl = column(gs);
		}
		/* general case */
		else 
		{
			vector<genericState> v = matchGeneric(chLabel,false);
			*theocl = 0.;
			for (size_t i=0; i<v.size(); ++i)
			{
				*theocl += column(v[i]);
			}
		}
	}
	else
	{
		/* this case, ionization stage of some element */
		/* find which element this is */
		nelem = 0;
		while( nelem < LIMELM && 
			strncmp(chLABEL_CAPS,elementnames.chElementNameShort[nelem],4) != 0 )
		{
			++nelem;
		}

		/* this is true if we have one of the first 30 elements in the label,
		 * nelem is on C scale */
		if( nelem < LIMELM )
		{

			/* sanity check - does this ionization stage exist?
			 * max2 is to pick up H2 as H 3 */
			if( ion > MAX2(3,nelem + 2) )
			{
				fprintf( ioQQQ, 
				  " cdColm asked to return ionization stage %ld for element \"%s\" but this is not physical.\n", 
				  ion, chLabel );
				return 1;
			}

			/* the column density, ion is on physics scale, but means are on C scale */
			*theocl = mean.xIonMean[0][nelem][ion-1][0];
			/*>>chng 06 jan 23, div by factor of two
			 * special case of H2 when being tricked as H 3 - this stores 2H_2 so that
			 * the fraction of H in H0 and H+ is correct - need to remove this extra
			 * factor of two here */
			if( nelem==ipHYDROGEN && ion==3 )
				*theocl /= 2.;
		}
		else
		{
			fprintf( ioQQQ, 
			  " cdColm did not understand this combination of label \"%s\" and ion %4ld.\n", 
						chLabel, ion );
			return 1;
		}
	}
	return 0;
}


/*************************************************************************
 *
 *cdErrors produce summary of all warnings, cautions, etc, on arbitrary io unit 
 *
 ************************************************************************/

void cdErrors(FILE *ioOUT)
{
	long int nc, 
	  nn, 
	  npe, 
	  ns, 
	  nte, 
	  nw ,
	  nIone,
	  nEdene;
	bool lgAbort_loc;

	DEBUG_ENTRY( "cdErrors()" );

	/* first check for number of warnings, cautions, etc */
	cdNwcns(&lgAbort_loc,&nw,&nc,&nn,&ns,&nte,&npe, &nIone, &nEdene );

	/* only say something is one of these problems is nonzero */
	if( nw || nc || nte || npe ||	nIone || nEdene || lgAbort_loc )
	{
		/* say the title of the model */
		fprintf( ioOUT, "%75.75s\n", input.chTitle );

		if( lgAbort_loc )
			fprintf(ioOUT," Calculation ended with abort!\n");

		/* print warnings on the io unit */
		if( nw != 0 )
		{
			cdWarnings(ioOUT);
		}

		/* print cautions on the io unit */
		if( nc != 0 )
		{
			cdCautions(ioOUT);
		}

		if( nte != 0 )
		{
			fprintf( ioOUT , "Te failures=%4ld\n", nte );
		}

		if( npe != 0 )
		{
			fprintf( ioOUT , "Pressure failures=%4ld\n", npe );
		}

		if( nIone != 0 )
		{
			fprintf( ioOUT , "Ionization failures=%4ld\n", nte );
		}

		if( nEdene != 0 )
		{
			fprintf( ioOUT , "Electron density failures=%4ld\n", npe );
		}
	}
	return;
}

/*************************************************************************
 *
 *cdDepth_depth get depth structure from previous iteration 
 *
 ************************************************************************/
void cdDepth_depth( double cdDepth[] )
{
	long int nz;

	DEBUG_ENTRY( "cdDepth_depth()" );

	for( nz = 0; nz<nzone; ++nz )
	{
		cdDepth[nz] = struc.depth[nz];
	}
	return;
}

/*************************************************************************
 *
 *cdPressure_depth routine to query results and return pressure of last iteration 
 *
 ************************************************************************/

/*
 * cdPressure_depth
 * This returns the pressure and its constituents for the last iteration. 
 * space was allocated in the calling routine for the vectors - 
 * before calling this, cdnZone should have been called to get the number of
 * zones, then space allocated for the arrays */
void cdPressure_depth(
	/* total pressure, all forms*/
	double TotalPressure[],			
	/* gas pressure */
	double GasPressure[],				
	/* line radiation pressure */
	double RadiationPressure[])
{
	long int nz;

	DEBUG_ENTRY( "cdPressure_depth()" );

	for( nz = 0; nz<nzone; ++nz )
	{
		TotalPressure[nz] = struc.pressure[nz];
		GasPressure[nz] = struc.GasPressure[nz];
		RadiationPressure[nz] = struc.pres_radiation_lines_curr[nz];
	}
	return;
}

/*************************************************************************
 *
 *cdPressure_last routine to query results and return pressure of last zone 
 *
 ************************************************************************/

void cdPressure_last(
		double *PresTotal,  /* total pressure, all forms, for the last computed zone*/
		double *PresGas,    /* gas pressure */
		double *PresRad)    /* line radiation pressure */
{

	DEBUG_ENTRY( "cdPressure_last()" );

	*PresGas = pressure.PresGasCurr;
	*PresRad = pressure.pres_radiation_lines_curr;
	*PresTotal = pressure.PresTotlCurr;
	return;
}

/*************************************************************************
 *
 *cdnZone gets number of zones
 *
 ************************************************************************/

/* returns number of zones */
long int cdnZone() 
{
	return nzone;
}

/*************************************************************************
 *
 *cdTemp_last routine to query results and return temperature of last zone 
 *
 ************************************************************************/


double cdTemp_last()
{
	return phycon.te;
}

/*************************************************************************
 *
 *cdIonFrac get ionization fractions for a constituent
 *
 ************************************************************************/

int cdIonFrac(
	/* four char string, null terminated, giving the element name */
	const char *chLabel, 
	/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
	 * 0 says special case */
	long int IonStage, 
	/* will be fractional ionization */
	double *fracin, 
	/* how to weight the average, must be "VOLUME" or "RADIUS" */
	const char *chWeight ,
	/* if true then weighting also has electron density, if false then only volume or radius */
	bool lgDensity ) 
	/* return value is 0 if element was found, non-zero if failed */
{
	long int ip, 
		ion, /* used as index within aaa vector*/
		nelem;
	realnum aaa[LIMELM + 1];
	/* use following to store local image of character strings */
	char chCARD[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "cdIonFrac()" );

	strcpy( chCARD, chWeight );
	/* make sure chWeight is all caps */
	caps(chCARD);/* convert to caps */

	int dim;
	if( strcmp(chCARD,"RADIUS") == 0 )
		dim = 0;
	else if( strcmp(chCARD,"AREA") == 0 )
		dim = 1;
	else if( strcmp(chCARD,"VOLUME") == 0 )
		dim = 2;
	else
	{
		fprintf( ioQQQ, " cdIonFrac: chWeight=%6.6s makes no sense to me, valid options are RADIUS, AREA, and VOLUME\n", 
		  chWeight );
		*fracin = 0.;
		return 1;
	}

	/* first ensure that chLabel is all caps */
	strcpy( chCARD, chLabel );
	/* make sure chLabel is all caps */
	caps(chCARD);/* convert to caps */
	trimTrailingWhiteSpace(chCARD);

	if( IonStage==0 )
	{
		/* special case */
		if( strcmp(chCARD,"H2" ) == 0 )
		{
			/* this will be request for H2, treated as third stage of hydrogen */
			nelem = 0;
			IonStage = 3;
		}
		else
		{
			fprintf( ioQQQ, " cdIonFrac: ion stage of zero and element %s makes no sense to me\n", 
			chCARD );
			*fracin = 0.;
			return 1;
		}
	}

	else 
	{
		/* find which element this is */
		nelem = 0;
		while( nelem < LIMELM && 
		       strcmp(chCARD,elementnames.chElementNameShort[nelem]) != 0 )
		{
			++nelem;
		}
	}

	/* if element not recognized and above loop falls through, nelem is LIMELM */
	if( nelem >= LIMELM )
	{
		fprintf( ioQQQ, " cdIonFrac called with unknown element chLabel, =%4.4s\n", 
		  chLabel );
		return 1;
	}

	/* sanity check - does this ionization stage exist? 
	 * IonStage is on spectroscopic scale and nelem is on C scale */

	/* ion will be used as pointer within the aaa array that contains average values,
	 * convert to C scale */
	ion = IonStage - 1;

	if( (ion > nelem+1 || ion < 0 ) && !(nelem==ipHYDROGEN&&ion==2))
	{
		fprintf( ioQQQ, " cdIonFrac asked to return ionization stage %ld for element %4.4s but this is not physical.\n", 
		  IonStage, chLabel );
		*fracin = -1.;
		return 1;
	}

	/* get average, aaa is filled in from 0 */
	/* aaa is dim'd LIMELM+1 so largest argument is LIMELM
	 * 'i' means ionization, not temperature */
	/* last argument says whether to include electron density */
	/* MeanIon uses c scale for nelem */
	mean.MeanIon('i',nelem,dim,&ip,aaa,lgDensity);
	*fracin = exp10((double)aaa[ion]);

	/* we succeeded - say so */
	return 0;
}

/*************************************************************************
 *
 * debugLine provides a debugging hook into the main line array 
 *
 ************************************************************************/

 /* debugLine provides a debugging hook into the main line array 
  * loops over whole array and finds every line that matches length,
  * the wavelength, the argument to the function
  * put breakpoint inside if test 
  * return value is number of matches, also prints all matches*/
long debugLine( realnum wavelength )
{
	long j, kount;
	realnum errorwave;

	kount = 0;

	/* get the error associated with specified significant figures */
	errorwave = WavlenErrorGet( wavelength, LineSave.sig_figs );

	for( j=0; j < LineSave.nsum; j++ )
	{
		/* check wavelength and chLabel for a match */
		/* if( fabs(LineSave.lines[j].wavelength- wavelength)/MAX2(DELTA,wavelength) < errorwave ) */
		if( fabs(LineSave.lines[j].wavelength()-wavelength) < errorwave )
		{
			LineSave.lines[j].prt(stdout);
			printf("\n");
			++kount;
		}
	}
	printf(" hits = %li\n", kount );
	return kount;
}

/*************************************************************************
 *
 *cdLine get the predicted line intensity, also index for line in stack 
 *
 ************************************************************************/

// returns array index for line in array stack if we found the line,
// return negative of total number of lines as debugging aid if line not found
// return <0 if problems
// emergent or intrinsic not specified - use intrinsic
long int cdLine(
		const char *chLabel, 
		/* wavelength of line in angstroms, not format printed by code */
		realnum wavelength, 
		/* linear intensity relative to normalization line*/
		double *relint, 
		/* log of luminosity or intensity of line */
		double *absint )
{
	DEBUG_ENTRY( "cdLine()" );
	long int i = cdLine( chLabel , wavelength , relint , absint, 0 );
	return i;
}
long int cdLine(
	  const char *chLabel, 
		/* wavelength of line in angstroms, not format printed by code */
	  realnum wavelength, 
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* log of luminosity or intensity of line */
	  double *absint ,
	  // 0 is intrinsic,
	  // 1 emergent
	  // 2 is intrinsic cumulative,
	  // 3 emergent cumulative
	  int LineType )
{
	DEBUG_ENTRY( "cdLine()" );

	*absint = 0.;
	*relint = 0.;

	long ipobs = LineSave.findline(chLabel, wavelength);
	if (ipobs > 0)
		cdLine_ip(ipobs,relint,absint,LineType);
	//	printf("%s\t %g\t ip = %ld\t absint = %g\n", chLabel, wavelength, ipobs, *absint);

	return ipobs;
}



/*cdLine_ip get the predicted line intensity, using index for line in stack */
void cdLine_ip(long int ipLine, 
			   /* linear intensity relative to normalization line*/
			   double *relint, 
			   /* log of luminosity or intensity of line */
			   double *absint )
{

	DEBUG_ENTRY( "cdLine_ip()" );
	cdLine_ip( ipLine , relint , absint , 0 );
}
void cdLine_ip(long int ipLine, 
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* log of luminosity or intensity of line */
	  double *absint ,
	  // 0 is intrinsic,
	  // 1 emergent
	  // 2 is intrinsic cumulative,
	  // 3 emergent cumulative
	  int LineType )
{

	DEBUG_ENTRY( "cdLine_ip()" );

	*relint = 0.;
	*absint = 0.;

	if( LineType<0 || LineType>3 )
	{
		fprintf( ioQQQ, " PROBLEM cdLine_ip called with insane nLineType - it must be between 0 and 3.\n");
		return;
	}

	/* this is zero when cdLine called with cdNoExec called too */
	if( LineSave.nsum == 0 )
	{
		return;
	}
	ASSERT( LineSave.ipNormWavL >= 0);
	ASSERT( LineSave.nsum > 0);

	/* does the normalization line have a positive intensity */
	if( LineSave.lines[LineSave.ipNormWavL].SumLine(LineType) > 0. )
		*relint = LineSave.lines[ipLine].SumLine(LineType)/
			LineSave.lines[LineSave.ipNormWavL].SumLine(LineType)*
			LineSave.ScaleNormLine;

	/* return log of current line intensity if it is positive */
	if( LineSave.lines[ipLine].SumLine(LineType) > 0. )
		*absint = LineSave.lines[ipLine].SumLine(LineType) *
			radius.Conv2PrtInten;

	return;
}



/*************************************************************************
 *
 *cdNwcns get the number of cautions and warnings, to tell if calculation is ok
 *
 ************************************************************************/

void cdNwcns(
	/* abort status, this better be false */
	bool *lgAbort_ret ,
	/* the number of warnings, cautions, notes, and surprises */
	long int *NumberWarnings, 
	long int *NumberCautions, 
	long int *NumberNotes, 
	long int *NumberSurprises, 
	/* the number of temperature convergence failures */
	long int *NumberTempFailures, 
	/* the number of pressure convergence failures */
	long int *NumberPresFailures,
	/* the number of ionization convergence failures */
	long int *NumberIonFailures, 
	/* the number of electron density convergence failures */
	long int *NumberNeFailures )
{

	DEBUG_ENTRY( "cdNwcns()" );

	/* this would be set true if code ended with abort - very very bad */
	*lgAbort_ret = lgAbort;
	/* this sub is called after comment, to find the number of various comments */
	*NumberWarnings = warnings.nwarn;
	*NumberCautions = warnings.ncaun;
	*NumberNotes = warnings.nnote;
	*NumberSurprises = warnings.nbang;

	/* these are counters that were incremented during convergence failures */
	*NumberTempFailures = conv.nTeFail;
	*NumberPresFailures = conv.nPreFail;
	*NumberIonFailures = conv.nIonFail;
	*NumberNeFailures = conv.nNeFail;
	return;
}

void cdOutput( const string& filename, const char *mode )
{
	DEBUG_ENTRY( "cdOutput()" );

	if( ioQQQ != stdout && ioQQQ != NULL )
		fclose(ioQQQ);
	FILE *fp = stdout;
	if( !filename.empty() )
		fp = open_data( filename.c_str(), mode, AS_LOCAL_ONLY );
	save.chOutputFile = filename;
	ioQQQ = fp;
}

void cdOutput( const string& filename, FILE* fp )
{
	DEBUG_ENTRY( "cdOutput()" );

	if( ioQQQ != stdout && ioQQQ != NULL )
		fclose(ioQQQ);
	save.chOutputFile = filename;
	ioQQQ = ( fp != NULL ) ? fp : stdout;
}

void cdInput( const string& filename, const char *mode )
{
	DEBUG_ENTRY( "cdInput()" );

	if( ioStdin != stdin && ioStdin != NULL )
		fclose(ioStdin);
	FILE *fp = stdin;
	if( !filename.empty() )
	{
		fp = open_data( filename.c_str(), mode, AS_LOCAL_ONLY_TRY );
		if( fp == NULL )
		{
			fprintf( ioQQQ, " input file \"%s\" not found\n", filename.c_str() );
			cdEXIT(ES_FAILURE);
		}
	}
	ioStdin = fp;
}

/*************************************************************************
 *
 * cdTalk tells the code whether to print results or be silent
 *
 ************************************************************************/

void cdTalk(bool lgTOn)
{

	DEBUG_ENTRY( "cdTalk()" );

	/* MPI talking rules must be obeyed, otherwise loss of output may result */
	called.lgTalk = lgTOn && cpu.i().lgMPI_talk();
	/* has talk been forced off? */
	called.lgTalkForcedOff = !lgTOn;
	return;
}

/*cdB21cm - returns B as measured by 21 cm */
double cdB21cm()
{
	double ret;

	DEBUG_ENTRY( "cdB21cm()" );

	// return average over radius
	if( mean.TempB_HarMean[0][1] > SMALLFLOAT )
	{
		ret = mean.TempB_HarMean[0][0]/mean.TempB_HarMean[0][1];
	}
	else
	{
		ret = 0.;
	}
	return ret;
}

/*************************************************************************
 *
 * cdTemp get mean electron temperature for any element 
 *
 ************************************************************************/

/* This routine finds the mean electron temperature for any ionization stage 
 * It returns 0 if it could find the species, 1 if it could not find the species.
 * The first argument is a null terminated 4 char string that gives the element
 * name as spelled by Cloudy.  
 * The second argument is ion stage, 1 for atom, 2 for first ion, etc
 * This third argument will be returned as result,
 * Last parameter is either "RADIUS", "AREA", or "VOLUME" to give weighting 
 *
 * if the ion stage is zero then the element label will have a special meaning.
 * The string "21CM" will return the 21 cm temperature.
 * The string "H2  " will return the temperature weighted wrt the H2 density
 * There are several other options as listed below */
/** \todo	2	this should have a last argument like cdIonFrac for whether or not weighting
 * is wrt electron density */

/* return value is o if things are ok and element was found, 
 * non-zero if element not found or there are problems */
int cdTemp(
	/* four char string, null terminated, giving the element name */
	const char *chLabel, 
	/* IonStage is ionization stage, 1 for atom, up to Z+1 where Z is atomic number,
	 * 0 means that chLabel is a special case */
	long int IonStage, 
	/* will be temperature */
	double *TeMean, 
	/* how to weight the average, must be "RADIUS", "AREA, or "VOLUME" */
	const char *chWeight ) 
{
	long int ip, 
		ion, /* used as pointer within aaa vector*/
		nelem;
	realnum aaa[LIMELM + 1];
	/* use following to store local image of character strings */
	char chWGHT[INPUT_LINE_LENGTH] , chELEM[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "cdTemp()" );

	/* make sure strings are all caps */
	strcpy( chWGHT, chWeight );
	caps(chWGHT);
	strcpy( chELEM, chLabel );
	caps(chELEM);
	trimTrailingWhiteSpace(chELEM);

	/* now see which weighting */
	int dim;
	if( strcmp(chWGHT,"RADIUS") == 0 )
		dim = 0;
	else if( strcmp(chWGHT,"AREA") == 0 )
		dim = 1;
	else if( strcmp(chWGHT,"VOLUME") == 0 )
		dim = 2;
	else
	{
		fprintf( ioQQQ, " cdTemp: chWeight=%6.6s makes no sense to me, the options are RADIUS, AREA, and VOLUME.\n", 
		  chWeight );
		*TeMean = 0.;
		return 1;
	}

	if( IonStage == 0 )
	{
		/* return atomic hydrogen weighted harmonic mean of gas kinetic temperature */
		if( strcmp(chELEM,"21CM") == 0 )
		{
			if( mean.TempHarMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempHarMean[dim][0]/mean.TempHarMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* return atomic hydrogen weighted harmonic mean 21 cm spin temperature,
		 * using actual level populations with 1s of H0 */
		else if( strcmp(chELEM,"SPIN") == 0 )
		{
			if( mean.TempH_21cmSpinMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempH_21cmSpinMean[dim][0] / mean.TempH_21cmSpinMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* return temperature deduced from ratio of 21 cm and Lya optical depths */
		else if( strcmp(chELEM,"OPTI") == 0 )
		{
			if( HFLines[0].Emis().TauCon() > SMALLFLOAT )
				*TeMean = 3.84e-7 * iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauCon() /
					HFLines[0].Emis().TauCon();
			else
				*TeMean = 0.;
		}
		/* mean temp weighted by H2 density */
		else if( strcmp(chELEM,"H2") == 0 )
		{
			if( mean.TempIonMean[dim][ipHYDROGEN][2][1] > SMALLFLOAT )
				*TeMean = mean.TempIonMean[dim][ipHYDROGEN][2][0] /
					mean.TempIonMean[dim][ipHYDROGEN][2][1];
			else
				*TeMean = 0.;
		}
		/* this is temperature weighted by eden */
		else if( strcmp(chELEM,"TENE") == 0 )
		{
			if( mean.TempEdenMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempEdenMean[dim][0]/mean.TempEdenMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* four spaces mean to return simple mean of Te */
		else if( strcmp(chELEM,"") == 0 )
		{
			if( mean.TempMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempMean[dim][0]/mean.TempMean[dim][1];
			else
				*TeMean = 0.;
		}
		else
		{
			fprintf( ioQQQ, " cdTemp called with ion=0 and unknown quantity, =%4.4s\n", 
			  chLabel );
			return 1;
		}

		/* say things are ok */
		return 0;
	}

	/* find which element this is */
	nelem = 0;
	while( nelem < LIMELM && 
	       strcmp(chELEM,elementnames.chElementNameShort[nelem]) != 0 )
	{
		++nelem;
	}

	/* if element not recognized and above loop falls through, nelem is LIMELM */
	if( nelem >= LIMELM )
	{
		fprintf( ioQQQ, " cdTemp called with unknown element chLabel, =%4.4s\n", 
		  chLabel );
		return 1;
	}

	/* sanity check - does this ionization stage exist? 
	 * IonStage on spectroscopic scale, nelem on c */

	/* ion will be used as pointer within the aaa array that contains average values,
	 * done this way to prevent lint from false problem in access to aaa array */
	ion = IonStage - 1;

	if( ion > nelem+1 || ion < 0 )
	{
		fprintf( ioQQQ, " cdTemp asked to return ionization stage %ld for element %4.4s but this is not physical.\n", 
		  IonStage, chLabel );
		return 1;
	}

	/* get average, aaa is filled in from 0 */
	/* aaa is dim'd LIMELM+1 so largest arg is LIMELM */
	/* MeanIon uses C scale for nelem */
	mean.MeanIon('t', nelem,dim,&ip,aaa,false);
	*TeMean = exp10((double)aaa[ion]);
	return 0;
}

/*************************************************************************
 *
 * cdRead routine to read in command lines 
 * called by user when cloudy used as subroutine 
 * called by cdMain and ParseInit when used as standalone program
 *
 ************************************************************************/

/* returns the number of lines available in command stack 
 * this is limit to how many more commands can be read */
int cdRead( const char *chInputLine ) /* the string containing the command */
{
	const bool COLUMN0 = true;
	const bool KEEP_VISIBLE = false;
	const bool STRIP_VISIBLE = true;

	DEBUG_ENTRY( "cdRead()" );

	/* this is set false when code loaded, set true when cdInit called,
	 * this is check that cdInit was called first */
	if( !lgcdInitCalled )
	{
		printf(" cdInit was not called first - this must be the first call.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* totally ignore hidden comment lines but want to include visible
	 * comments since they are copied to output */
	if( lgIsCommentSeq( chInputLine, COLUMN0, KEEP_VISIBLE ) )
	{
		/* return value is number of lines that can still be stuffed in */
		return NKRD - input.nSave;
	}

	/***************************************************************************
	* validate a location to store this line image, then store the version     *
	* that has been truncated from special end of line characters              *
	* stored image will have < INPUT_LINE_LENGTH valid characters              *
	***************************************************************************/

	/* this will now point to the next free slot in the line image save array 
	 * this is where we will stuff this line image */
	++input.nSave;

	if( input.nSave >= NKRD )
	{
		/* too many input commands were entered - bail */
		fprintf( ioQQQ, " Too many line images entered to cdRead.  The limit is %d.\n", NKRD );
		fprintf( ioQQQ, " This limit is set by the variable NKRD which appears in input.h\n" );
		cdEXIT(EXIT_FAILURE);
	}

	string chLocal( chInputLine );

	/* this indicates an overlong input string */
	if( chLocal.length() >= INPUT_LINE_LENGTH )
	{
		chLocal.erase(INPUT_LINE_LENGTH-1);
		fprintf(ioQQQ," PROBLEM cdRead, while parsing the following input line:\n %s\n",
			chInputLine);
		fprintf(ioQQQ," found that the line is longer than %i characters, the longest possible line.\n",
			INPUT_LINE_LENGTH-1);
		fprintf(ioQQQ," Please make the command line no longer than this limit.\n");
	}

	/* now kill hidden comments
	 * this stops info user wants ignored from entering after here
	 * also remove eol character and possible underscores, etc */
	StripComment( chLocal, KEEP_VISIBLE );

	/* save string in master array for later use in output */
	strcpy( input.chCardSav[input.nSave], chLocal.c_str() );

	/* now kill visible comment part of the input line */
	StripComment( chLocal, STRIP_VISIBLE );

	/* pad with two spaces to avoid spurious problems with parsing */
	size_t npad = min(2,INPUT_LINE_LENGTH-1-chLocal.length());
	chLocal.append( npad, ' ' );

	/* save string in master array for later use in parser */
	strcpy( input.chCardStrip[input.nSave], chLocal.c_str() );

	/* keep track of include level of this line */
	input.InclLevel[input.nSave] = input.curInclLevel;

	Parser p;
	p.setline(chLocal.c_str());

	/* check whether this is a trace command, turn on printout if so */
	if( p.hasCommand("TRACE") )
		trace.lgTrace = true;

	if( p.hasCommand("NO ") && p.nMatch("VARY") )
		optimize.lgNoVary = true;

	/* prt.lgPrintTime can be used by check_data() before regular parsing starts */
	if( p.hasCommand("NO ") && p.nMatch("TIME") )
		prt.lgPrintTime = false;

	if( p.hasCommand("GRID") )
	{
		optimize.lgOptimr = true;
		grid.lgGrid = true;
		++grid.nGridCommands;
	}

	/* NO BUFFERING turn off buffered io for standard output, 
	 * used to get complete output when code crashes */
	if( p.hasCommand("NO ") && p.nMatch("BUFF") )
	{
		/* if output has already been redirected (e.g. using cdOutput()) then
		 * ignore this command, a warning will be printed in ParseDont() */
		/* stdout may be a preprocessor macro, so lets be really careful here */
		FILE *test = stdout;
		if( ioQQQ == test )
		{
			fprintf( ioQQQ, " cdRead found NO BUFFERING command, redirecting output to stderr now.\n" );
			/* make sure output is not lost */
			fflush( ioQQQ );
			/* stderr is always unbuffered */
			/* don't use setvbuf() here since we cannot guarantee
			 * that no operations have been performed on stdout */
			ioQQQ = stderr;
			/* will be used to generate comment at end */
			input.lgSetNoBuffering = true;
		}
		else if( !save.chOutputFile.empty() )
		{
			fprintf( ioQQQ, " cdRead found NO BUFFERING command, reopening file %s now.\n",
				 save.chOutputFile.c_str() );
			fclose( ioQQQ );
			// using AS_SILENT_TRY assures that open_data will not write to ioQQQ
			ioQQQ = open_data( save.chOutputFile.c_str(), "a", AS_SILENT_TRY );
			if( ioQQQ == NULL )
			{
				// ioQQQ is no longer valid, so switch to stderr and abort...
				ioQQQ = stderr;
				fprintf( ioQQQ, " cdRead failed to reopen %s, aborting!\n",
					 save.chOutputFile.c_str() );
				cdEXIT(EXIT_FAILURE);
			}
			if( setvbuf( ioQQQ, NULL, _IONBF, 0 ) != 0 )
				fprintf( ioQQQ, " PROBLEM -- cdRead failed to set NO BUFFERING mode.\n" );
			else
				input.lgSetNoBuffering = true;
		}
	}

	/* optimize command read in */
	if( p.hasCommand("OPTI") )
		optimize.lgOptimr = true;

	if( p.hasCommand("SET ") && p.nMatch("CONT") && p.nMatch("RESO") )
	{
		double factor = p.getNumberDefault("",1.);
		/* negative numbers are always logs */
		if( factor <= 0. )
			factor = exp10(factor);
		rfield.setResolutionScaleFactor(factor);
	}

	/* if the command is an init command, process it immediately
	 * so that the lines in the init file can be inserted in situ */
	if( p.hasCommand("INIT") )
	{
		ParseInit(p);
		input.lgInitPresent = true;
	}

	/* here we check for specific keywords, avoiding the TITLE command since
	 * that contains free text, which could lead to spurious matches... */
	if( ! p.hasCommand("TITL") )
	{
		/* now check whether VARY is specified */
		if( p.nMatch("VARY") )
			/* a optimize flag was on the line image */
			optimize.lgVaryOn = true;
	}

	/* print input lines if trace specified */
	if( trace.lgTrace )
		fprintf( ioQQQ,"cdRead [%d] =%s=\n",input.curInclLevel,input.chCardSav[input.nSave] );

	return NKRD - input.nSave;
}

/* wrapper to close all save files */
void cdClosePunchFiles()
{

	DEBUG_ENTRY( "cdClosePunchFiles()" );

	CloseSaveFiles( true );
	return;
}
