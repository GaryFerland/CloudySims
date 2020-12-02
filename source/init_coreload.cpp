/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitCoreload one time initialization of core load, called from cdDrive, this sets
* minimum set of values needed for the code to start - called after
* input lines have been read in and checked for VARY or GRID - so
* known whether single or multiple sims will be run */
#include "cddefines.h"
#include "init.h"
#include "parse.h"
#include "struc.h"
#include "atmdat.h"
#include "hcmap.h"
#include "h2.h"
#include "hextra.h"
#include "heavy.h"
#include "ionbal.h"
#include "iso.h"
#include "taulines.h"
#include "cosmology.h"
#include "broke.h"
#include "dense.h"
#include "optimize.h"

// //////////////////////////////////////////////////////////////////////////
//
//
// NB DO NOT ADD VARIABLES TO THIS FILE!  THE GOAL IS TO REMOVE THIS FILE
// initialization of variables done one time per coreload should be done in
// a constructor for the data
//
//
// //////////////////////////////////////////////////////////////////////////

/*InitCoreload one time initialization of core load, called from cdDrive, this sets
* minimum set of values needed for the code to start - called after
* input lines have been read in and checked for VARY or GRID - so
* known whether single or multiple sims will be run */
void InitCoreload( void )
{
	static int nCalled=0;
	long int nelem;

	DEBUG_ENTRY( "InitCoreload()" );

	/* return if already called */
	if( nCalled )
		return;

	++nCalled;

	hcmap.lgMapOK = true;
	hcmap.lgMapDone = false;

	/* will be reset to positive number when map actually done */
	hcmap.nMapAlloc = 0;
	hcmap.nmap = 0;
	hcmap.lgMapBeingDone = false;

	/* name of output file from optimization run */
	strncpy( chOptimFileName , "optimal.in" , sizeof( chOptimFileName ) );

	/* number of electrons in valence shell that can Compton recoil ionize */
	long int nCom[LIMELM] = 
	{
		1 , 2 ,								/* K 1s shell */
		1 , 2 ,								/* L 2s shell */
		1 , 2 , 3 , 4 , 5 , 6 ,				/* L 2p shell */
		1 , 2 ,								/* M 3s shell */
		1 , 2 , 3 , 4 , 5 , 6 ,				/* M 3p shell */
		1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 ,		/* M 3d shell */
		1 , 2 ,								/* N 4s shell */
		1 , 2								/* N 4p shell */
	};

	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		ionbal.nCompRecoilElec[nelem] = nCom[nelem];
	}

	/* list of shells, 1 to 7 */
	strcpy( Heavy.chShell[0], "1s" );
	strcpy( Heavy.chShell[1], "2s" );
	strcpy( Heavy.chShell[2], "2p" );
	strcpy( Heavy.chShell[3], "3s" );
	strcpy( Heavy.chShell[4], "3p" );
	strcpy( Heavy.chShell[5], "3d" );
	strcpy( Heavy.chShell[6], "4s" );
	
	/* variables for H-like sequence */
	/* default number of levels for hydrogen iso sequence */
	for( nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		iso_sp[ipH_LIKE][nelem].n_HighestResolved_max = 5;
		iso_sp[ipH_LIKE][nelem].nCollapsed_max = 2;
	}

	/* add more collapsed levels for these abundant elements.  Tests
	 * show collapsed levels have very little impact on time */
	iso_sp[ipH_LIKE][ipCARBON].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipNITROGEN].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipOXYGEN].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipNEON].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipSILICON].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipMAGNESIUM].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipSULPHUR].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipIRON].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][LIMELM-1].nCollapsed_max = 5;

	/* H and He are special cases since very high resolution, S/N
	 * optical / IR data are common
	 */
	iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max = 10;
	iso_sp[ipH_LIKE][ipHYDROGEN].nCollapsed_max = 15;

	iso_sp[ipH_LIKE][ipHELIUM].n_HighestResolved_max = 10;
	iso_sp[ipH_LIKE][ipHELIUM].nCollapsed_max = 15;

	/* variables for He-like sequence */
	/* "he-like" hydrogen (H-) is treated elsewhere */
	iso_sp[ipHE_LIKE][ipHYDROGEN].n_HighestResolved_max = -LONG_MAX;
	iso_sp[ipHE_LIKE][ipHYDROGEN].numLevels_max = -LONG_MAX;
	iso_sp[ipHE_LIKE][ipHYDROGEN].nCollapsed_max = -LONG_MAX;

	for( nelem=ipHELIUM; nelem < LIMELM; ++nelem )
	{
		/* put at least three resolved and 1 collapsed in every element for he-like
		 * An n shell for He has twice the number of resolved levels due to singlet/triplet */
		iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max = 3;
		iso_sp[ipHE_LIKE][nelem].nCollapsed_max = 2;
	}

	/* And n=5 for these because they are most abundant */
	iso_sp[ipHE_LIKE][ipCARBON].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipCARBON].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipNITROGEN].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipNITROGEN].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipOXYGEN].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipOXYGEN].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipNEON].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipNEON].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipSILICON].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipSILICON].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipMAGNESIUM].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipMAGNESIUM].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipSULPHUR].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipSULPHUR].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipIRON].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipIRON].nCollapsed_max = 5;
	/* also set this, for exercising any possible issues with highest charge models */
	iso_sp[ipHE_LIKE][LIMELM-1].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][LIMELM-1].nCollapsed_max = 5;

	/* He I is very special case - we want to do a good job, lots of levels */
	iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max = 6;
	iso_sp[ipHE_LIKE][ipHELIUM].nCollapsed_max = 20;

	iso_ctrl.chISO[ipH_LIKE] = "H-like ";
	iso_ctrl.chISO[ipHE_LIKE] = "He-like";

	max_num_levels = 0;
	for( long ipISO = ipH_LIKE; ipISO < NISO; ipISO++ )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* set this to LONG_MAX, reduce to actual number later,
			 * then verify number of levels is not increased after initial coreload */
			iso_sp[ipISO][nelem].numLevels_malloc = LONG_MAX;
			iso_update_num_levels( ipISO, nelem );
		}
	}

	lgStatesAdded = false;
	lgLinesAdded = false;

	/* turn element on if this is first call to this routine,
	* but do not reset otherwise since would clobber how elements were set up */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		/* never turn on element if turned off on first iteration */
		dense.lgElmtOn[nelem] = true;
		/* no elements yet explicitly turned off */
		dense.lgElmtSetOff[nelem] = false;

		/* option to set ionization with element ioniz cmnd,
		* default (false) is to solve for ionization */
		dense.lgSetIoniz[nelem] = false;
		
		// are we using Chianti for this ion?
		for( int ion=0; ion<=nelem+1; ++ion )
		{
			dense.lgIonChiantiOn[nelem][ion] = false;
			dense.lgIonStoutOn[nelem][ion] = false;
			dense.maxWN[nelem][ion] = 0.;
		}
	}

	//set all sources to blank and false.
	//Sources will be added from masterlist files in species.cpp
	for( int nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		for( int ion=0; ion<nelem+2; ++ion )
		{
			strcpy(atmdat.chdBaseSources[nelem][ion]," ");
			atmdat.lgdBaseSourceExists[nelem][ion] = false;
		}
	}

	/** set atmdat.ipSpecIon to uninit flag **/
	for( int nelem = 0; nelem < LIMELM; nelem++)
	{
		vector<long> row;
		for(int ion = 0; ion < nelem+1; ion++)
		{
			row.push_back(-1);
		}
		atmdat.ipSpecIon.push_back(row);
	}

	/* smallest density to permit in any ion - if smaller then set to zero */
	dense.density_low_limit = SMALLFLOAT * 1e3;
	dense.density_low_limit = MAX2( dense.density_low_limit, 1e-50 );

	/* default cosmic ray background */
	/* >>chng 99 jun 24, slight change to value
	* quoted by 
	* >>refer	cosmic ray	ionization rate	McKee, C.M., 1999, astro-ph 9901370
	* this will produce a total
	* secondary ionization rate of 2.5e-17 s^-1, as tested in 
	* test suite cosmicray.in.  If each ionization produces 2.4 eV of heat,
	* the background heating rate should be 9.6e-29 * n*/
	/* >>chng 04 jan 26, update cosmic ray ionization rate for H0 to
	* >>refer	cosmic ray	ionization	Williams, J.P., Bergin, E.A., Caselli, P., 
	* >>refercon	Myers, P.C., & Plume, R. 1998, ApJ, 503, 689,
	* H0 ionization rate of 2.5e-17 s-1 and a H2 ionization rate twice this
	* >>chng 04 mar 15, comment said 2.5e-17 which is correct, but code produce 8e-17,
	* fix back to correct value 
	*/
	/* NB - the rate is derived from the density.  these two are related by the secondary
	* ionization efficiency problem.  background_rate is only here to provide the relationship
	* for predominantly neutral gas.  the background_density is the real rate. 
	hextra.background_density = 1.99e-9f;*/
	/* >>chng 05 apr 16, to get proper ionization rate in ism_set_cr_rate, where
	* H is forced to be fully atomic, no molecules, density from 1.99 to 2.15 */
	/* >>chng 02 apr 05, update to
	 * >>refer	cosmic ray 	ionization	Indriolo, N., Geballe, T., Oka, T., & McCall, B.J. 2007, ApJ, 671, 1736
	 */
	hextra.background_density = 2.15e-9f*7.9f;
	hextra.background_rate = 2.5e-17f*7.9f;

	for( long i=0; i < LIMPAR; ++i )
		optimize.lgOptimizeAsLinear[i] = false;

	/* limit on ionization we check for zoning and prtcomment */
	struc.dr_ionfrac_limit = 1e-3f;

	SaveFilesInit();

	diatoms_init();

	/* initialize cosmological information */
	cosmology.redshift_current = 0.f;
	cosmology.redshift_start = 0.f;
	cosmology.redshift_step = 0.f;
	cosmology.omega_baryon = 0.04592f;
	cosmology.omega_rad = 8.23e-5f;
	cosmology.omega_lambda = 0.7299177f;
	cosmology.omega_matter = 0.27f;
	cosmology.omega_k = 0.f;
	/* the Hubble parameter in 100 km/s/Mpc */
	cosmology.h = 0.71f;
	/* the Hubble parameter in km/s/Mpc */
	cosmology.H_0 = 100.f*cosmology.h;
	cosmology.lgDo = false;

	// the code is ok at startup; only init here so that code remains broken
	// in a grid if any single model announces broken
	broke.lgBroke = false;
	broke.lgFixit = false;
	broke.lgCheckit = false;
	broke.lgPrintFixits = false;

	return;
}
