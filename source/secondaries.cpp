/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "secondaries.h"
t_secondaries secondaries;

void t_secondaries::zero()
{
	DEBUG_ENTRY( "t_secondaries::zero()" );
	/**********************************************************************
	 * all parameters having to do with secondary ionization 
	 * by suprathermal electrons 
	 **********************************************************************/
	SetCsupra = 0.;
	lgCSetOn = false;
	lgSecOFF = false;
	SecHIonMax = 0.;

	HeatEfficPrimary = 1.;
	SecIon2PrimaryErg = 0.;
	SecExcitLya2PrimaryErg = 0.;
	x12tot = 0.;
	sec2total = 0.;

	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		for( long ion=0; ion<nelem+1; ++ion )
		{
			/* secondary ionization rate for each species */
			csupra[nelem][ion] = 0.;
			/* the rate of each species relative to H0 */
			csupra_effic[nelem][ion] = 1.f;
		}
	}
	/* this scale factor is from table 10 of Tielens & Hollenbach 1985 */
	csupra_effic[ipHELIUM][0] = 1.08f;
	
}

void t_secondaries::alloc()
{
	/* malloc space for supra[nelem][ion] */
	csupra = (realnum **)MALLOC( (unsigned)LIMELM*sizeof(realnum *) );
	csupra_effic = (realnum **)MALLOC( (unsigned)LIMELM*sizeof(realnum *) );
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		csupra[nelem] = (realnum *)MALLOC( (unsigned)(nelem+1)*sizeof(realnum) );
		csupra_effic[nelem] = (realnum *)MALLOC( (unsigned)(nelem+1)*sizeof(realnum) );
	}
}
