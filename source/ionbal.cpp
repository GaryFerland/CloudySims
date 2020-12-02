/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "ionbal.h"
#include "dense.h"
#include "warnings.h"

t_ionbal ionbal;

void t_ionbal::alloc()
{
	/* these will save bound electron recoil information data */
	ipCompRecoil = 
		(long**)MALLOC(sizeof(long*)*(unsigned)LIMELM );
	CompRecoilIonRate = 
		(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
	CompRecoilIonRateSave = 
		(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
	CompRecoilHeatRate = 
		(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
	CompRecoilHeatRateSave = 
		(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
	PhotoRate_Shell = 
		(double****)MALLOC(sizeof(double***)*(unsigned)LIMELM );
	CollIonRate_Ground = 
		(double***)MALLOC(sizeof(double**)*(unsigned)LIMELM );
	ExcitationGround = 
		(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
	UTA_ionize_rate = 
		(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
	UTA_heat_rate = 
		(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
	/* space for ionization recombination arrays */
	RateIoniz = (double ***)MALLOC(sizeof(double **)*(unsigned)LIMELM );
	RateRecomTot = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
	RateRecomIso = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
	RR_rate_coef_used = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
	RR_Verner_rate_coef = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
	
	/* rate coefficients [cm3 s-1] for Badnell DR recombination */
	DR_Badnell_rate_coef = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
	RR_Badnell_rate_coef = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
	CX_recomb_rate_used = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );

	DR_Badnell_suppress_fact.reserve( LIMELM );
	
	/* create arrays for ions */
	for( long nelem=0; nelem<LIMELM; ++nelem )
	{
		DR_Badnell_rate_coef[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		RR_Badnell_rate_coef[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		CX_recomb_rate_used[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		
		RateIoniz[nelem] = (double **)MALLOC(sizeof(double *)*(unsigned)(nelem+1) );
		RateRecomTot[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );					
		
		for( long ion=0; ion<nelem+1; ++ion )
		{
			RateIoniz[nelem][ion] = (double *)MALLOC(sizeof(double )*(unsigned)(nelem+2) );
			for( long ion2=0; ion2<nelem+2; ++ion2 )
				RateIoniz[nelem][ion][ion2] = 0.;
		}
		
		RR_rate_coef_used[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		RR_Verner_rate_coef[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		UTA_ionize_rate[nelem] = 
			(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		UTA_heat_rate[nelem] = 
			(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		ipCompRecoil[nelem] = 
			(long*)MALLOC(sizeof(long)*(unsigned)(nelem+1) );
		CompRecoilIonRate[nelem] = 
			(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		CompRecoilIonRateSave[nelem] = 
			(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		CompRecoilHeatRate[nelem] = 
			(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		CompRecoilHeatRateSave[nelem] = 
			(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		PhotoRate_Shell[nelem] = 
			(double***)MALLOC(sizeof(double**)*(unsigned)(nelem+1) );
		CollIonRate_Ground[nelem] = 
			(double**)MALLOC(sizeof(double*)*(unsigned)(nelem+1) );
		ExcitationGround[nelem] = 
			(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
		RateRecomIso[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(NISO) );
		for( long ipISO=0; ipISO<NISO; ++ipISO )
		{
			RateRecomIso[nelem][ipISO] = 0.;	
		}

		for( long ion=0; ion<nelem+1; ++ion )
		{
			/* >>chng 03 aug 09, set these to impossible values */
			RateRecomTot[nelem][ion] = -1.;
			UTA_ionize_rate[nelem][ion] = -1.;
			UTA_heat_rate[nelem][ion] = -1.;
			ipCompRecoil[nelem][ion] = -99;
			CompRecoilIonRate[nelem][ion] = -1.;
			CompRecoilIonRateSave[nelem][ion] = -1.;
			CompRecoilHeatRate[nelem][ion] = -1.;
			CompRecoilHeatRateSave[nelem][ion] = -1.;
			
			/* finish mallocing space */
			PhotoRate_Shell[nelem][ion] = 
				(double**)MALLOC(sizeof(double*)*(unsigned)NSHELLS );
			CollIonRate_Ground[nelem][ion] = 
				(double*)MALLOC(sizeof(double)*(unsigned)2 );
			for( long ns=0; ns<NSHELLS; ++ns )
			{
				PhotoRate_Shell[nelem][ion][ns] = 
					(double*)MALLOC(sizeof(double)*(unsigned)3 );
			}
			
			/* now set to impossible values */
			ipCompRecoil[nelem][ion] = -100000;
			DR_Badnell_rate_coef[nelem][ion] = 0.;
			RR_Badnell_rate_coef[nelem][ion] = 0.;
		}
		
		set_NaN( RR_rate_coef_used[nelem], nelem+1 );
		set_NaN( RR_Verner_rate_coef[nelem], nelem+1 );
		set_NaN( CX_recomb_rate_used[nelem], nelem+1 );

		DR_Badnell_suppress_fact.reserve( nelem, nelem+1 );
	}

	DR_Badnell_suppress_fact.alloc();
}

void t_ionbal::zero()
{
	/* now zero out these arrays */
	for( long nelem=0; nelem< LIMELM; ++nelem )
	{
		for( long ion=0; ion<nelem+1; ++ion )
		{

			CompRecoilHeatRate[nelem][ion] = 0.;
			CompRecoilIonRate[nelem][ion] = 0.;
			UTA_ionize_rate[nelem][ion] = 0.;
			UTA_heat_rate[nelem][ion] = 0.;
			CollIonRate_Ground[nelem][ion][0] = 0.;
			CollIonRate_Ground[nelem][ion][1] = 0.;
			ExcitationGround[nelem][ion] = 0.;
			RateRecomTot[nelem][ion] = 0.;
			for( long ns=0; ns < NSHELLS; ++ns )
			{
				/* must be zero since ion routines use these when
				 * not yet defined */
				PhotoRate_Shell[nelem][ion][ns][0] = 0.;
				PhotoRate_Shell[nelem][ion][ns][1] = 0.;
				PhotoRate_Shell[nelem][ion][ns][2] = 0.;
			}
		}
	}

	/* limits for highest and lowest stages of ionization in TrimStage */
	lgTrimhiOn = true;
	trimhi = 1e-6;
	lgTrimloOn = true;
	trimlo = 1e-10;
	lgNewTrim = false;
	
	lgPhotoIoniz_On = true;
	lgCompRecoil = true;

	lgDRsup = true;

	lgNoCota = false;
	for( long nelem = 0; nelem < LIMELM; ++nelem )
	{
		CotaRate[nelem] = 0.;
	}
	ilt = 0;
	iltln = 0;
	ilthn = 0;
	ihthn = 0;
	ifail = 0;
	lgGrainIonRecom = true;

	/* option to print recombination coefficient then exit */
	lgRecom_Badnell_print = false;
	guess_noise = 0.;

} 

void t_ionbal::comment(t_warnings& w)
{
	char chLine[INPUT_LINE_LENGTH];

	if( ionbal.CompHeating_Max > 0.05 )
	{
		sprintf( chLine, 
			"  !Bound Compton heating reached %.2f%% of the local heating.", 
		  ionbal.CompHeating_Max*100. );
		w.bangin(chLine);
	}
	else if( ionbal.CompHeating_Max > 0.01 )
	{
		sprintf( chLine, 
			"   Bound Compton heating reached %.2f%% of the local heating.", 
		  ionbal.CompHeating_Max*100. );
		w.notein(chLine);
	}

	/* say if 3-body recombination coefficient function outside range of validity
	 * tripped if te/z**2 < 100 or approx 10**13: => effect >50% of radiative
	 * other integers defined in source for da */
	if( ionbal.ifail > 0 && ionbal.ifail <= 10 )
	{
		sprintf( chLine, 
			"   3 body recombination coefficient outside range %ld", ionbal.ifail );
		w.notein(chLine);
	}
	else if( ionbal.ifail > 10 )
	{
		sprintf( chLine, 
			" C-3 body recombination coefficient outside range %ld", ionbal.ifail );
		w.caunin(chLine);
	}
}

double t_ionbal::RateIonizTot( long nelem, long ion ) const
{
	double sum = 0.;
	
	for( long ion_to=ion+1; ion_to<=dense.IonHigh[nelem]; ion_to++ )
		sum += RateIoniz[nelem][ion][ion_to];
	
	return sum;
}
