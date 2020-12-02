/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HydroEinstA calculates Einstein A's from  osillator strengths*/
#include "cddefines.h"
#include "hydroeinsta.h"
#include "hydro_bauman.h"
#include "hydrooscilstr.h"
#include "iso.h"

STATIC double hydro_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo );
STATIC double hydro_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo );

double HydroEinstA(long int n1, 
	  long int n2)
{
	long int lower, iupper;
	double EinstA_v, 
	  ryd, 
	  xl, 
	  xmicron, 
	  xu;

	DEBUG_ENTRY( "HydroEinstA()" );
	/* (lower,upper) of Johnson 1972.  */

	/* strictly n -> n' transition probabilities
	 * no attempt to distribute according to l,l' */

	/* sort out the order of upper and lower, so can be called either way */
	lower = MIN2( n1 , n2 );
	iupper = MAX2( n1, n2 );
	if( lower < 1 || lower == iupper )
	{
		fprintf(ioQQQ," HydroEinstA called with impossible ns, =%li %li\n", lower, iupper);
		cdEXIT(EXIT_FAILURE);
	}

	xl = (double)lower;
	xu = (double)iupper;
	ryd = 1./POW2(xl) - 1./POW2(xu);
	xmicron = 1.E4/(ryd*RYD_INF);
	EinstA_v = HydroOscilStr(xl,xu)*TRANS_PROB_CONST*1e8f/(POW2(xmicron))*xl*xl/xu/xu;
	return( EinstA_v );
}

realnum hydro_transprob( long nelem, long ipHi, long ipLo )
{
	double Aul;
	long ipISO = ipH_LIKE;
	/* charge to 4th power, needed for scaling laws for As*/
	double z4 = POW4((double)nelem+1.);
	DEBUG_ENTRY( "hydro_transprob()" );

	if( ipHi >= iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max )
	{
		if( ipLo >= iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max )
		{
			Aul = hydro_transprob_collapsed_to_collapsed( nelem, N_(ipHi), N_(ipLo) );
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);

			ASSERT( Aul > 0.);
		}
		else 
		{
			Aul = hydro_transprob_collapsed_to_resolved( nelem, N_(ipHi), N_(ipLo), L_(ipLo) );
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.01f,0.01f);
		}
	}
	else
	{
		if(  N_(ipHi) == N_(ipLo)  )
		{	
			/** \todo 1 define quantum defects and use scqdri to 
			 * calculate A's if levels are not exactly degenerate. */
			Aul = SMALLFLOAT;
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else if( ipLo == 0 && ipHi == 1 )
		{
			// >> refer	H-like	As	Marrus, E. \& Mohr, P. J. Advances in Atomic and Molecular Physics, Vol. 14, Academic, New York, 1978, p. 181
			Aul =  2.46e-6*powi((double)(nelem+1.),10);
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else if( ipLo == 0 && ipHi == 2 )
		{
			Aul = 6.265e8*z4;
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else if( abs( L_(ipLo) - L_(ipHi) )== 1 )
		{
			Aul = H_Einstein_A( N_(ipHi), L_(ipHi), N_(ipLo), L_(ipLo), nelem+1 );
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else
		{
			ASSERT( N_(ipHi) > N_(ipLo) );
			ASSERT( (L_(ipHi) == L_(ipLo)) || 
				( abs(L_(ipHi)-L_(ipLo)) > 1) );
			Aul = SMALLFLOAT;
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
	}

	return (realnum)Aul;
}

STATIC double hydro_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo )
{
	DEBUG_ENTRY( "hydro_transprob_collapsed_to_collapsed()" );

	double Aul = 0.;
	ASSERT( nHi > nLo );

	for( long lLo=0; lLo < nLo; ++lLo )
	{
		double Aul1 = hydro_transprob_collapsed_to_resolved( nelem, nHi, nLo, lLo );
		Aul += Aul1;
	}

	return Aul;
}

STATIC double hydro_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo )
{
	DEBUG_ENTRY( "hydro_transprob_collapsed_to_resolved()" );
		
	long ipISO = ipH_LIKE;	
	t_iso_sp *sp = &iso_sp[ipISO][nelem];
	long nResMax = sp->n_HighestResolved_max;
	long ipLoRes = sp->IndexIfAllResolved[nLo][lLo][2];
	ASSERT( nLo > nResMax || ipLoRes <= sp->numLevels_max - sp->nCollapsed_max );

	/* Lower level resolved, upper not. First calculate Aul
	 * from upper level with ang mom one higher.	*/
	double Aul = H_Einstein_A( nHi, lLo+1, nLo, lLo, nelem+1 );

	sp->CachedAs[ nHi-nResMax-1 ][ ipLoRes ][0] = (realnum)Aul;

	Aul *= (2.*(lLo+1.)+1.) * 2. / (2.*(double)nHi*(double)nHi);

	if( lLo != 0 )
	{
		/* For all l>0, add in transitions from upper level with ang mom one lower.	*/
		double Aul1 = H_Einstein_A( nHi, lLo-1, nLo, lLo, nelem+1 );
	
		sp->CachedAs[ nHi-nResMax-1 ][ ipLoRes ][1] = (realnum)Aul1;

		/* now add in this part after multiplying by stat weight for lHi = lLo-1.	*/
		Aul += Aul1*(2.*(lLo-1.)+1.) * 2. / (2.*(double)nHi*(double)nHi);
	}
	else
		sp->CachedAs[ nHi-nResMax-1 ][ ipLoRes ][1] = 0.f;

	ASSERT( Aul > 0.);

	return Aul;
}

