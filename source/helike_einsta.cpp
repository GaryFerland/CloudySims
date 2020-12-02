/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*he_1trans compute Aul for given line	*/
/*ritoa - converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
/*ForbiddenAuls calculates transition probabilities for forbidden transitions.	*/
/*scqdri - stands for Semi-Classical Quantum Defect Radial Integral	*/
/*Jint - used by scqdri	*/
/*AngerJ - used by scqdri */
/*DoFSMixing - applies a fine structure mixing approximation to A's.  To be replaced by 
 * method that treats the entire rate matrix.	*/

#include "cddefines.h" 
#include "dense.h"
#include "trace.h"
#include "hydro_bauman.h"
#include "iso.h"
#include "helike.h"
#include "helike_einsta.h"
#include "freebound.h"
#include "lines_service.h"
#include "integrate.h"

/* the array of transitions probabilities read from data file.  */
static double ***TransProbs;

STATIC double helike_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo, double Enerwn );
STATIC double helike_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo, long sLo, long jLo, double Enerwn );

/*ritoa converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
STATIC double ritoa( long li, long lf, long nelem, double k, double RI2 );

/* handles all forbidden transition probabilities. */
STATIC double ForbiddenAuls( long ipHi, long ipLo, long nelem );

/*
static long RI_ipHi, RI_ipLo, RI_ipLev;
static long globalZ;
*/

/* used as parameters in qg32 integration */
static double vJint , zJint;

/* the integrand used in the AngerJ function below. */
STATIC double Jint( double theta )
{
	/*	[ cos[vx - z sin[x]] ]  */
	double 
		d1 = vJint * theta,
		d2 = zJint * sin(theta),
		d3 = (d1 - d2),
		d4 = cos(d3);

	return( d4 );
}

/* AngerJ function. */
// See also http://dlmf.nist.gov/11.10 for series expansion and recurrences
STATIC double AngerJ_mac( double vv, double zz )
{
	// Maclaurin series for J_nu(z)
	// http://dlmf.nist.gov/11.10#iii
	// Series is absolutely convergent, but can take
	// numerous terms for large z
	// Problems when vv=2 -- should fall back to Bessel functions?
	// Worth trying the acceleration techinques of http://dlmf.nist.gov/3.9
	double g1=tgamma(1+0.5*vv), g2=tgamma(1-0.5*vv),
		g3=tgamma(1.5+0.5*vv), g4=tgamma(1.5-0.5*vv);
	double y1 = 0., y2=0., zp = 1.;
	for (long k=0; k<100; ++k)
	{
		double y10 = y1;
		y1 += zp/(g1*g2);
		if (y1 == y10)
			break;
		zp *= 0.5*zz;
		// Gamma(t+1) = t*Gamma(t)
		g1 *= k+0.5*vv+1;
		g2 *= k-0.5*vv+1;
		y2 += zp/(g3*g4);
		zp *= -0.5*zz;
		g3 *= k+0.5*vv+1.5;
		g4 *= k-0.5*vv+1.5;
	}
	double cosfac=cos(0.5*PI*vv), sinfac=sin(0.5*PI*vv);
	return y1*cosfac+y2*sinfac;
}

STATIC double bessjnu( double vv, double zz )
{
	// DLMF 10.17.3
	// Asymptotic series for Bessel function for non-integer index
	double rz=1./zz, akrz=1., bff=0., bgg=0., bfp=0.;
	for (long k=0; k<100; ++k)
	{
		if (k > 0 && fabs(akrz) > fabs(bfp))
			break;
		bfp = akrz;
		bff += akrz;
		akrz *= (4*vv*vv-(4*k+1)*(4*k+1))/(8*(2*k+1)*zz);
		bgg += akrz;
		akrz *= -(4*vv*vv-(4*k+3)*(4*k+3))/(8*(2*k+2)*zz);
	}
	double omega = zz-(2*vv+1)*PI/4;
	return sqrt(2/PI*rz)*(cos(omega)*bff-sin(omega)*bgg);
}
STATIC double AngerJ_asymp( double vv, double zz )
{
	double ff = 0., gg=0., fk=1., gk=1., rz=1./zz, rzsq = rz*rz, 
		fp=0.;
	for (long k=0; k<100; ++k)
	{
		if (k > 0)
		{
			fk *= (vv*vv-(2*k-1)*(2*k-1))*rzsq;
			gk *= (vv*vv-4*k*k)*rzsq;
			if (fk > fp)
				break;
		}
		fp = fk;
		ff += fk;
		gg += gk;
	}
	return bessjnu(vv,zz)+sin(PI*vv)/(PI*zz)*(ff-(vv*rz)*gg);
}
STATIC double AngerJ( double vv, double zz )
{
	long int rep = 0, ddiv, divsor;

	double y = 0.0;

	DEBUG_ENTRY( "AngerJ()" );

	/* Estimate number of peaks in integrand.  */                            
	/* Divide region of integration by number  */
	/*  peaks in region.                       */
	if( (fabs(vv)) - (int)(fabs(vv)) > 0.5 )
		ddiv = (int)(fabs(vv)) + 1;
	else 
		ddiv = (int)(fabs(vv));

	divsor  = ((ddiv == 0) ? 1 : ddiv);
	vJint = vv;
	zJint = zz;

	for( rep = 0; rep < divsor; rep++ )
	{
		double
			rl = (((double) rep)/((double) divsor)),
			ru = (((double) (rep+1))/((double) divsor)),
			x_low = (PI * rl),
			x_up  = (PI * ru);      

		y += qg32( x_low, x_up, Jint );
	}

	return( y / PI );
}

/******************************************************************************/
/******************************************************************************/
/*                                                                            */
/*    Semi-Classical Quantum Defect Radial Integral                           */      
/*                                                                            */
/*   See for example                                                          */
/*     Atomic, Molecular & Optical Physics Handbook                           */
/*     Gordon W. F. Drake; Editor                                             */
/*     AIP Press                                                              */
/*     Woddbury, New York.                                                    */
/*     1996                                                                   */
/*                                                                            */
/* NOTE:: we do not include the Bohr Radius a_o in the                        */
/*           definition of of  R(n,L;n'L') as per Drake.                      */
/*                                                                            */
/*                                                                            */
/*                   1  (n_c)^2 | {      D_l max(L,L') }                      */
/*    R(n,L;n'L') = --- ------- | { 1 -  ------------- } J( D_n-1; -x )  -    */
/*                   Z   2 D_n  | {          n_c       }                      */
/*                                                                            */
/*                                                                            */
/*                    {      D_L max(L,L') }                                  */
/*                -   { 1 +  ------------- } J( D_n+1; -x )                   */
/*                    {          n_c       }                                  */
/*                                                                            */
/*                                                                            */
/*                     2                    |                                 */
/*                  + ---  sin(Pi D_n)(1-e) |                                 */
/*                     Pi                   |                                 */
/*                                                                            */
/*  where                                                                     */
/*        n_c = (2n*'n*)/(n*'+n*)                                             */
/*                                                                            */
/*       Here is the quantity Drake gives...                                  */
/*            n_c = ( 2.0 * nstar * npstar ) / ( nstar + npstar );            */
/*                                                                            */
/*       while V.A. Davidkin uses                                             */
/*            n_c = sqrt(  nstar * npstar  );                                 */
/*                                                                            */
/*        D_n = n*' - n*                                                      */
/*                                                                            */      
/*        D_L = L*' - L*                                                      */
/*                                                                            */
/*        x = e D_n                                                           */
/*                                                                            */
/*        Lmx  = max(L',L)                                                    */
/*                                                                            */
/*        e = sqrt( 1 - {Lmx/n_c}^2 )                                         */
/*                                                                            */
/*                                                                            */
/*       Here n* = n - qd where qd is the quantum defect                      */
/*                                                                            */
/******************************************************************************/
/******************************************************************************/
STATIC double scqdri(/* upper and lower quantum numbers...n's are effective	*/
			  double nstar, long int l,
			  double npstar, long int lp,
              double iz
              )
{
	double n_c = ((2.0 * nstar * npstar ) / ( nstar + npstar ));

	double D_n = (nstar - npstar);
	double D_l = (double) ( l - lp );
	double lg  = (double) ( (lp > l) ? lp : l);

	double h = (lg/n_c);
	double g = h*h;
	double f = ( 1.0 - g );
	double e = (( f >= 0.0) ? sqrt( f ) : 0.0 );

	double x  = (e * D_n);
	double z  = (-1.0 * x);
	double v1 = (D_n + 1.0);
	double v2 = (D_n - 1.0); 

	double d1,d2,d7,d8,d9,d34,d56,d6_1; 

	DEBUG_ENTRY( "scqdri()" );

	if( iz == 0.0 )
		iz += 1.0;

	if( D_n == 0.0 )
	{
		return( -1.0 );
	}

	if( D_n < 0.0 )
	{
		return( -1.0 );
	}

	if( f < 0.0 )
	{
		/* This can happen for certain  quantum defects   */
		/* in the lower n=1:l=0 state. In which case you  */
		/* probably should be using some other alogrithm  */
		/* or theory to calculate the dipole moment.      */
		return( -1.0 );
	}

	d1 = ( 1.0 / iz );

	d2 = (n_c * n_c)/(2.0 * D_n);

	if (0)
	{
		// Confirm accuracy of result of AngerJ function

		// relative accuracy is ~1e-10 in the worst case, where this error
		// may well be in summing series representation
		double vv = v2, zz = z, y = AngerJ(vv,zz);
		double yt = AngerJ_mac(vv,zz);
		double vp = (zz > 0.) ? vv : -vv;
		double ya = fabs(zz) < fabs(vv) ? 0. : AngerJ_asymp(vp,fabs(zz));
		printf("ANGER %15.8g %15.8g: %15.8g %15.8g (%15.8g) asymp %15.8g\n",
				 zz,vv,yt,y,fabs(yt-y)/(1e-37+fabs(yt)),ya);
	}

	d34 = (1.0 - ((D_l * lg)/n_c)) * AngerJ( v1, z );

	d56 = (1.0 + ((D_l * lg)/n_c)) * AngerJ( v2, z );

	d6_1 = PI * D_n;

	d7 = (2./PI) * sin( d6_1 ) * (1.0 - e);

	d8 = d1 * d2 * ( (d34) - (d56) + d7 );

	d9 = d8 * d8;

	ASSERT( D_n  > 0.0 );
	ASSERT( l  >= 0  );
	ASSERT( lp >= 0 );
	ASSERT( (l == lp + 1) || ( l == lp - 1) );
	ASSERT( n_c != 0.0 );
	ASSERT( f >= 0.0 );
	ASSERT( d9  > 0.0 );

	return( d9 );
}

STATIC double ForbiddenAuls( long ipHi, long ipLo, long nelem )
{
	double A;

	DEBUG_ENTRY( "ForbiddenAuls()" );

	int ipISO = ipHE_LIKE;

	if( (ipLo == ipHe1s1S) && (N_(ipHi) == 2) )
	{
		if( nelem == ipHELIUM )
		{
			// All of these but the second and third one (with values 1E-20) are from
			// >>refer	HeI	As	Lach, G, & Pachucki, K, 2001, Phys. Rev. A 64, 042510
			// 1E-20 is made up
			double ForbiddenHe[5] = { 1.272E-4, 1E-20, 1E-20, 177.58, 0.327 };

			A = ForbiddenHe[ipHi - 1];
			iso_put_error(ipHE_LIKE,nelem,ipHe2s3S ,ipHe1s1S,IPRAD, 0.01f, 0.20f);
			iso_put_error(ipHE_LIKE,nelem,ipHe2s1S ,ipHe1s1S,IPRAD, 0.01f, 0.05f);
			iso_put_error(ipHE_LIKE,nelem,ipHe2p3P0,ipHe1s1S,IPRAD, 0.00f, 0.00f);
			iso_put_error(ipHE_LIKE,nelem,ipHe2p3P1,ipHe1s1S,IPRAD, 0.01f, 0.05f);
			iso_put_error(ipHE_LIKE,nelem,ipHe2p3P2,ipHe1s1S,IPRAD, 0.01f, 0.01f);
		}
		else
		{
			switch ( (int)ipHi )
			{
			case 1: /* Parameters for 2^3S to ground transition.	*/
				/* >>refer	Helike	As	Lin, C.D., Johnson, W.R., and Dalgarno, A. 1977, 
				 * >>refercon	Phys. Rev. A 15, 1, 010015	*/
				// 2nu is treated elsewhere
				A = (3.9061E-7) * pow( (double)nelem+1., 10.419 );
				break;
			case 2: /* Parameters for 2^1S to ground transition.	*/
				A = iso_ctrl.SmallA; // 2nu is treated elsewhere
				break;
			case 3: /* Parameters for 2^3P0 to ground transition.	*/
				A = iso_ctrl.SmallA;
				break;
			case 4: /* Parameters for 2^3P1 to ground transition.	*/
				/* >>chng 06 jul 25, only use the fit to Johnson et al. values up to and
				 * including argon, where there calculation stops.  For higher Z use below */
				if( nelem <= ipARGON )
				{
					A = ( 11.431 * pow((double)nelem, 9.091) );
				}
				else
				{
					/* a fit to the Lin et al. 1977. values, which go past zinc. */
					A = ( 383.42 * pow((double)nelem, 7.8901) );
				}
				break;
			case 5: /* Parameters for 2^3P2 to ground transition.	*/
				/* fit to Lin et al. 1977 values.  This fit is good
				 * to 7% for the range from carbon to iron. The Lin et al. values
				 * differ from the Hata and Grant (1984) values (only given from 
				 * calcium to copper) by less than 2%. */
				A = ( 0.11012 * pow((double)nelem, 7.6954) ); 
				break;
			default:
				TotalInsanity();
			}
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
		}

		return A;
	}

	/* The next two cases are fits to probabilities given in 
	 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
	 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
	/* originally astro.ph. 0201454 */
	/* The involve Triplet P and Singlet S.  Rates for Triplet S to Singlet P 
	 * do not seem to be available.	*/

	/* Triplet P to Singlet S...Delta n not equal to zero!	*/
	else if( nelem>ipHELIUM && L_(ipHi)==1 && S_(ipHi)==3 && 
		L_(ipLo)==0 && S_(ipLo)==1 && N_(ipLo) < N_(ipHi) )
	{
		A = 8.0E-3 * exp(9.283/sqrt((double)N_(ipLo))) * pow((double)nelem,9.091) /
			pow((double)N_(ipHi),2.877);
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	/* Singlet S to Triplet P...Delta n not equal to zero!	*/
	else if( nelem > ipHELIUM && L_(ipHi)==0 && S_(ipHi)==1 && 
		L_(ipLo)==1 && S_(ipLo)==3 && N_(ipLo) < N_(ipHi) )
	{
		A = 2.416 * exp(-0.256*N_(ipLo)) * pow((double)nelem,9.159) / pow((double)N_(ipHi),3.111);

		if( ( (ipLo == ipHe2p3P0) || (ipLo == ipHe2p3P2) ) )
		{
			/* This is divided by 3 instead of 9, because value calculated is specifically for 2^3P1.
			 * Here we assume statistical population of the other two.	*/
			A *= (2.*(ipLo-3)+1.0)/3.0;
		}
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else if( ( ipLo == ipHe2s3S ) && ( ipHi == ipHe3s3S ) )
	{
		double As_3TripS_to_2TripS[9] = {
			7.86E-9, 4.59E-6, 1.90E-4, 2.76E-3, 2.25E-2,
			1.27E-1, 5.56E-1, 2.01E0, 6.26E0 };

		/* These M1 transitions given by 
		 * >>refer He-like As Savukov, I.M., Labzowsky, and Johnson, W.R. 2005, PhysRevA, 72, 012504 */
		if( nelem <= ipNEON )
		{
			A = As_3TripS_to_2TripS[nelem-1];
			/* 20% error is based on difference between Savukov, Labzowsky, and Johnson (2005)
			 * and Lach and Pachucki (2001) for the helium transition. */
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.2f,0.2f);
		}
		else
		{
			/* This is an extrapolation to the data given above.  The fit reproduces
			 * the above values to 10% or better. */
			A= 7.22E-9*pow((double)nelem, 9.33);
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.3f,0.3f);
		}
	}

	else if( ( ipLo == ipHe2s3S ) && ( ipHi == ipHe2p1P ) )
	{
		/* This transition,1.549 , given by Lach and Pachucki, 2001 for the atom */
		if( nelem == ipHELIUM )
		{
			A = 1.549;
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
		}
		else
		{
			/* This is a fit to data given in
			 * >>refer	He-like	As	Savukov, I.M., Johnson, W.R., & Safronova, U.I. 
			 * >>refercon	astro-ph 0205163	*/
			A= 0.1834*pow((double)nelem, 6.5735);
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
		}
	}

	else if( nelem==ipHELIUM && ipHi==4 && ipLo==3 )
	{
		/* >>refer	He	As	Bulatov, N.N. 1976, Soviet Astronomy, 20, 315 */
		fixit("find a transition probability for this 2^3P0 - 2^3P1 transition.");
		/* This is the 29.6 GHz line that can be seen in radio galaxies. */
		/** \todo	2	find a transition probability for this 2^3P0 - 2^3P1 transition.
		 * It will require a bit of trickery to insert into the rate matrix, 
		 * because of the fact that the lower level has a higher index.  
		 * See discussion "Energy order within 2 3P" near the top of helike.c */
		A = 5.78E-12;
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);

	}

	else if( nelem==ipHELIUM && ipHi==5 && ipLo==4 )
	{
		fixit("find a transition probability for this 2^3P1 - 2^3P2 transition.");
		/* This is the 3 GHz line that can be seen in radio galaxies. */
		/** \todo	2	find a transition probability for this 2^3P1 - 2^3P2 transition.
		 * It will require a bit of trickery to insert into the rate matrix, 
		 * because of the fact that the lower level has a higher index.  
		 * See discussion "Energy order within 2 3P" near the top of helike.c */
		A = 3.61E-15;
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem==ipHELIUM && ipHi==ipHe3p3P && ipLo==ipHe1s1S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 44.326 * (1./3.) + 0.1146547 * (5./9.);
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem==ipHELIUM && ipHi==ipHe3p3P && ipLo==ipHe2s1S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 1.459495 * (1./3.) + 3.6558e-5 * (5./9.);
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem==ipHELIUM && ipHi==ipHe3p3P && ipLo==ipHe3s1S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 2.2297e-3 * (1./3.);
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem==ipHELIUM && ipHi==ipHe3p1P && ipLo==ipHe2s3S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 0.1815436;
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem==ipHELIUM && ipHi==ipHe3p1P && ipLo==ipHe3s3S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 0.14895852;
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else if( ( ipLo==0 && L_(ipHi)==2 && S_(ipHi)==1 ) ||
		( ipLo==ipHe2s3S && L_(ipHi)==2 && S_(ipHi)==3 ) ||
		( ipLo==ipHe2s1S && L_(ipHi)==2 && S_(ipHi)==1 ) ||
		( N_(ipLo)==2 && L_(ipLo)==1 && S_(ipLo)==3 && N_(ipHi)>=3 && L_(ipHi)==1 && S_(ipHi)==3 ) ||
		( ipLo==ipHe2p1P && L_(ipHi)==1 && S_(ipHi)==1 ) )
	{
		// >>refer       Helike  QOS     Cann, N. M. \& Thakkar, A. J. 2002, J.Phys.B, 35, 421

		static const double f_params[5][4][3] = {
			{
		/* 1,0,3,2,1 */ {9.360591E-007, -3.1937E-006,  3.5186E-006},
		/* 1,0,4,2,1 */ {4.631435E-007, -1.4973E-006,  1.4848E-006},
		/* 1,0,5,2,1 */ {2.493912E-007, -7.8658E-007,  7.3994E-007},
		/* 1,0,6,2,1 */ {1.476742E-007, -4.5953E-007,  4.1932E-007}},
			{
		/* 2,0,3,2,3 */ {1.646733E-006, -2.0028E-006, -1.3552E-006},
		/* 2,0,4,2,3 */ {9.120593E-008,  3.1301E-007, -3.2190E-007},
		/* 2,0,5,2,3 */ {1.360965E-008,  1.1467E-007,  8.6977E-008},
		/* 2,0,6,2,3 */ {3.199421E-009,  4.5485E-008,  1.1016E-007}},
			{
		/* 2,0,3,2,1 */ {1.646733E-006, -2.9720E-006,  1.5367E-006},
		/* 2,0,4,2,1 */ {9.120593E-008, -3.9037E-008,  3.9156E-008},
		/* 2,0,5,2,1 */ {1.360965E-008,  1.4671E-008,  1.5626E-008},
		/* 2,0,6,2,1 */ {3.199421E-009,  8.9806E-009,  1.2436E-008}},
			{
		/* 2,1,3,1,3 */ {1.543812E-007, -2.8220E-007,  3.0261E-008},
		/* 2,1,4,1,3 */ {3.648237E-008, -6.6824E-008,  4.5879E-009},
		/* 2,1,5,1,3 */ {1.488556E-008, -2.7304E-008,  1.6628E-009},
		/* 2,1,6,1,3 */ {7.678610E-009, -1.4112E-008,  6.8453E-010}},
			{
		/* 2,1,3,1,1 */ {1.543812E-007, -2.8546E-007,  4.6605E-008},
		/* 2,1,4,1,1 */ {3.648237E-008, -6.8422E-008,  1.7038E-008},
		/* 2,1,5,1,1 */ {1.488556E-008, -2.8170E-008,  8.5012E-009},
		/* 2,1,6,1,1 */ {7.678610E-009, -1.4578E-008,  4.6686E-009}}
				};

		ASSERT( ipLo <= 6 );
		long index_lo;
		if( ipLo==ipHe2p1P )
			index_lo = 4;
		else if( ipLo==ipHe2p3P1 || ipLo==ipHe2p3P2 )
			index_lo = 3;
		else
			index_lo = ipLo;

		ASSERT( N_(ipHi) >= 3 );
		long index_hi = MIN2( N_(ipHi), 6 ) - 3;
		double f_lu = POW2(nelem+1) * ( 
			f_params[index_lo][index_hi][0] + 
			f_params[index_lo][index_hi][1]/(nelem+1) + 
			f_params[index_lo][index_hi][2]/POW2(nelem+1) );

		// extrapolate for higher upper n
		if( N_(ipHi) > 6 )
			f_lu *= pow3( 6. / N_(ipHi) );

		double gLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].g();
		double gHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].g();
		double eWN = iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo).EnergyWN();

		A = eina( gLo * f_lu, eWN, gHi );
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}
	
	else
	{
		/* Current transition is not supported.	*/
		A = iso_ctrl.SmallA;
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	ASSERT( A > 0.);

	return A;   
}

/* Calculates Einstein A for a given transition.	*/
double he_1trans( 
		/* charge on c scale, Energy is wavenumbers, Einstein A	*/
		long nelem , double Enerwn ,
		/* quantum numbers of upper level:	*/
		double Eff_nupper, long lHi, long sHi, long jHi,
		/* and of lower level: */
		double Eff_nlower, long lLo, long sLo, long jLo,
		realnum *error1, realnum *error2 )
		/* Note j is only necessary for 2 triplet P...for all other n,l,s,
		 * j is completely ignored.	*/
{
	double RI2, Aul;
	long nHi, nLo, ipHi, ipLo;

	DEBUG_ENTRY( "he_1trans()" );

	ASSERT(nelem > ipHYDROGEN);

	/* Since 0.4 is bigger than any defect, adding that to the effective principle quantum number,
	 * and truncating to an integer will produce the principal quantum number.	*/
	nHi = (int)(Eff_nupper + 0.4);
	nLo = (int)(Eff_nlower + 0.4);

	/* Make sure this worked correctly.	*/
	ASSERT( fabs(Eff_nupper-(double)nHi) < 0.4 );
	ASSERT( fabs(Eff_nlower-(double)nLo) < 0.4 );

	ipHi = iso_sp[ipHE_LIKE][nelem].QuantumNumbers2Index[nHi][lHi][sHi];
	if( (nHi==2) && (lHi==1) && (sHi==3) )
	{
		ASSERT( (jHi>=0) && (jHi<=2) );
		ipHi -= (2 - jHi);
	}

	ipLo = iso_sp[ipHE_LIKE][nelem].QuantumNumbers2Index[nLo][lLo][sLo];
	if( (nLo==2) && (lLo==1) && (sLo==3) )
	{
		ASSERT( (jLo>=0) && (jLo<=2) );
		ipLo -= (2 - jLo);
	}

	ASSERT( iso_sp[ipHE_LIKE][nelem].st[ipHi].n() == nHi );
	if( nHi <= iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
	{
		ASSERT( iso_sp[ipHE_LIKE][nelem].st[ipHi].l() == lHi );
		ASSERT( iso_sp[ipHE_LIKE][nelem].st[ipHi].S() == sHi );
	}
	ASSERT( iso_sp[ipHE_LIKE][nelem].st[ipLo].n() == nLo );
	if( nLo <= iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
	{
		ASSERT( iso_sp[ipHE_LIKE][nelem].st[ipLo].l() == lLo );
		ASSERT( iso_sp[ipHE_LIKE][nelem].st[ipLo].S() == sLo );
	}

	/* First do allowed transitions	*/
	if( (sHi == sLo) && (abs((int)(lHi - lLo)) == 1) )
	{
		Aul = -2.;

		/* For clarity, let's do this in two separate chunks...one for helium, one for everything else.	*/
		if( nelem == ipHELIUM )
		{
			/* Retrieve transition probabilities for Helium.	*/
			/* >>refer He	As	Drake, G.W.F., Atomic, Molecular, and Optical Physics Handbook */
			if( ipHi <= MAX_TP_INDEX && nHi <= iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max )
			{
				/*Must not be accessed by collapsed levels!	*/
				ASSERT( ipHi < ( iso_sp[ipHE_LIKE][ipHELIUM].numLevels_max - iso_sp[ipHE_LIKE][ipHELIUM].nCollapsed_max ) );
				ASSERT( ipLo < ( iso_sp[ipHE_LIKE][ipHELIUM].numLevels_max - iso_sp[ipHE_LIKE][ipHELIUM].nCollapsed_max ) );
				ASSERT( ipHi > 2 );

				Aul = TransProbs[nelem][ipHi][ipLo];

				*error1 = 0.0001f;
				*error2 = 0.002f;
			}
			else
			{
				ASSERT( nHi > 2 );
				long ipHiTemp = iso_sp[ipHE_LIKE][nelem].IndexIfAllResolved[nHi][lHi][sHi];
				long ipLoTemp = iso_sp[ipHE_LIKE][nelem].IndexIfAllResolved[nLo][lLo][sLo];
				if( (nLo==2) && (lLo==1) && (sLo==3) )
				{
					ASSERT( (jLo>=0) && (jLo<=2) );
					ipLoTemp -= (2 - jLo);
				}
				ASSERT( ipLoTemp < ipHiTemp );

				if( ipHiTemp <= MAX_TP_INDEX )
				{
					Aul = TransProbs[nelem][ipHiTemp][ipLoTemp];
					*error1 = 0.0001f;
					*error2 = 0.002f;
				}
			}

			if( Aul < 0. )
			{
				/* Here are the Lyman transitions.	*/
				if( nLo==1 )
				{
					ASSERT( (lHi == 1) && (sHi == 1) );

					/* these fits calculated from Drake A's (1996) */
					Aul = (1.59208e10) / pow3(Eff_nupper);
					ASSERT( Aul > 0.);
					*error1 = 0.0002f;
					*error2 = 0.002f;
				}

				/* last resort for transitions involving significant defects, 
				 * except that highest lLo are excluded */
				else if( lHi>=2 && lLo>=2 && nHi>nLo )
				{
					Aul = H_Einstein_A(nHi ,lHi , nLo , lLo , nelem);
					ASSERT( Aul > 0.);
					*error1 = 0.006f;
					*error2 = 0.04f;
				}
				/* These fits come from extrapolations of Drake's oscillator strengths
				 * out to the series limit.  We also use this method to obtain threshold
				 * photoionization cross-sections for the lower level of each transition here. */
				/* See these two references :
				 * >>refer	He	As	Hummer, D. G. \& Storey, P. J. 1998, MNRAS 297, 1073 
				 * >>refer	Seaton's Theorem	Seaton, M. J. 1958, MNRAS 118, 504 */
				else if( nHi>10 && nLo<=5 && lHi<=2 && lLo<=2 )
				{
					int paramSet=0;
					double emisOscStr, x, a, b, c;
					static const double extrapol_Params[2][4][4][3] = {
						/* these are for singlets */
						{	
							{	/* these are P to S */
								{	0.8267396	,	1.4837624	,	-0.4615955	},
								{	1.2738405	,	1.5841806	,	-0.3022984	},
								{	1.6128996	,	1.6842538	,	-0.2393057	},
								{	1.8855491	,	1.7709125	,	-0.2115213	}},
							{	/* these are S to P */
								{	-1.4293664	,	2.3294080	,	-0.0890470	},
								{	-0.3608082	,	2.3337636	,	-0.0712380	},
								{	0.3027974	,	2.3326252	,	-0.0579008	},
								{	0.7841193	,	2.3320138	,	-0.0497094	}},
							{	/* these are D to P */
								{	1.1341403	,	3.1702435	,	-0.2085843	},
								{	1.7915926	,	2.4942946	,	-0.2266493	},
								{	2.1979400	,	2.2785377	,	-0.1518743	},
								{	2.5018229	,	2.1925720	,	-0.1081966	}},
							{	/* these are P to D */
								{	0.0000000	,	0.0000000	,	0.0000000	},
								{	-2.6737396	,	2.9379143	,	-0.3805367	},
								{	-1.4380124	,	2.7756396	,	-0.2754625	},
								{	-0.6630196	,	2.6887253	,	-0.2216493	}},
						},
						/* these are for triplets */
						{	
							{	/* these are P to S */
								{	0.3075287	,	0.9087130	,	-1.0387207	},
								{	0.687069	,	1.1485864	,	-0.6627317	},
								{	0.9776064	,	1.3382024	,	-0.5331906	},
								{	1.2107725	,	1.4943721	,	-0.4779232	}},
							{	/* these are S to P */
								{	-1.3659605	,	2.3262253	,	-0.0306439	},
								{	-0.2899490	,	2.3279391	,	-0.0298695	},
								{	0.3678878	,	2.3266603	,	-0.0240021	},
								{	0.8427457	,	2.3249540	,	-0.0194091	}},
							{	/* these are D to P */
								{	1.3108281	,	2.8446367	,	-0.1649923	},
								{	1.8437692	,	2.2399326	,	-0.2583398	},
								{	2.1820792	,	2.0693762	,	-0.1864091	},
								{	2.4414052	,	2.0168255	,	-0.1426083	}},
							{	/* these are P to D */
								{	0.0000000	,	0.0000000	,	0.0000000	},
								{	-1.9219877	,	2.7689624	,	-0.2536072	},
								{	-0.7818065	,	2.6595150	,	-0.1895313	},
								{	-0.0665624	,	2.5955623	,	-0.1522616	}},
						}
					};

					if( lLo==0 )
					{
						paramSet = 0;
					}
					else if( lLo==1 && lHi==0 )
					{
						paramSet = 1;
					}
					else if( lLo==1 && lHi==2 )
					{
						paramSet = 2;
					}
					else if( lLo==2 )
					{
						paramSet = 3;
						ASSERT( lHi==1 );
					}

					ASSERT( (int)((sHi-1)/2) == 0 || (int)((sHi-1)/2) == 1 );
					a = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][0];
					b = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][1];
					c = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][2];
					ASSERT( Enerwn > 0. );
					x = log( iso_sp[ipHE_LIKE][nelem].fb[ipLo].xIsoLevNIonRyd*RYD_INF/Enerwn );

					emisOscStr = exp(a+b*x+c*x*x)/pow3(Eff_nupper)*
						(2.*lLo+1)/(2.*lHi+1);

					Aul = TRANS_PROB_CONST*Enerwn*Enerwn*emisOscStr;

					if( nLo==2 && lLo==1 && sLo==3 )
					{
						ASSERT( jLo == ipLo-3 );
						Aul *= (2.*jLo+1.0)/9.0;
					}

					ASSERT( Aul > 0. );
					*error1 = 0.0002f;
					*error2 = 0.002f;
				}
				else
				{
					/* Calculate the radial integral from the quantum defects.	*/
					RI2 = scqdri(Eff_nupper,lHi,Eff_nlower,lLo,(double)(ipHELIUM));
					ASSERT( RI2 > 0. );
					/* Convert radial integral to Aul.	*/
					Aul = ritoa(lHi,lLo,ipHELIUM,Enerwn,RI2);
					/* radial integral routine does not recognize fine structure.
					 * Here we split 2^3P.	*/
					if( nLo==2 && lLo==1 && sLo==3 )
					{
						ASSERT( jLo == ipLo-3 );
						Aul *= (2.*jLo+1.0)/9.0;
					}

					ASSERT( Aul > 0. );
					*error1 = 0.01f;
					*error2 = 0.07f;
				}
			}
		}

		/* Heavier species	*/
		else
		{
			/* Retrieve transition probabilities for Helike ions.	*/
			/* >>refer He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
			 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J, originally astro.ph. 0201454 */
			if( ipHi <= MAX_TP_INDEX && nHi <= iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
			{
				/*Must not be accessed by collapsed levels!	*/
				ASSERT( ipHi < ( iso_sp[ipHE_LIKE][nelem].numLevels_max - iso_sp[ipHE_LIKE][nelem].nCollapsed_max ) );
				ASSERT( ipLo < ( iso_sp[ipHE_LIKE][nelem].numLevels_max - iso_sp[ipHE_LIKE][nelem].nCollapsed_max ) );
				ASSERT( ipHi > 2 );

				Aul = TransProbs[nelem][ipHi][ipLo];
				*error1 = 0.1f;
				*error2 = 0.1f;
			}
			else
			{
				ASSERT( nHi > 2 );
				long ipHiTemp = iso_sp[ipHE_LIKE][nelem].IndexIfAllResolved[nHi][lHi][sHi];
				long ipLoTemp = iso_sp[ipHE_LIKE][nelem].IndexIfAllResolved[nLo][lLo][sLo];
				if( (nLo==2) && (lLo==1) && (sLo==3) )
				{
					ASSERT( (jLo>=0) && (jLo<=2) );
					ipLoTemp -= (2 - jLo);
				}
				ASSERT( ipLoTemp < ipHiTemp );

				if( ipHiTemp <= MAX_TP_INDEX )
				{
					Aul = TransProbs[nelem][ipHiTemp][ipLoTemp];
					*error1 = 0.1f;
					*error2 = 0.1f;
				}
			}

			if( Aul < 0. )
			{
				/* Do same-n transitions. */
				if( nLo==nHi && lHi<=2 && lLo<=2 )
				{
					/* These are 2p3Pj to 2s3S fits to (low Z) Porquet & Dubau (2000) & 
					 * (high Z) NIST Atomic Spectra Database.	*/
					if( nLo==2 && lLo==0 && sLo==3 )
					{
						ASSERT( nHi==2 && lHi==1 && sHi==3 );
						ASSERT( jHi==ipHi-3 );

						if( jHi==0 )
							Aul = 3.31E7 + 1.13E6 * pow((double)nelem+1.,1.76);
						else if( jHi==1 )
							Aul = 2.73E7 + 1.31E6 * pow((double)nelem+1.,1.76);
						else if( jHi==2 )
							Aul = 3.68E7 + 1.04E7 * exp(((double)nelem+1.)/5.29);
						else
						{
							fprintf(ioQQQ," PROBLEM Impossible value for jHi in he_1trans.\n");
							TotalInsanity();
						}
					}

					/* These are 2p1P to 2s1S fits to data from TOPbase.	*/
					else if( ( nLo==2 && lLo==0 && sLo==1 ) && ( nHi==2 && lHi==1 && sHi==1 ) )
					{
						Aul = 5.53E6 * exp( 0.171*(nelem+1.) );
					}

					else 
					{
						/* This case should only be entered if n > 2.  Those cases were done above.	*/
						ASSERT( nLo > 2 );

						/* Triplet P to triplet S, delta n = 0	*/
						if( (lHi == 1) && (sHi == 3) && (lLo == 0) && (sLo == 3))
						{
							Aul = 0.4 * 3.85E8 * pow((double)nelem,1.6)/pow((double)nHi,5.328);
						}
						/* Singlet P to singlet D, delta n = 0	*/
						else if( (lHi == 1) && (sHi == 1) && (lLo == 2) && (sLo == 1))
						{
							Aul = 1.95E4 * pow((double)nelem,1.6) / pow((double)nHi, 4.269);
						}
						/* Singlet P to singlet S, delta n = 0	*/
						else if( (lHi == 1) && (sHi == 1) && (lLo == 0) )
						{
							Aul = 6.646E7 * powpq((double)nelem,3,2) / pow((double)nHi, 5.077);
						}
						else 
						{
							ASSERT( (lHi == 2) && (sHi == 3)  && (lLo == 1) );
							Aul = 3.9E6 * pow((double)nelem,1.6) / pow((double)nHi, 4.9);
							if( (lHi >2) || (lLo > 2) )
								Aul *= (lHi/2.);
							if(lLo > 2)
								Aul *= (1./9.);
						}
					}
					ASSERT( Aul > 0.);
				}

				/* assume transitions involving F and higher orbitals are hydrogenic.	*/
				else if( (nHi > nLo) && ((lHi > 2) || (lLo > 2)) )
				{
					Aul = H_Einstein_A(nHi ,lHi , nLo , lLo , nelem);
					ASSERT( Aul > 0.);
				}

				/* These transitions are of great importance, but the below radial integral 
				 * routine fails to achieve desirable accuracy, so these are fits as produced 
				 * from He A's for nupper through 9.  They are transitions to ground and 
				 * 2, 3, and 4 triplet S.	*/
				else if( ( nLo==1 ) || ( nLo==2 && lLo==0 ) ||
					( nLo==3 && lLo==0 && sLo==3 ) ||
					( nLo==4 && lLo==0 && sLo==3 ) )
				{
					/* Here are the Lyman transitions.	*/
					if( nLo==1 )
					{
						ASSERT( (lHi == 1) && (sHi == 1) );

						/* In theory, this Z dependence should be Z^4, but values from TOPbase 
						 * suggest 3.9 is a more accurate exponent.	Values from 
						 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
						 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
						/* originally astro.ph. 0201454  */
						Aul = 1.375E10 * pow((double)nelem, 3.9) / pow((double)nHi,3.1);
					}

					/* Here are the Balmer transitions.	*/
					else if( nLo==2 && lLo==0 && sLo==1 )
					{
						ASSERT( (lHi == 1) && (sHi == 1) );

						/* More fits, as above. */
						Aul = 5.0e8 * pow4((double)nelem) / pow((double)nHi, 2.889);
					}

					/* Here are transitions down to triplet S	*/
					else
					{
						ASSERT( (lLo==0) && (sLo==3) );
						ASSERT( (lHi==1) && (sHi==3) );

						/* These fits also as above. */
						if( nLo == 2 )
							Aul = 1.5 * 3.405E8 * pow4((double)nelem) / pow((double)nHi, 2.883);
						else if( nLo == 3 )
							Aul = 2.5 * 4.613E7 * pow4((double)nelem) / pow((double)nHi, 2.672);
						else 
							Aul = 3.0 * 1.436E7 * pow4((double)nelem) / pow((double)nHi, 2.617);
					}

					ASSERT( Aul > 0.);
				}

				/* Every other allowed transition is calculated as follows.	*/
				else
				{
					/* Calculate the radial integral from the quantum defects.	*/
					RI2 = scqdri(Eff_nupper,lHi,Eff_nlower,lLo,(double)(nelem));
					/* Convert radial integral to Aul.	*/
					Aul = ritoa(lHi,lLo,nelem,Enerwn,RI2);
					/* radial integral routine does not recognize fine structure.
					 * Here we split 2^3P.	*/
					if( nLo==2 && lLo==1 && sLo==3 && Aul > iso_ctrl.SmallA )
					{
						ASSERT( jLo == ipLo-3 );
						Aul *= (2.*jLo+1.0)/9.0;
					}

				}
				*error1 = 0.1f;
				*error2 = 0.1f;
			}
		}
	}

	/* Now do forbidden transitions from 2-1 ... */
	/* and those given by  
	 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
	 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
	/* originally astro.ph. 0201454  
	 * for heavy elements. These are triplet P to singlet S, 
	 * going either up or down...Triplet S to Singlet P are not included, as they are far weaker.	*/
	else 
	{
		ASSERT( (sHi != sLo) || (abs((int)(lHi - lLo)) != 1) );
		Aul = ForbiddenAuls(ipHi, ipLo, nelem);
		ASSERT( Aul > 0. );
	}

	Aul = MAX2( Aul, iso_ctrl.SmallA );
	ASSERT( Aul >= iso_ctrl.SmallA );

	/* negative energy for a transition with substantial transition probability
	 * would be major logical error - but is ok for same-n l transitions */
	if( Enerwn < 0. && Aul > iso_ctrl.SmallA )
	{
		fprintf( ioQQQ," he_1trans hit negative energy, nelem=%li, val was %f \n", nelem ,Enerwn );
	}

	return Aul;
}

STATIC double helike_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo, double Enerwn )
{
	DEBUG_ENTRY( "helike_transprob_collapsed_to_collapsed()" );

	double Aul = 0.;
	ASSERT( Enerwn > 0. );
	ASSERT( nHi > nLo );

	for( long lLo=0; lLo < nLo; ++lLo )
	{
		double Aul_sing = helike_transprob_collapsed_to_resolved( nelem, nHi, nLo, lLo, 1, -1, Enerwn );
		double Aul_trip = helike_transprob_collapsed_to_resolved( nelem, nHi, nLo, lLo, 3, -1, Enerwn );
		Aul += Aul_sing + Aul_trip;
	}

	return Aul;
}

STATIC double helike_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo, long sLo, long jLo, double Enerwn )
{
	DEBUG_ENTRY( "helike_transprob_collapsed_to_resolved()" );

	double n_eff_hi = nHi - helike_quantum_defect(nelem,nHi,-1,-1,-1);
	double n_eff_lo = nLo - helike_quantum_defect(nelem,nLo,lLo,sLo,jLo);
	realnum error1, error2;

	long ipISO = ipHE_LIKE;
	t_iso_sp *sp = &iso_sp[ipISO][nelem];
	long nResMax = sp->n_HighestResolved_max;
	long ipLoRes = sp->IndexIfAllResolved[nLo][lLo][sLo];
	// make sure 2^3Pj terms are fully specified
	if( nLo==2 && lLo==1 && sLo==3 )
	{
		ASSERT( jLo>=0 && jLo<=2 );
		ipLoRes -= (2 - jLo);
	}
	ASSERT( nLo > nResMax || ipLoRes <= sp->numLevels_max - sp->nCollapsed_max );
	
	/* Lower level resolved, upper not. First calculate Aul
	 * from upper level with ang mom one higher.	*/
	double Aul = he_1trans( nelem, Enerwn, 
			n_eff_hi, lLo+1, sLo, -1, n_eff_lo, lLo, sLo, jLo,
			&error1, &error2);

	sp->CachedAs[ nHi-nResMax-1 ][ ipLoRes ][0] = (realnum)Aul;

	Aul *= (2.*(lLo+1.)+1.) * sLo / (4.*(double)nHi*(double)nHi);

	if( lLo != 0 )
	{
		/* For all l>0, add in transitions from upper level with ang mom one lower.	*/
		double Aul1 = he_1trans( nelem, Enerwn, 
				n_eff_hi, lLo-1, sLo, -1, n_eff_lo, lLo, sLo, jLo,
				&error1, &error2);

		sp->CachedAs[ nHi-nResMax-1 ][ ipLoRes ][1] = (realnum)Aul1;

		/* now add in this part after multiplying by stat weight for lHi = lLo-1.	*/
		Aul += Aul1*(2.*(lLo-1.)+1.) * sLo / (4.*(double)nHi*(double)nHi);
	}
	else
		sp->CachedAs[ nHi-nResMax-1 ][ ipLoRes ][1] = 0.f;

	ASSERT( Aul > 0.);

	return Aul;
}

void DoFSMixing( long nelem, long ipLoSing, long ipHiSing )
{
	/* Every bit of this routine is based upon the singlet-triplet mixing formalism given in
	 * >>refer He	FSM	Drake, G. W. F. 1996, in Atomic, Molecular, \& Optical Physics Handbook,
	 * >>refercon	ed. G. W. F. Drake (New York: AIP Press).	
	 * That formalism mixes the levels themselves, but since this code is not J-resolved, we simulate
	 * that by mixing only the transition probabilities.  We find results comparable to those calculated
	 * in the fully J-resolved model spearheaded by Rob Bauman, and described in
	 * >>refer	He	FSM	Bauman, R. P., Porter, R. L., Ferland, G. J., \& MacAdam, K. B. 2005, ApJ, accepted */
	long int nHi, lHi, sHi, nLo, lLo, sLo, ipHiTrip, ipLoTrip;
	double Ass, Att, Ast, Ats;
	double SinHi, SinLo, CosHi, CosLo;
	double HiMixingAngle, LoMixingAngle , error;
	double Kss, Ktt, Kts, Kst, fss, ftt, fssNew, fttNew, ftsNew, fstNew, temp;

	DEBUG_ENTRY( "DoFSMixing()" );

	nHi = iso_sp[ipHE_LIKE][nelem].st[ipHiSing].n();
	lHi = iso_sp[ipHE_LIKE][nelem].st[ipHiSing].l();
	sHi = iso_sp[ipHE_LIKE][nelem].st[ipHiSing].S();
	nLo = iso_sp[ipHE_LIKE][nelem].st[ipLoSing].n();
	lLo = iso_sp[ipHE_LIKE][nelem].st[ipLoSing].l();
	sLo = iso_sp[ipHE_LIKE][nelem].st[ipLoSing].S();

	if( ( sHi == 3 || sLo == 3 ) ||
	    ( abs(lHi - lLo) != 1 ) ||
	    ( nLo < 2 ) ||
	    ( lHi <= 1 || lLo <= 1 ) ||
	    ( nHi == nLo && lHi == 1 && lLo == 2 ) ||
	    ( nHi > nLo && lHi != 1 && lLo != 1 ) )
	{
		return;
	}

	ASSERT( lHi > 0 );
	/*ASSERT( (lHi > 1) && (lLo > 1) );*/

	ipHiTrip = iso_sp[ipHE_LIKE][nelem].QuantumNumbers2Index[nHi][lHi][3];
	ipLoTrip = iso_sp[ipHE_LIKE][nelem].QuantumNumbers2Index[nLo][lLo][3];

	if( lHi == 2 )
	{
		HiMixingAngle = 0.01;
	}
	else if( lHi == 3 )
	{
		HiMixingAngle = 0.5;
	}
	else
	{
		HiMixingAngle = PI/4.;
	}

	if( lLo == 2 )
	{
		LoMixingAngle = 0.01;
	}
	else if( lLo == 3 )
	{
		LoMixingAngle = 0.5;
	}
	else
	{
		LoMixingAngle = PI/4.;
	}

	/* These would not work correctly if l<=1 were included in this treatment!	*/
	ASSERT( ipHiTrip > ipLoTrip );
	ASSERT( ipHiTrip > ipLoSing );
	ASSERT( ipHiSing > ipLoTrip );
	ASSERT( ipHiSing > ipLoSing );

	SinHi = sin( HiMixingAngle );
	SinLo = sin( LoMixingAngle );
	CosHi = cos( HiMixingAngle );
	CosLo = cos( LoMixingAngle );

	Kss = iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).EnergyWN();
	Ktt = iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).EnergyWN();
	Kst = iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoTrip).EnergyWN();
	Kts = iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoSing).EnergyWN();

	fss = iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).Emis().Aul()/TRANS_PROB_CONST/POW2(Kss)/
		(*iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).Lo()).g()*(*iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).Hi()).g();

	ftt = iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).Emis().Aul()/TRANS_PROB_CONST/POW2(Ktt)/
		(*iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).Lo()).g()*(*iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).Hi()).g();

	temp = sqrt(fss/Kss)*CosHi*CosLo + sqrt(ftt/Ktt)*SinHi*SinLo;
	fssNew = Kss*POW2( temp );
	temp = sqrt(fss/Kss)*SinHi*SinLo + sqrt(ftt/Ktt)*CosHi*CosLo;
	fttNew = Ktt*POW2( temp );
	temp = sqrt(fss/Kss)*CosHi*SinLo - sqrt(ftt/Ktt)*SinHi*CosLo;
	fstNew = Kst*POW2( temp );
	temp = sqrt(fss/Kss)*SinHi*CosLo - sqrt(ftt/Ktt)*CosHi*SinLo;
	ftsNew = Kts*POW2( temp );

	Ass = TRANS_PROB_CONST*POW2(Kss)*fssNew*(*iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).Lo()).g()/
		(*iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).Hi()).g();

	Att = TRANS_PROB_CONST*POW2(Ktt)*fttNew*(*iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).Lo()).g()/
		(*iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).Hi()).g();

	Ast = TRANS_PROB_CONST*POW2(Kst)*fstNew*(*iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoTrip).Lo()).g()/
		(*iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoTrip).Hi()).g();

	Ats = TRANS_PROB_CONST*POW2(Kts)*ftsNew*(*iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoSing).Lo()).g()/
		(*iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoSing).Hi()).g();

	error = fabs( ( iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).Emis().Aul()+
		iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).Emis().Aul() )/
		(Ass+Ast+Ats+Att) - 1.f );

	if( error > 0.001 )
	{
		fprintf( ioQQQ, "FSM error: %e\t Element: %ld\t States: LS %li HS %li LT %li HT %li\t New values: Ass %e Att %e Ast %e Ats %e\n", error,
			nelem, ipLoSing, ipHiSing, ipLoTrip, ipHiTrip, Ass, Att, Ast, Ats );
	}

	if( !iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).hasEmis() )
	{
		iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).AddLine2Stack();
		iso_setRedisFun(ipHE_LIKE, nelem, ipLoSing, ipHiSing);
		iso_setOpacity (ipHE_LIKE, nelem, ipLoSing, ipHiSing);
	}
	if( !iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).hasEmis() )
	{
		iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).AddLine2Stack();
		iso_setRedisFun(ipHE_LIKE, nelem, ipLoTrip, ipHiTrip);
		iso_setOpacity (ipHE_LIKE, nelem, ipLoTrip, ipHiTrip);
	}
	if( !iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoTrip).hasEmis() )
	{
		iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoTrip).AddLine2Stack();
		iso_setRedisFun(ipHE_LIKE, nelem, ipLoTrip, ipHiSing);
		iso_setOpacity (ipHE_LIKE, nelem, ipLoTrip, ipHiSing);
	}
	if( !iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoSing).hasEmis() )
	{
		iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoSing).AddLine2Stack();
		iso_setRedisFun(ipHE_LIKE, nelem, ipLoSing, ipHiTrip);
		iso_setOpacity (ipHE_LIKE, nelem, ipLoSing, ipHiTrip);
	}

	iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoSing).Emis().Aul() = (realnum)Ass;
	iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoTrip).Emis().Aul() = (realnum)Att;
	iso_sp[ipHE_LIKE][nelem].trans(ipHiSing,ipLoTrip).Emis().Aul() = (realnum)Ast;
	iso_sp[ipHE_LIKE][nelem].trans(ipHiTrip,ipLoSing).Emis().Aul() = (realnum)Ats;

	return;
}

/*ritoa converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
STATIC double ritoa(long li, long lf, long nelem, double k, double RI2)
{
	/*  Variables are as follows:               */
	/*  lg = larger of li and lf                */
	/*  fmean = mean oscillator strength        */
	/*    for a given level.                    */
	/*  mu = reduced mass of optical electron.  */
	/*  EinsteinA = Einstein emission coef.     */
	/*  w = angular frequency of transition.    */
	/*  RI2_cm = square of rad. int. in cm^2.   */
	long lg;
	double fmean,mu,EinsteinA,w,RI2_cm;

	DEBUG_ENTRY( "ritoa()" );

	mu = ELECTRON_MASS/(1+ELECTRON_MASS/(dense.AtomicWeight[nelem]*ATOMIC_MASS_UNIT));

	w = 2. * PI * k * SPEEDLIGHT;

	RI2_cm = RI2 * BOHR_RADIUS_CM * BOHR_RADIUS_CM;

	lg = max(lf, li);

	fmean = 2.0*mu*w*lg*RI2_cm/((3.0*H_BAR) * (2.0*li+1.0));

	EinsteinA = TRANS_PROB_CONST*k*k*fmean;

	/* ASSERT( EinsteinA > SMALLFLOAT ); */

	return EinsteinA;
}

realnum helike_transprob( long nelem, long ipHi, long ipLo )
{
	double Aul;
	double Enerwn, n_eff_hi, n_eff_lo;
	long ipISO = ipHE_LIKE;

	DEBUG_ENTRY( "helike_transprob()" );

	Enerwn = iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN();
	n_eff_hi = N_(ipHi) - helike_quantum_defect(nelem,N_(ipHi),L_(ipHi),S_(ipHi),J_(ipHi));
	n_eff_lo = N_(ipLo) - helike_quantum_defect(nelem,N_(ipLo),L_(ipLo),S_(ipLo),J_(ipLo));

	if( ipHi >= iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max )
	{
		if( ipLo >= iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max )
		{
			/* Neither upper nor lower is resolved...Aul()'s are easy.	*/
			Aul = helike_transprob_collapsed_to_collapsed( nelem, N_(ipHi), N_(ipLo), Enerwn );
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);

			ASSERT( Aul > 0.);
		}
		else 
		{
			Aul = helike_transprob_collapsed_to_resolved( nelem, N_(ipHi), N_(ipLo), L_(ipLo), S_(ipLo), J_(ipLo), Enerwn );
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.01f,0.01f);
			ASSERT( Aul > 0.);
		}
	}
	else
	{
		realnum error1, error2;
		/* Both levels are resolved...do the normal bit.	*/
		if( Enerwn < 0. )
		{
			Aul = he_1trans( nelem, -1.*Enerwn,  n_eff_lo,
				L_(ipLo), S_(ipLo), ipLo-3, n_eff_hi, L_(ipHi), S_(ipHi), ipHi-3, &error1, &error2);
		}
		else
		{
			Aul = he_1trans( nelem, Enerwn,  n_eff_hi,
				L_(ipHi), S_(ipHi), ipHi-3, n_eff_lo, L_(ipLo), S_(ipLo), ipLo-3, &error1, &error2);
		}
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.01f,0.01f);
	}

	return (realnum)Aul;
}

/* This routine is pure infrastructure and bookkeeping, no physics. */
void HelikeTransProbSetup( void )
{

	const int chLine_LENGTH = 1000;
	char chLine[chLine_LENGTH];

	FILE *ioDATA;
	bool lgEOL;

	long nelem, ipLo, ipHi, i, i1, i2, i3;

	DEBUG_ENTRY( "HelikeTransProbSetup()" );

	TransProbs = (double ***)MALLOC(sizeof(double **)*(unsigned)LIMELM );

	for( nelem=ipHELIUM; nelem < LIMELM; ++nelem )
	{

		TransProbs[nelem] = (double**)MALLOC(sizeof(double*)*(unsigned)(MAX_TP_INDEX+1) );

		for( ipLo=ipHe1s1S; ipLo <= MAX_TP_INDEX;++ipLo )
		{
			TransProbs[nelem][ipLo] = (double*)MALLOC(sizeof(double)*(unsigned)MAX_TP_INDEX );
		}
	}

	/********************************************************************/
	/*************** Read in data from he_transprob.dat	*****************/
	if( trace.lgTrace )
		fprintf( ioQQQ," HelikeTransProbSetup opening he_transprob.dat:");

	ioDATA = open_data( "he_transprob.dat", "r" );

	/* check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " HelikeTransProbSetup could not read first line of he_transprob.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	i1 = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	i2 = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	if( i1 !=TRANSPROBMAGIC || i2 != N_HE1_TRANS_PROB )
	{
		fprintf( ioQQQ, 
			 " HelikeTransProbSetup: the version of he_transprob.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			 " HelikeTransProbSetup: I expected to find the number %i %i and got %li %li instead.\n" ,
			 TRANSPROBMAGIC, N_HE1_TRANS_PROB, i1, i2 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* Initialize TransProbs[nelem][ipLo][ipHi] before filling it in.	*/
	for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
	{	
		for( ipHi=0; ipHi <= MAX_TP_INDEX; ipHi++ )
		{
			for( ipLo=0; ipLo < MAX_TP_INDEX; ipLo++ )
			{
				TransProbs[nelem][ipHi][ipLo] = -1.;
			}
		}
	}

	for( ipLo=1; ipLo <= N_HE1_TRANS_PROB; ipLo++ )
	{
		char *chTemp;

		/* get next line image */
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();

		while( chLine[0]=='#' )
		{
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				BadRead();
		}

		i3 = 1;
		i1 = (long)FFmtRead(chLine,&i3,sizeof(chLine),&lgEOL);
		i2 = (long)FFmtRead(chLine,&i3,sizeof(chLine),&lgEOL);
		/* check that these numbers are correct */
		if( i1<0 || i2<=i1 )
		{
			fprintf( ioQQQ, " HelikeTransProbSetup detected insanity in he_transprob.dat.\n");
			cdEXIT(EXIT_FAILURE);
		}

		chTemp = chLine;

		/* skip over 2 tabs to start of data */
		for( i=0; i<1; ++i )
		{
			if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
			{
				fprintf( ioQQQ, " HelikeTransProbSetup could not init he_transprob\n" );
				cdEXIT(EXIT_FAILURE);
			}
			++chTemp;
		}

		/* now read in the data */
		for( nelem = ipHELIUM; nelem < LIMELM; nelem++ )
		{
			if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
			{
				fprintf( ioQQQ, " HelikeTransProbSetup could not scan he_transprob\n" );
				cdEXIT(EXIT_FAILURE);
			}
			++chTemp;

			sscanf( chTemp , "%le" , &TransProbs[nelem][i2][i1] );

			/** \todo	2	this test is out of place, where should it go? */
			if( lgEOL )
			{
				fprintf( ioQQQ, " HelikeTransProbSetup detected insanity in he_transprob.dat.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* check that ending magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " HelikeTransProbSetup could not read last line of he_transprob.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	i1 = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	i2 = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	if( i1 !=TRANSPROBMAGIC || i2 != N_HE1_TRANS_PROB )
	{
		fprintf( ioQQQ, 
			 " HelikeTransProbSetup: the version of he_transprob.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			 " HelikeTransProbSetup: I expected to find the number %i %i and got %li %li instead.\n" ,
			 TRANSPROBMAGIC, N_HE1_TRANS_PROB, i1, i2 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* close the data file */
	fclose( ioDATA );
	return;
}
