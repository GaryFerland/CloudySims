/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HeRecom - do recomb coef for He, called by HeLike */
/*cross_section - calculates the photoionization cross_section for a given level and photon energy*/
/*radrecomb - calculates radiative recombination coefficients. */
/*He_cross_section returns cross section (cm^-2), 
 * given EgammaRyd, the photon energy in Ryd,
 * ipLevel, the index of the level, 0 is ground, 3 within 2 3P,
 * nelem is charge, equal to 1 for Helium 
 * this is a wrapper for cross_section */
/*Recomb_Seaton59 - find recombination for given n,using Seaton 59 approximation.
 * The following three are needed by Recomb_Seaton59:
 *     ExponentialInt
 *     X1Int
 *     X2Int	*/

#include "cddefines.h" 
#include "helike_recom.h"
#include "hydro_bauman.h"
#include "iso.h"
#include "thirdparty.h"
#include "atmdat.h"
#include "freebound.h"
#include "integrate.h"

/* The three of these are used in the computation of Recomb_Seaton59	*/
STATIC double ExponentialInt( double v );
STATIC double X1Int( double u );
STATIC double X2Int( double u );

STATIC double cross_section(double EgammaRyd, double EthRyd, long nelem, long n, long l, long s);
STATIC double GetHS98CrossSection( long n, long l, long s, double EgammaRyd );

static double Xn_S59; 

/*He_cross_section returns cross section (cm^-2), 
 * given EgammaRyd, the photon energy in Ryd,
 * ipLevel, the index of the level, 0 is ground, 3 within 2 3P,
 * nelem is charge, equal to 1 for Helium 
 * this is a wrapper for cross_section */
double He_cross_section( double EgammaRyd , double EthRyd, long n, long l, long S, long nelem )
{
	// get cross section in megabarns
	double cs = cross_section( EgammaRyd, EthRyd, nelem, n, l, S );

	// rescale low-lying He values to Hummer & Storey 98, Table 1 Extrapolated
	if( nelem==ipHELIUM && n <=5 && l<=2 )
	{
		static const double rescaled[31] = {
			7.394, 
			5.485,  9.219, 15.985, 15.985, 15.985, 13.504,
			8.018, 14.417, 28.501, 18.486, 18.132, 27.009,
			10.721, 20.235, 41.568, 36.717, 35.766, -1.000, -1.000, 41.787, // };
			13.527, 26.539, 55.692, 55.010, 53.514, -1.000, -1.000, -1.000, -1.000, 58.120 };
		long ipLev = iso_sp[ipHE_LIKE][nelem].QuantumNumbers2Index[n][l][S]; 
		ASSERT( rescaled[ipLev] > 0. ); 
		cs *= rescaled[ipLev]/cross_section( EthRyd, EthRyd, nelem, n, l, S ); 
	}

	// convert to cm^-2
	return cs * (1.e-18);
}

/*cross_section calculates the photoionization cross_section for a given level and photon energy
 * this routine returns megabarns */
STATIC double cross_section(double EgammaRyd, double EthRyd, long nelem, long n, long l, long S)
{
	/* These fit parameters (E0, sigma, y_a, P, y_w, yzero, and yone) all come from the following work: */
	/* >>refer	He	pcs	Verner, D. A., Verner, E. M., \& Ferland , G. J. 1996,
	 * >>refercon	Atomic Data and Nuclear Data Tables, Vol. 64, p.1 */
	static const double E0[29] = {
	1.36E+01,2.01E+01,1.76E+01,3.34E+01,4.62E+01,6.94E+01,8.71E+01,1.13E+02,1.59E+02,2.27E+02,
	2.04E+02,2.74E+02,2.75E+02,3.38E+02,4.39E+02,4.17E+02,4.47E+02,5.18E+02,6.30E+02,6.27E+02,
	8.66E+02,7.67E+02,9.70E+02,9.66E+02,1.06E+03,1.25E+03,1.35E+03,1.43E+03,1.56E+03};
	static const double sigma[29] = {
	9.49E+02,3.20E+02,5.46E+02,2.85E+02,2.34E+02,1.52E+02,1.33E+02,1.04E+02,6.70E+01,4.00E+01,
	6.14E+01,4.04E+01,4.75E+01,3.65E+01,2.45E+01,3.14E+01,3.11E+01,2.59E+01,1.94E+01,2.18E+01,
	1.23E+01,1.76E+01,1.19E+01,1.31E+01,1.20E+01,9.05E+00,8.38E+00,8.06E+00,7.17E+00};
	static const double y_a[29] = {
	1.47E+00,7.39E+00,1.72E+01,2.16E+01,2.18E+01,2.63E+01,2.54E+01,2.66E+01,3.35E+01,5.32E+01,
	2.78E+01,3.57E+01,2.85E+01,3.25E+01,4.41E+01,3.16E+01,3.04E+01,3.28E+01,3.92E+01,3.45E+01,
	5.89E+01,3.88E+01,5.35E+01,4.83E+01,5.77E+01,6.79E+01,7.43E+01,7.91E+01,9.10E+01};
	static const double P[29] = {
	3.19E+00,2.92E+00,3.16E+00,2.62E+00,2.58E+00,2.32E+00,2.34E+00,2.26E+00,2.00E+00,1.68E+00,
	2.16E+00,1.92E+00,2.14E+00,2.00E+00,1.77E+00,2.04E+00,2.09E+00,2.02E+00,1.86E+00,2.00E+00,
	1.62E+00,1.93E+00,1.70E+00,1.79E+00,1.72E+00,1.61E+00,1.59E+00,1.58E+00,1.54E+00};
	static const double y_w[29] = 
	{2.039,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	static const double yzero[29] =
	{0.4434,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	static const double yone[29] =
	{2.136,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	double pcs,Egamma,y,F,x;
	double rel_photon_energy;
	
	Egamma = EgammaRyd * EVRYD;

	/* >>chng 02 apr 24, more protection against calling with too small an energy */
	/* evaluating H-like photo cs at He energies, may be below threshold -
	 * prevent this from happening */
	rel_photon_energy = EgammaRyd / EthRyd;
	rel_photon_energy = MAX2( rel_photon_energy , 1. + FLT_EPSILON*2. );

	long s=0;//portland group 11.2 trips without =0, does not recognized that TotalInsanity does not return
	if( S == 1 )
		s=0;
	else if( S == 3 )
		s=1;
	else
		TotalInsanity();

	if( nelem==ipHELIUM && n<=25 && l<=4 )
	{
		// use Hummer & Storey 1998 cross-sections
		pcs = GetHS98CrossSection( n, l, s, EgammaRyd );
	}
	else if( nelem==ipHELIUM && n>25 && l<=2 )
	{
		static const double scale[3][2] = {
			{1.4673,1.3671},
			{1.5458,1.5011},
			{1.4912,1.5144}};

		long ipLev = iso_sp[ipHE_LIKE][nelem].QuantumNumbers2Index[25][l][S];
		double EthRyd_25 = iso_sp[ipHE_LIKE][nelem].fb[ipLev].xIsoLevNIonRyd;
		pcs = GetHS98CrossSection( 25, l, s, EthRyd_25 * rel_photon_energy );
		pcs *= pow((double)n/25., scale[l][s]);
	}
	else if( n==1 )
	{
		/* >>refer Helike	PCS	Verner, D.A., Ferland, G.J., Korista, K.T., & Yakovlev, D.G.
		 * >>refercon	1996a, ApJ 465,487	*/
		/* All helike (but not Helium itself) ground state cross-sections calculated here.	*/
		x = Egamma/E0[nelem-1] - yzero[nelem-1];
		y = sqrt(x*x + yone[nelem-1]*yone[nelem-1]);
		F = ((x-1)*(x-1)+y_w[nelem-1]*y_w[nelem-1])
			* pow(y,0.5*P[nelem-1]-5.5) * pow((1+sqrt(y/y_a[nelem-1])),-P[nelem-1]);
		pcs = sigma[nelem-1]*F;
	}
	else if( nelem>=ipLITHIUM && nelem<=ipCALCIUM && n<11 && OP_Helike_NumPts[nelem][n][l][s]>0 )
	{
		// use TOPbase cross-sections
		long numDataPoints = OP_Helike_NumPts[nelem][n][l][s];
		ASSERT( numDataPoints > 0 );
		ASSERT( OP_Helike_Xsectn[nelem][n][l][s] != NULL );

		if( EgammaRyd < OP_Helike_Energy[nelem][n][l][s][numDataPoints-1] )
		{
			pcs = linint( OP_Helike_Energy[nelem][n][l][s], OP_Helike_Xsectn[nelem][n][l][s], numDataPoints, EgammaRyd );
		}
		else
		{
			// use a E^-3 tail 
			pcs = OP_Helike_Xsectn[nelem][n][l][s][numDataPoints-1] * POW3( OP_Helike_Energy[nelem][n][l][s][numDataPoints-1]/EgammaRyd );
		}
	}
	else 
	{
		/* To everything else we apply a hydrogenic routine.	*/
		pcs = (1.e18)*H_photo_cs(rel_photon_energy , n, l, nelem);
	}

	ASSERT( pcs > 0. && pcs < 1.E10 );

	return pcs;
}

STATIC double GetHS98CrossSection( long n, long l, long s, double EgammaRyd )
{
	double pcs;
	ASSERT( n<=25 );
	ASSERT( l<=4 );
	ASSERT( s==0 || s==1 );

	// use Hummer & Storey 1998 cross-sections
	if( EgammaRyd < HS_He1_Energy[n][l][s][NUM_HS98_DATA_POINTS-1] )
	{
		pcs = linint( HS_He1_Energy[n][l][s], HS_He1_Xsectn[n][l][s], NUM_HS98_DATA_POINTS, EgammaRyd );
	}
	else
	{
		// use a E^-3 tail 
		pcs = HS_He1_Xsectn[n][l][s][NUM_HS98_DATA_POINTS-1] * POW3( HS_He1_Energy[n][l][s][NUM_HS98_DATA_POINTS-1]/EgammaRyd );
	}

	return pcs;
}

#if 1
/* >>refer	He-like	RR	Seaton, M.J. 1959, MNRAS 119, 81S */
double Recomb_Seaton59( long nelem, double temp, long n)
{
	double lambda = TE1RYD * nelem * nelem / temp;
	/* smallest x ever used here should be lowest Z, highest T, highest n...
	 * using helium, logt = 10., and n = 1000, we get xmin = 1.5789E-11.	*/
	double x = lambda / n / n;
	double AlphaN;
	double SzeroOfX = 0.;
	double SoneOfX = 0.;
	double StwoOfX = 0.;
	double SnOfLambda = 0.;
	double lowerlimit, upperlimit, step;

	fixit("the variant below should be faster and more accurate, but needs more testing");
	Xn_S59 = x;

	/* Equation 12	*/
	lowerlimit = x;
	step = 3. * x;
	upperlimit = lowerlimit + step;
	SzeroOfX = qg32( lowerlimit, upperlimit, ExponentialInt);

	do
	{
		lowerlimit = upperlimit;
		step *= 2;
		upperlimit = lowerlimit + step;
		SzeroOfX += qg32( lowerlimit, upperlimit, ExponentialInt);
	} while ( upperlimit < 20. );

	/* This must be placed inside integral...too big to be 
	 * handled separately.	
	SzeroOfX *= exp( x );	*/

	/* Equations 13 and 14 */
	lowerlimit = 0.;
	step = 0.5;
	upperlimit = lowerlimit + step;
	SoneOfX = qg32(lowerlimit, upperlimit, X1Int);
	StwoOfX = qg32(lowerlimit, upperlimit, X2Int);

	do
	{
		lowerlimit = upperlimit;
		step *= 2;
		upperlimit = lowerlimit + step;
		SoneOfX += qg32( lowerlimit, upperlimit, X1Int);
		StwoOfX += qg32( lowerlimit, upperlimit, X2Int);
	} while ( upperlimit < 200. );

	SoneOfX *= 0.1728 * powpq( x, 1, 3 );
	StwoOfX *= -0.0496 * powpq( x, 2, 3 );

	/* Equation 11	*/
	SnOfLambda = SzeroOfX + powpq(1./lambda, 1, 3)*SoneOfX + powpq(1./lambda, 2, 3)*StwoOfX;

	if( false )
	{
		double SSzeroOfX = e1_scaled(x);
		double gamma13 = tgamma(1./3.);
		double gamma23 = PI2/(sqrt(3.)*gamma13);
		double x13 = cbrt(x);
		double x23 = pow2(x13);
		double SSoneOfX = 0.1728*((3.*x+1)*igamc_scaled(1./3.,x)*gamma13 - 3.*x13);
		double SStwoOfX = -0.0496*(((1.5*x+2.)*x+1.)*igamc_scaled(2./3.,x)*gamma23 - 1.5*x23*(x+1.));
		double SSnOfLambda = SSzeroOfX + powpq(1./lambda, 1, 3)*SSoneOfX + powpq(1./lambda, 2, 3)*SStwoOfX;

		dprintf( ioQQQ, "%e %e %e %e %e old %e new %e\n", x, SzeroOfX/SSzeroOfX-1., SoneOfX/SSoneOfX-1.,
			 StwoOfX/SStwoOfX-1., SnOfLambda/SSnOfLambda-1., SnOfLambda, SSnOfLambda );
	}

	AlphaN = 5.197E-14 * nelem * powpq(x, 3, 2) * SnOfLambda;

	return AlphaN;
}

STATIC double ExponentialInt( double v )
{
	double Integrand;

	Integrand = exp( -v + Xn_S59) / v;

	return Integrand;
}

STATIC double X1Int( double u )
{
	double Integrand;

	Integrand = powpq(1./(u + 1.), 5, 3) * (u - 1.) * exp( -Xn_S59 * u );

	return Integrand;
}

STATIC double X2Int( double u )
{
	double Integrand;

	Integrand = powpq(1./(u + 1.), 7, 3) * (u*u + 4./3.*u + 1.) * exp( -Xn_S59 * u );

	return Integrand;
}
#else
/* >>refer	He-like	RR	Seaton, M.J. 1959, MNRAS 119, 81S */
double Recomb_Seaton59(long nelem, double temp, long n)
{
	double lambda = TE1RYD * nelem * nelem / temp;
	/* smallest x ever used here should be lowest Z, highest T, highest n...
	 * using helium, logt = 10., and n = 1000, we get xmin = 1.5789E-11.	*/
	double x = lambda / n / n;

	/* Equation 12	*/
	double SzeroOfX = e1_scaled(x);

	/* Equations 13 and 14 */
	double gamma13 = tgamma(1./3.);
	double gamma23 = PI2/(sqrt(3.)*gamma13);
	double x13 = cbrt(x);
	double x23 = pow2(x13);
	double SoneOfX = 0.172826*((3.*x+1)*igamc_scaled(1./3.,x)*gamma13 - 3.*x13);
	double StwoOfX = -0.0495957*(((1.5*x+2.)*x+1.)*igamc_scaled(2./3.,x)*gamma23 - 1.5*x23*(x+1.));

	/* Equation 11	*/
	double z13 = cbrt(1./lambda);
	double SnOfLambda = (StwoOfX*z13 + SoneOfX)*z13 + SzeroOfX;

	return 5.197e-14 * nelem * powpq(x,3,2) * SnOfLambda;
}
#endif
