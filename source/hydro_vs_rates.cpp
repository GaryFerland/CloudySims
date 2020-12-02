/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* CS_VS80 compute thermally averaged collision strength for collisional deexcitation 
 * of hydrogenic atoms, from
 * >>refer	H1	collision	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
 *hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients */
/*Hion_coll_ioniz_ratecoef calculate hydrogenic ionization rates for ions 
 * with all n, and Z*/
#include "cddefines.h"
#include "dense.h"
#include "phycon.h"
#include "iso.h"
#include "hydro_vs_rates.h"
#include "thirdparty.h"
#include "lines_service.h"
#include "integrate.h"

STATIC double hydro_vs_coll_str( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul, long nelem, long Collider, double energy );

namespace {
	class my_Integrand : public std::unary_function<double, double> 
	{
	public:
		long nelem, Collider;
		double Aul;
		double temp;
		long nHi, gHi, nLo, gLo;
		double IP_Ryd_Hi, IP_Ryd_Lo;
		
		double operator() (double EOverKT) const
			{
				double col_str = hydro_vs_coll_str( nHi, gHi, IP_Ryd_Hi, nLo, gLo, IP_Ryd_Lo, 
					Aul, nelem, Collider, EOverKT * EVRYD * temp/TE1RYD );
				return exp( -1.*EOverKT ) * col_str;
			}
	};
}

/*
 Neither of these rates can be modified to account for impact by non-electrons because they 
 are fits to thermally-averaged rates for electron impact...it can not be unravelled to 
 expose a projectile energy that can then be scaled to account for other projectiles.
 Instead, we have included the original cross section formula (eq 14) from 
 Vriens & Smeets (1980) below.
*/

/* VS80 stands for Vriens and Smeets 1980 */
/* This routine calculates thermally-averaged collision strengths. */
double CS_VS80( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul, long nelem, long Collider, double temp )
{
	double coll_str;

	if( Collider == ipELECTRON )
	{
		coll_str = hydro_vs_deexcit( nHi, gHi, IP_Ryd_Hi, nLo, gLo, IP_Ryd_Lo, Aul );
	}
	else
	{
		/* evaluate collision strength, two options,
		 * do thermal average (very slow) if set with
		 * SET COLLISION STRENGTH AVERAGE command,
		 * default just do point at kT
		 * tests show no impact on test suite, much faster */
		if( iso_ctrl.lgCollStrenThermAver )
		{
			my_Integrand func;

			func.nHi = nHi;
			func.gHi = gHi;
			func.IP_Ryd_Hi = IP_Ryd_Hi;
			func.nLo = nLo;
			func.gLo = gLo;
			func.IP_Ryd_Lo = IP_Ryd_Lo;
			func.nelem = nelem;
			func.temp = temp;
			func.Collider = Collider;
			func.Aul = Aul;

			Integrator<my_Integrand,Gaussian32> VS80( func );
			/* do average over Maxwellian */
			coll_str = VS80.sum( 0., 1. );
			coll_str += VS80.sum( 1., 10. );
		}
		else
		{
			/* the argument to the function is equivalent to evaluating the function at
			 * hnu = kT */
			coll_str = hydro_vs_coll_str( nHi, gHi, IP_Ryd_Hi, nLo, gLo, IP_Ryd_Lo, 
					Aul, nelem, Collider, EVRYD*temp/TE1RYD );
		}
	}

	ASSERT( coll_str >= 0. );
	return coll_str;
}

/*hydro_vs_coll_str compute collision strength for collisional deexcitation for hydrogen atom, 
 * from Vriens and Smeets */
STATIC double hydro_vs_coll_str( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul, long nelem, long Collider, double energy )
{
	DEBUG_ENTRY( "hydro_vs_coll_str()" );

	// number of colliders is 4 in def, should be macro
	ASSERT( Collider >= 0.&& Collider <4 );
	double reduced_mass = dense.AtomicWeight[nelem]*colliders.list[Collider].mass_amu/
		(dense.AtomicWeight[nelem]+colliders.list[Collider].mass_amu)*ATOMIC_MASS_UNIT;

	/* This comes from equations 14,15, and 16 of Vriens and Smeets. */
	/* >>refer he-like cs Vriens, L. \& Smeets, A. H. M. Phys. Rev. A, 22, 940 */ 
	/* Computes the Vriens and Smeets
	 * rate coeff. for collisional dexcitation between any two levels of H.
	 * valid for all nhigh, nlow
	 * at end converts this excitation rate to collision strength	*/

	double n = (double)nHi;
	double p = (double)nLo;

	double g_n = (double)gHi;
	double g_p = (double)gLo;

	double ryd = EVRYD;
	double s = fabs(n-p);
	ASSERT( s > 0. );

	double Epn = EVRYD * (IP_Ryd_Lo - IP_Ryd_Hi);
	double Epi = EVRYD * IP_Ryd_Lo;

	/* This is an absorption oscillator strength. */
	double abs_osc_str = GetGF( Aul, Epn*RYD_INF/EVRYD, g_n)/g_p;
	double rEpn = 1./Epn;
	double Apn = 2.*ryd*rEpn*abs_osc_str;

	double rp = 1./p;
	double bp = rp*(1.4*log(p) - .7 + rp*(- .51 + rp*(1.16*rp - 0.55*rp)));

	double Bpn = 4.*ryd*ryd/pow3(n)*rEpn*rEpn*(1. + Epi*rEpn*(4./3. + bp*Epi*rEpn));

	double delta = exp(-Bpn/Apn) - 0.4*Epn/ryd;


	/* Scale the energy to get an equivalent electron energy. */
	energy *= colliders.list[ipELECTRON].mass_amu/colliders.list[Collider].mass_amu;

	double cross_section;
	/* cross section in units of PI*a_o^2 */
	if( energy/2./ryd+delta <= 0 /*|| energy < Epn*/ )
	{
		cross_section = 0.;
	}
	else
	{
		double gamma = ryd*(8. + 23.*POW2(s*rp))*s*s/
			( s*s*(8. + 1.1*n*s) + 0.8 + 0.4*powpq(n*s,3,2)*fabs(s-1.0) );
		cross_section = 2.*ryd/(energy + gamma)*(Apn*log(energy/2./ryd+delta) + Bpn);
		cross_section = MAX2( cross_section, 0. );
	}

	/* convert to collision strength */
	double coll_str = ConvCrossSect2CollStr( cross_section * PI*BOHR_RADIUS_CM*BOHR_RADIUS_CM, g_n, energy/EVRYD, reduced_mass );

	ASSERT( coll_str >= 0. );

	return( coll_str );
}

/*hydro_vs_coll_recomb generate hydrogenic collisional recombination rate coefficients */
double hydro_vs_coll_recomb( double ionization_energy_Ryd, double Te, double stat_level, double stat_ion )
{
	double coef, 
	  denom, 
	  epi, 
	  t_eV;

	DEBUG_ENTRY( "hydro_vs_coll_recomb()" );

	/* This is equation 9 of
	 * >>refer	H1	collision recomb	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940 */

	/* this is kT in eV */
	t_eV = Te/EVDEGK;

	/* this is the ionization energy relative to kT, dimensionless */
	epi = ionization_energy_Ryd * EVRYD / t_eV;

	/* this is the denominator of equation 8 of VS80. */
	denom = pow(epi,2.33) + 4.38*pow(epi,1.72) + 1.32*epi;

	/* this is equation 9 of VS80 */
	coef = 3.17e-27 / pow3(t_eV) * stat_level / stat_ion / denom;

	ASSERT( coef >= 0. );
	return( coef );
}


/*hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients */
double hydro_vs_ioniz( double ionization_energy_Ryd, double Te )
{
	double coef, 
	  denom, 
	  epi, 
	  t_eV;

	DEBUG_ENTRY( "hydro_vs_ioniz()" );

	/* a function written to calculate the rate coefficients 
	 * for hydrogen collisional ionizations from
	 * Jason Ferguson, summer 94
	 * valid for all n
	 * >>refer	H1	collision	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
	 * */

	/* this is kT in eV */
	t_eV = Te/EVDEGK;

	/* this is the ionization energy relative to kT, dimensionless */
	epi = ionization_energy_Ryd * EVRYD / t_eV;

	/* this is the denominator of equation 8 of VS80. */
	denom = pow(epi,2.33) + 4.38*pow(epi,1.72) + 1.32*epi;

	/* this is equation 8 of VS80 */
	coef = 9.56e-6 / powpq(t_eV,3,2) * dsexp( epi ) / denom;

	ASSERT( coef >= 0. );
	return( coef );
}

/*Hion_coll_ioniz_ratecoef calculate hydrogenic ionization rates for all n, and Z*/
double Hion_coll_ioniz_ratecoef(
		/* the isoelectronic sequence */
		long int ipISO ,
		/* element, >=1 since only used for ions 
		 * nelem = 1 is helium the least possible charge */
		long int nelem,
		/* principal quantum number, > 1
		 * since only used for excited states */
		long int n,
		double ionization_energy_Ryd,
		double Te )
{
	double H, 
	  HydColIon_v, 
	  Rnp, 
	  charge,
	  chim, 
	  eone, 
	  etwo, 
	  efour,
	  g, 
	  rate, 
	  rate2, 
	  boltz,
	  t1, 
	  t2, 
	  t3, 
	  t4, 
	  tev, 
	  xn, 
	  y;
	static const double arrH[4] = {1.48,3.64,5.93,8.32};
	static const double arrRnp[8] = {2.20,1.90,1.73,1.65,1.60,1.56,1.54,1.52};
	static const double arrg[10] = {0.8675,0.932,0.952,0.960,0.965,0.969,0.972,0.975,0.978,0.981};

	static double small = 0.;

	DEBUG_ENTRY( "Hion_coll_ioniz_ratecoef()" );
	/*calculate hydrogenic ionization rates for all n, and Z
	 * >>refer	HI	cs	Allen 1973, Astro. Quan. for low Te.
	 * >>refer	HI	cs	Sampson and Zhang 1988, ApJ, 335, 516 for High Te.
	 * */

	charge = nelem - ipISO;
	/* this routine only for ions, nelem=0 is H, nelem=1 he, etc */
	ASSERT( charge > 0);
	ASSERT( n>1 );

	if( n > 4 )
	{
		H = 2.15*n;
	}
	else
	{
		H = arrH[n-1];
	}

	if( n > 8 )
	{
		Rnp = 1.52;
	}
	else
	{
		Rnp = arrRnp[n-1];
	}

	if( n > 10 )
	{
		g = arrg[9];
	}
	else
	{
		g = arrg[n-1];
	}

	tev = EVRYD/TE1RYD*Te;
	xn = (double)n;
	chim = EVRYD * ionization_energy_Ryd;
	y = chim/tev;
	boltz = dsexp( chim/tev );

	eone = e1(y);
	etwo = e2(y);
	efour = expn(4,y);

	t1 = 1./xn*eone;
	t2 = 1./xn*efour;
	t3 = 3.*H/xn/(3. - Rnp)*(y*etwo - 2.*y*eone + boltz);
	t4 = 3.36*y*(eone - etwo);
	rate = 7.69415e-9*sqrt(Te)*9.28278e-3*powi(xn/(charge+1),4)*g*y;
	rate *= t1 - t2 + t3 + t4;
	rate2 = 2.1e-8*sqrt(Te)/chim/chim*dsexp(2.302585*5040.*chim/Te);

	/* don't let the rates go negative */
	rate = MAX2(rate,small);
	rate2 = MAX2(rate2,small);

	/* Take the lowest of the two, they fit nicely together... */
	/* >>chng 10 Sept 02, sometimes one of these is zero and the other is positive.
	 * in that case take the bigger one.	*/
	if( rate==0. || rate2==0. )
		HydColIon_v = MAX2(rate,rate2);
	else
		HydColIon_v = MIN2(rate,rate2);

	ASSERT( HydColIon_v >= 0. );
	return( HydColIon_v );
}

/*hydro_vs_deexcit compute collisional deexcitation rates for hydrogen atom, 
 * from Vriens and Smeets 1980 */
double hydro_vs_deexcit( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul )
{
	DEBUG_ENTRY( "hydro_vs_deexcit()" );

	double kT_eV = EVRYD * phycon.te / TE1RYD;

	/* This comes from equations 24 of Vriens and Smeets. */
	/* >>refer he-like cs Vriens, L. \& Smeets, A. H. M. Phys. Rev. A, 22, 940 */ 
	/* Computes the Vriens and Smeets
	 * rate coeff. for collisional dexcitation between any two levels of H.
	 * valid for all nhigh, nlow
	 * at end converts this excitation rate to collision strength	*/

	double n = (double)nLo;
	double p = (double)nHi;

	ASSERT( n!=p );

	double g_n = (double)gLo;
	double g_p = (double)gHi;

	double ryd = EVRYD;
	double s = fabs(n-p);

	double Enp = EVRYD * (IP_Ryd_Lo - IP_Ryd_Hi);
	double rEnp = 1./Enp;
	double Eni = EVRYD * IP_Ryd_Hi;

	ASSERT( Enp > 0. );

	/* This is an absorption oscillator strength. */
	double abs_osc_str = GetGF( Aul, Enp*RYD_INF/EVRYD, g_p)/g_n;
	double Anp = 2.*ryd*rEnp*abs_osc_str;

	double rn = 1./n;
	double bn = rn*( 1.4*log(n) - .7 +rn*(- .51 + rn*(1.16 - 0.55*rn)));
	
	double Bnp = 4.*ryd*ryd/pow3(p)*rEnp*rEnp*(1. + Eni*rEnp*(4./3. + bn*Eni*rEnp));

	double delta_np = exp(-Bnp/Anp) + 0.06 * s*s / (p*n*n);

	double rate;

	if( 0.3*kT_eV/ryd+delta_np <= 0 )
	{
		rate = 0.;
	}
	else
	{
		double Gamma_np = ryd*log(1. + n*n*n*kT_eV/ryd) * (3. + 11.*s*s*rn*rn) * s * s /
			( s*s* (6. + 1.6*p*s) + 0.3 + 0.8*powpq(p*s,3,2)*fabs(s-0.6) );

		rate = 1.6E-7 * sqrt(kT_eV) * g_n / (g_p * ( kT_eV + Gamma_np ) ) *
			( Anp * log(0.3*kT_eV/ryd + delta_np) + Bnp );
	}

	double col_str = rate / COLL_CONST * phycon.sqrte * g_p;

	return col_str ;
}

