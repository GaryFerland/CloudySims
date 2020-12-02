/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*rt_continuum_shield_fcn computing continuum shielding due to single line */
/*conpmp local continuum pumping rate radiative transfer for all lines */
/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*DrvContPump local continuum pumping rate radiative transfer for all lines */
/*con_pump_op  routine used to get continuum pumping of lines 
 * used in DrvContPump in call to qg32 */

#include "cddefines.h"
#include "rt.h"
#include "rt_escprob.h"
#include "transition.h"
#include "thirdparty.h"
#include "integrate.h"

#include "cosmology.h"

/*conpmp local continuum pumping rate radiative transfer for all lines */
STATIC double conpmp(double tau, double damp);
STATIC double conpmp_romb( double tau, double damp);
STATIC double growth_romb( double tau, double damp);
STATIC inline double fitted( double t );

namespace {
	
/*con_pump_op  routine used to get continuum pumping of lines 
 * used in DrvContPump in call to qg32 */
class my_Integrand_con_pump_op
{
public:
	/* variable used for inward optical depth for pumping */
	realnum PumpTau;
	/* damping constant used for pumping */
	realnum damp;
	my_Integrand_con_pump_op( realnum tau, realnum damp) :
		PumpTau(tau), damp(damp)
		{}
private:
	// Do not implement
	my_Integrand_con_pump_op( void );
public:
	double operator() (double x) const
	{
		realnum v, rx = realnum(x);
		VoigtH(damp,&rx,&v,1);
		double opfun_v = sexp(PumpTau*v)*v;
		return opfun_v;
	}
};

class con_pump_op_conv
{
public:
	/* variable used for inward optical depth for pumping */
	realnum PumpTau;
	/* damping constant used for pumping */
	realnum damp;
	con_pump_op_conv( realnum tau, realnum damp) :
		PumpTau(tau), damp(damp)
		{}
private:
	// Do not implement
	con_pump_op_conv( void );
public:
	double operator() (double y) const
	{
		// Substitute y = 1/(1+x) to map integral onto a finite domain
		double x = 1./y - 1.;
		realnum v, rx = realnum(x);
		if (damp >= 0.0)
			VoigtU(damp,&rx,&v,1);
		else
			v = 1.0/(PI*(1.0+x*x));
		double opfun_v = sexp(PumpTau*v)*v;
		return opfun_v/(y*y);
	}
};

class profile_pow
{
public:
	/* power of profile function */
	int npow;
	/* damping constant used for pumping */
	realnum damp;
	profile_pow( int npow, realnum damp) :
		npow(npow), damp(damp)
		{}
private:
	// Do not implement
	profile_pow( void );
public:
	double operator() (double y) const
	{
		// Substitute y = 1/(1+x) to map integral onto a finite domain
		double x = 1./y - 1.;
		realnum v, rx = realnum(x);
		if (damp >= 0.)
		{
			VoigtU(damp,&rx,&v,1);
		}
		else
		{
			// Normalized Lorentz for negative a
			v = 1./(PI*(1.0+x*x));
		}
		double opfun_v = powi(v,npow);
		return opfun_v/(y*y);
	}
};

class growth
{
public:
	/* variable used for inward optical depth for pumping */
	realnum PumpTau;
	/* damping constant used for pumping */
	realnum damp;
	growth( realnum tau, realnum damp) :
		PumpTau(tau), damp(damp)
		{}
private:
	// Do not implement
	growth( void );
public:
	double operator() (double y) const
	{
		// Substitute y = 1/(1+x) to map integral onto a finite domain
		double x = 1./y - 1.;
		realnum v, rx = realnum(x);
		if (damp >= 0.0)
			VoigtU(damp,&rx,&v,1);
		else
			v = 1.0/(PI*(1.0+x*x));
		double opfun_v = PumpTau*v;
		if (opfun_v > 1e-6)
			opfun_v = 1-sexp(opfun_v);
		else
			opfun_v *= 1-0.5*opfun_v;
		return opfun_v/(y*y);
	}
};

}

static double shieldFederman(double tau, double damp, bool lgBug);

static double shieldRodgers(double tau, double damp);
static double growthRodgers(double tau, double damp);

/*rt_continuum_shield_fcn computing continuum shielding due to single line */
STATIC double RT_continuum_shield_fcn_point( const TransitionProxy& t, double tau )
{
	DEBUG_ENTRY( "rt_continuum_shield_fcn_point()" );
	double value = -1.;
	
	ASSERT( t.Emis().damp() > 0. );
	if( cosmology.lgDo )
	{
		if( tau < 1e-5 )
			value = (1. - tau/2.);
		else
			value = (1. - dsexp( tau ) )/ tau;
	}
	else if( rt.nLineContShield == LINE_CONT_SHIELD_PESC )
	{
		/* set continuum shielding pesc - shieding based on escape probability */
		if( t.Emis().iRedisFun() == ipPRD )
		{
			value =  esc_PRD_1side(tau,t.Emis().damp());
		}
		else if( t.Emis().iRedisFun() == ipCRD )
		{
			value =  esca0k2(tau);
		}
		else if( t.Emis().iRedisFun() == ipCRDW )
		{
			value =  esc_CRDwing_1side(tau,t.Emis().damp());
		}
		else if( t.Emis().iRedisFun() == ipLY_A )
		{
			value = esc_PRD_1side(tau,t.Emis().damp());
		}
		else
			TotalInsanity();
	}
	else if( rt.nLineContShield == LINE_CONT_SHIELD_FEDERMAN )
	{
		// Last arg is whether to use buggy terms compatible
		// with the implementation at 13.02 and previously, and 
		// as published by Federman et al.  These have
		// been shown to be incorrect, but may not be fixed in
		// other codes.  false means to use the corrected version.
		value = shieldFederman(tau,t.Emis().damp(),false);
	}
	else if( rt.nLineContShield == LINE_CONT_SHIELD_FEDERMAN_BUG )
	{
		value = shieldFederman(tau,t.Emis().damp(),true);
	}
	else if( rt.nLineContShield == LINE_CONT_SHIELD_FERLAND )
	{
		/* set continuum shielding ferland */
		value = conpmp( tau , t.Emis().damp() );
	}
	else if( rt.nLineContShield == LINE_CONT_SHIELD_RODGERS )
	{
		value = shieldRodgers(tau*SQRTPI,t.Emis().damp());
	}
	else if( rt.nLineContShield == LINE_CONT_SHIELD_INTEGRAL )
	{
		value = conpmp_romb(tau*SQRTPI,t.Emis().damp());
		value = MIN2(1.,value);
	}
	else if( rt.nLineContShield == 0 )
	{
		/* set continuum shielding none */
		value = 1.;
	}
	else
	{
		TotalInsanity();
	}

	/* the returned pump shield function must be greater than zero,
	 * and less than 1 if a maser did not occur */
	ASSERT( value>=0 && (value<=1.||tau<0.) );
	return value;
}

// Function to collect a range of options to handle finite optical
// depth effects for continuum shielding.  s1 is the shielding factor
// evaluated at the front face, s2 is the shielding factor evaluated at
// the rear face.
STATIC double avg_shield(double s1, double s2, double dTau, double tau)
{
	enum options { AVG_OLD, AVG_DEPTH, AVG_NO, AVG_BACK, AVG_ARITH, 
						AVG_GEOM, AVG_COUNT };
	const enum options opt = AVG_DEPTH;
	if (opt == AVG_OLD)
	{
		return s1*log(1.+dTau)/dTau;
	}
	else if (opt == AVG_DEPTH)
	{
		// Average (1+tau)/(1+tau+dTau) over zone...
		double dTauRel = dTau/(1.+tau);
		ASSERT( 1.+dTauRel >= 0. );
		return s1*log(1.+dTauRel)/dTauRel;
	}
	else if (opt == AVG_NO)
	{
		return s1; // no shielding
	}
	else if (opt == AVG_BACK)
	{
		return s2; // implicit
	}
	else if (opt == AVG_ARITH)
	{
		return 0.5*(s1+s2); // Arithmetic average
	}
	else if (opt == AVG_GEOM)
	{
		return sqrt(s1*s2); // Geometric average
	}
	else if (opt == AVG_COUNT)
	{
		return (s1-s2)/dTau; // photon counting estimate [?]
	}
	TotalInsanity(); // invalid option
}

STATIC double shieldFederman(double tau, double damp, bool lgBug)
{
	DEBUG_ENTRY( "shieldFederman()" );
		/* set continuum shielding Federman - this is the default */
	double core, wings, value;

	/* these expressions implement the appendix of
	 * >>refer	line	shielding	Federman, S.R., Glassgold, A.E., & 
	 * >>refercon	Kwan, J. 1979, ApJ, 227, 466 */
	/* doppler core - equation A8 */
	if( tau < 2. )
	{
		core = sexp( tau * 0.66666 );
	}
	else if( tau < 10. )
	{
		core = 0.638 * pow(tau,(realnum)-1.25f );
	}
	else if( tau < 100. )
	{
		core = 0.505 * pow(tau,(realnum)-1.15f );
	}
	else
	{
		core = 0.344 * pow(tau,(realnum)-1.0667f );
	}
	
	/* do we add damping wings? */
	wings = 0.;
	if( damp>0. )
	{
		/* equation A6 */
		double t1 = 3.02*pow(damp*1e3,-0.064 );
		double tauwing = lgBug ? tau : tau*SQRTPI;
		double u1 = sqrt(MAX2(tauwing,0.)*damp )/SDIV(t1);
		wings = damp/SDIV(t1)/sqrt( 0.78540 + POW2(u1) );
		/* add very large optical depth tail to converge this with respect
		 * to escape probabilities - if this function falls off more slowly
		 * than escape probability then upper level will become overpopulated.
		 * original paper was not intended for this regime */
		if( lgBug && tau>1e7 )
			wings *= pow( tau/1e7,-1.1 );
	}
	value = core + wings;
	/* some x-ray lines have vastly large damping constants, greater than 1.
	 * in these cases the previous wings value does not work - approximation
	 * is for small damping constant - do not let pump efficiency exceed unity
	 * in this case */
	if( tau>=0. )
		value = MIN2(1., value );
	return value;
}

double RT_continuum_shield_fcn( const TransitionProxy& t, 
				bool lgShield_this_zone, double dTau )
{
	DEBUG_ENTRY( "rt_continuum_shield_fcn()" );

	double tau = t.Emis().TauCon();
	double value = RT_continuum_shield_fcn_point(t,tau);

	if( lgShield_this_zone && dTau > 1e-3 )
	{
		if (0 && t.ipCont() == 3627)
			fprintf(ioQQQ,"?shield %ld %15.8g %15.8g %15.8g %15.8g %s\n",
					  nzone,tau,dTau,1./(1. + dTau ),
					  RT_continuum_shield_fcn_point(t,tau+dTau)/value,
					  chLineLbl(t).c_str());
		/* correction for line optical depth across zone */
		value = avg_shield(value,RT_continuum_shield_fcn_point(t,tau+dTau),
								 dTau,tau);
	}

	return value;
}

/*conpmp local continuum pumping rate radiative transfer for all lines */
STATIC double conpmp_qg32( double tau, double damp)
{
	DEBUG_ENTRY( "conpmp_qg32()" );

	my_Integrand_con_pump_op func(tau, damp);
	Integrator<my_Integrand_con_pump_op,Gaussian32> opfun(func);
	
	const double BREAK = 3.;
	
	double yinc1 = opfun.sum( 0., BREAK );
	double yinc2 = opfun.sum( BREAK, 100. );
	
	double a0 = 0.5*SQRTPI;
	return (yinc1 + yinc2)/a0;
}

/*conpmp local continuum pumping rate radiative transfer for all lines */
STATIC double conpmp_romb( double tau, double damp)
{
	DEBUG_ENTRY( "conpmp_romb()" );

	double tauint = tau;

	con_pump_op_conv func(tauint, damp);

	// Ignore underflowing values at small offsets x, i.e. y~=1
	double top=1.0;
	for (;;)
	{
		if (func(0.5*top) > 0.0)
			break;
		top *= 0.5;
	}

	class integrate::Romberg<integrate::Midpoint<con_pump_op_conv> >
		intl(integrate::Midpoint<con_pump_op_conv>(func,0.0,top));
	intl.update(1e-10);
	// 2.0 because only integrating RHS
	return 2.0*intl.sum();
}

//*conpmp local continuum pumping rate radiative transfer for all lines */
STATIC double growth_romb( double tau, double damp)
{
	DEBUG_ENTRY( "growth_romb()" );

	double tauint = tau;

	growth func(tauint, damp);

	class integrate::Romberg<integrate::Midpoint<growth> >
		intl(integrate::Midpoint<growth>(func,0.0,1.0));
	intl.update(1e-10);
	// 2.0 because only integrating RHS
	return 2.0*intl.sum();
}

/*conpmp local continuum pumping rate radiative transfer for all lines */
STATIC double conpmp( double tau, double damp)
{
	DEBUG_ENTRY( "conpmp()" );
	double conpmp_v;

	/* tau required is optical depth in center of next zone */
	/* compute pumping probability */
	if( tau <= 10+2.5*damp)
	{
		double tausc;
		if (damp < 1e-2)
		{
			tausc = tau*(1.-1.5*damp);
		}
		else
		{
			tausc = tau*(1+damp)/(1.+2.5*damp*(1.+damp));
		}
		conpmp_v = fitted(tausc);
	}
	else if( tau > 1e6 )
	{
		/* this far in winds line opacity well below electron scattering
		 * so ignore for this problem */
		conpmp_v = 0.;
	}
	else
	{
		conpmp_v = conpmp_qg32(tau, damp);
	}
	return conpmp_v;
}

/* fit to results for tau less than 10 */
inline double fitted(double t)
{
	return (0.98925439 + 0.084594094*t)/(1. + t*(0.64794212 + t*0.44743976));
}

double DrvContPump(double tau, double damp)
{
	DEBUG_ENTRY( "DrvContPump()" );
	return conpmp(tau, damp);
}

void DrivePump(double damp)
{
	double taufac = SQRTPI;
	double dampa = damp > 0.0 ? damp : 0.0;
	fprintf(ioQQQ,"#tau damp QG32 conpmp Federman Romb\n");
	for (long i=0; i<49; ++i)
	{
		double tau = 1e-4*exp10(i/4.);
		fprintf(ioQQQ,"%13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",
				  tau, damp, conpmp_qg32(tau/taufac,dampa), 
				  conpmp(tau/taufac,dampa),
				  shieldFederman(tau/taufac,dampa,true),conpmp_romb(tau,damp),
				  shieldRodgers(tau,damp)
			);
	}
	fprintf(ioQQQ,"\n");

	if (0)
	{
		int NSAMP=10;
		FILE* ioVALID = open_data("growth.dat","w");
		fprintf(ioVALID,"#tau damp growth_romb growthRodgers\n");
		for (long j=0; j<1+7*NSAMP; ++j)
		{
			double damp1 = 1e-5*exp10(j/double(NSAMP));
			for (long i=0; i<1+12*NSAMP; ++i)
			{
				double tau = 1e-4*exp10(i/double(NSAMP));
				double gromb = growth_romb(tau,damp1),
					grodg = growthRodgers(tau,damp1),
					err = safe_div(grodg,gromb,1.)-1.;
				fprintf(ioVALID,"%13.8e %13.8e %13.8e %13.8e %13.8e\n",
						  tau, damp1,gromb, grodg, err	);
			}
		}
		fclose(ioVALID);
	}

	if (0)
	{
		int NSAMP=10;
		FILE* ioVALID = open_data("shield.dat","w");
		fprintf(ioVALID,"#tau damp shield_romb shieldRodgers\n");
		for (long j=0; j<1+7*NSAMP; ++j)
		{
			double damp1 = 1e-5*exp10(j/double(NSAMP));
			for (long i=0; i<1+12*NSAMP; ++i)
			{
				double tau = 1e-4*exp10(i/double(NSAMP));
				double gromb = conpmp_romb(tau,damp1),
					grodg = shieldRodgers(tau,damp1), // shieldFederman(tau/taufac,damp1,false), 
					err = safe_div(grodg,gromb,1.)-1.;
				fprintf(ioVALID,"%13.8e %13.8e %13.8e %13.8e %13.8e\n",
						  tau, damp1,gromb, grodg, err	);
			}
		}
		fclose(ioVALID);
	}

	if (damp >= 0.0)
	{
		double tau = 1.0;
		if (damp != 0.0 && damp < 0.5)
		{
			tau = 1.0/damp;
			for (;;)
			{
				double tau0 = tau;
				tau = 1.0/(damp*log(tau));
				if (fabs(tau-tau0) < 1e-4*tau)
				{
					break;
				}
			}
		}
		for (long i=0; i<49; ++i)
		{
			double damp1 = 1e-6*damp*exp10(0.25*i);
			fprintf(ioQQQ,"%13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",
					  tau, damp1, conpmp_qg32(tau,damp1), conpmp(tau,damp1),
					  shieldFederman(tau,damp1,true),conpmp_romb(tau,damp1)
				);
		}
		fprintf(ioQQQ,"\n");
	}
	
	fprintf(ioQQQ,"Integrals of phi^n\n");
	fprintf(ioQQQ,"%13.8e", damp);
	for (int n=1; n<10; ++n)
	{
		class integrate::Romberg<integrate::Midpoint<profile_pow> >
			intg(integrate::Midpoint<profile_pow>(profile_pow(n,damp),0.0,1.0));
		intg.update(1e-10);
		fprintf(ioQQQ," %13.8e",2.0*intg.sum());
	}
	fprintf(ioQQQ,"\n");
}

// Rodgers & Williams JQSRT 14, 319 (1974)
// tau phi(x) = -k(v) m
//   a -> y = alpha_L / alphaD
//   tau -> S m / alphaD
//   x -> (nu-nu0) / alphaD
//   W = int 1-exp(-k m) dnu = alphaD int 1-exp(-phi tau) dx
// So can divide all W's by alphaD to normalize
//   WM = [WL^2 + WD^2 - (WL WD/WW)^2]^{1/2}
// stands.
//   WL = 2 pi a L ( tau / 2 pi a)
//   WD = D ( tau / sqrtpi )
//   WW = S m / alphaD = tau

static void shieldRodgersLorentz(double tau, double damp, 
											double& w, double &dw)
{
	// z = S m / 2 pi alphaL = tau / 2 sqrt(pi) y 
	double z = tau/(2 * PI * damp ); // Check scaling
	if (z < 1.98) 
	{
		static const double a[] = {9.99997697674e-1,
											-4.99882252233e-1,
											2.49005331874e-1,
											-1.00956020660e-1,
											3.13580341312e-2,
											-6.50144170817e-3,
											6.48325819427e-4};
		w = z*(a[0]+z*(
					 a[1]+z*(
						 a[2]+z*(
							 a[3]+z*(
								 a[4]+z*(
									 a[5]+z*a[6]))))));
		dw = a[0]+z*(
			2*a[1]+z*(
				3*a[2]+z*(
					4*a[3]+z*(
						5*a[4]+z*(
							6*a[5]+z*7*a[6])))));
	}
	else
	{
		// Note that this is an expansion in z^-1, not z, a typo
		// in the original paper
		static const double b[] = {7.97885095473e-1,
											-9.97555137671e-2,
											-1.97623661059e-2,
											1.05865022723e-2,
											-1.32467496350e-1,
											1.47947878333e-1};
		double sz = sqrt(z),rz = 1./z;
		w = sz*(b[0]+rz*(b[1]+rz*(b[2]+rz*(b[3]+rz*(b[4]+rz*b[5])))));
		dw = (b[0]+rz*(
					-1*b[1]+rz*(
						-3*b[2]+rz*(
							-5*b[3]+rz*(
								-7*b[4]-rz*9*b[5])))))/
			(2.*sz);
	}
	w *= 2.*PI*damp;
	//dw *= damp;
}

static void shieldRodgersDoppler(double tau,
											double& w, double &dw)
{
	// z = S m / pi^1/2 alphaD
	double z = tau/SQRTPI; // Check scaling
	if (z < 5) 
	{
		static const double c[] = {9.99998291698e-1,
											-3.53508187098e-1,
											9.60267807976e-2,
											-2.04969011013e-2,
											3.43927368627e-3,
											-4.27593051557e-4,
											3.42209457833e-5,
											-1.28380804108e-6};
		w = z*SQRTPI*(c[0]+z*(
							  c[1]+z*(
								  c[2]+z*(
									  c[3]+z*(
										  c[4]+z*(
											  c[5]+z*(
												  c[6]+z*c[7])))))));
		dw = (c[0]+z*(
					2*c[1]+z*(
						3*c[2]+z*(
							4*c[3]+z*(
								5*c[4]+z*(
									6*c[5]+z*(
										7*c[6]+z*8*c[7])))))));
	}
	else
	{
		static const double d[] = {1.99999898289e0,
											5.77491987800e-1,
											-5.05367549898e-1,
											8.21896973657e-1,
											-2.52226724530e0,
											6.10070274810e0,
											-8.51001627836e0,
											4.65351167650e0};
		double lz = log(z), slz=sqrt(lz), rlz=1./lz;
		w = slz*(d[0]+rlz*(
						d[1]+rlz*(
							d[2]+rlz*(
								d[3]+rlz*(
									d[4]+rlz*(
										d[5]+rlz*(
											d[6]+rlz*d[7]))))))); //SQRTPI;
		dw = (d[0]+rlz*(
					-1*d[1]+rlz*(
						-3*d[2]+rlz*(
							-5*d[3]+rlz*(
								-7*d[4]+rlz*(
									-9*d[5]+rlz*(
										-11*d[6]-rlz*13*d[7])))))))/
			(2.*slz*z*SQRTPI);
	}
}

// Rodgers & Williams curve of growth
static double growthRodgers(double tau, double damp)
{
	if (damp < 0)
	{
		double wl, dwl;
		shieldRodgersLorentz(tau,1.0,wl,dwl);
		return wl;
	}

	double wd, dwd;
	shieldRodgersDoppler(tau,wd,dwd);
	if (damp == 0.0)
		return wd;

	double wl, dwl;
	shieldRodgersLorentz(tau,damp,wl,dwl);
	double wls = wl*wl, wds = wd*wd, rtaus=1./(tau*tau);
	double wv = sqrt(wls+wds-wls*wds*rtaus);
	return wv;
}

// Rodgers & Williams derived shielding function
static double shieldRodgers(double tau, double damp)
{
	if (damp < 0)
	{
		double wl, dwl;
		shieldRodgersLorentz(tau,1.0,wl,dwl);
		return dwl;
	}

	double wd, dwd;
	shieldRodgersDoppler(tau,wd,dwd);
	if (damp == 0.0)
		return dwd;

	double wl, dwl;
	shieldRodgersLorentz(tau,damp,wl,dwl);
	double wls = wl*wl, wds = wd*wd, rtaus=1./(tau*tau);
	double wv = sqrt(wls+wds-wls*wds*rtaus);
	double dwv = (wl*dwl*(1.0-wds*rtaus)+
					  wd*dwd*(1.0-wls*rtaus)+
					  wls*wds*rtaus/tau)/wv;
	return dwv;
}
