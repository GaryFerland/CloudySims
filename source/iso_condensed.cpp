/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

/*
 * iso_condensed.cpp
 *
 *  Created on: Sep 19, 2019
 *      Author: fran
 */
#include "cddefines.h"
#include "dense.h"
#include "freebound.h"
#include "helike.h"
#include "helike_cs.h"
#include "hydroeinsta.h"
#include "hydrogenic.h"
#include "hydro_vs_rates.h"
#include "iso.h"
#include "iso_condensed.h"
#include "phycon.h"
#include "thirdparty.h"
#include "thirdparty_quadpack.h"

bool lgpop = true;
bool lgprt = false;

/* This class gets the integrand of the matrix element times the population for
 * calculating the integral given by Brokclehurst MNRAS148 (1970) 417-434
 */
class integrand_matrix_element : public D_fp
{
	public:
	long nelem, ipISO, iplevel, ntop, num_bj, gLo;

	double ELoRyd;
	vector<double> srcterm;
	vector<long> mj;
	bool lgsnk;

	t_iso_sp* sp=&iso_sp[ipISO][nelem];

	multi_arr<double,2,C_TYPE> einji2;

	integrand_matrix_element(long nelem, long ipISO, long iplevel ,long ntop,long num_bj, long gLo,
			double ELoRyd, vector<double> srcterm, vector<long> mj, bool lgsnk, multi_arr<double,2,C_TYPE> einji):
		nelem(nelem), ipISO(ipISO), iplevel(iplevel), ntop(ntop),num_bj(num_bj), gLo(gLo), ELoRyd(ELoRyd),
		srcterm(srcterm), mj(mj),lgsnk(lgsnk), einji2(einji)
	{

	}

	double operator()(double nHi ) const
	{
		double gHi = 2*pow2(nHi);
		double K;
		double IP_Ryd_Hi = 0.;
		if( ipISO == ipH_LIKE )
		{
			IP_Ryd_Hi = hydro_energy(nelem, nHi, 0,0,0)*WAVNRYD;
		}
		else if( ipISO == ipHE_LIKE )
		{
			IP_Ryd_Hi = helike_energy(nelem, nHi, 0,0,0) * WAVNRYD;
		}
		else
		{
			/* Other iso sequences don't exist yet. */
			TotalInsanity();
		}


		//1/pow2(nHi);
		double srct1 = 0., snkt1 =0.;

		/* matrix element */
		K = Kbn( nelem, ipISO, nHi, gHi, iplevel, ntop, IP_Ryd_Hi);

		if(lgsnk)
		{
			for (long j=0 ; j<4 ; j++)
			{
				snkt1 += bfn(j,IP_Ryd_Hi)*einji2[j][num_bj];
				if(lgpop)
					snkt1/=pow2(sp->fb[mj[num_bj]].ipm);
			}


			K *= snkt1;

			/*print integrant values*/
			/*
			fprintf(ioQQQ,"integrant %g %li %g %g %g %g %g %g %g %g %g %g %g\n",
					nHi,iplevel,K,K/snkt1,snkt1,einji2[0][num_bj],einji2[1][num_bj],einji2[2][num_bj],einji2[3][num_bj],
					bfn(0, IP_Ryd_Hi ),bfn(1, IP_Ryd_Hi ),bfn(2, IP_Ryd_Hi ),bfn(3, IP_Ryd_Hi));
			//*/
		}
		else
		{


			for (long j=0; j<4 ; j++)
				srct1 += srcterm[j] * bfn(j, IP_Ryd_Hi );
			srct1 = 1-srct1;
				if(lgpop)
					srct1 += IP_Ryd_Hi/phycon.te_ryd;
			K*=srct1;

			/*print integrant values*/
			/*
			fprintf(ioQQQ,"integrant %g %li %g %g %g %g %g %g %g %g %g %g %g\n %g %g %g %g %g\n",
					nHi,iplevel,K,K/srct1,srct1,srcterm[0],srcterm[1],srcterm[2],srcterm[3],
					bfn(0, IP_Ryd_Hi ),bfn(1, IP_Ryd_Hi ),bfn(2, IP_Ryd_Hi ),bfn(3, IP_Ryd_Hi),
					srcterm[0]*bfn(0, IP_Ryd_Hi ),srcterm[1]*bfn(1, IP_Ryd_Hi ),srcterm[2]*bfn(2, IP_Ryd_Hi ),srcterm[3]*bfn(3, IP_Ryd_Hi ),
					srcterm[0]*bfn(0, IP_Ryd_Hi )+srcterm[1]*bfn(1, IP_Ryd_Hi )+srcterm[2]*bfn(2, IP_Ryd_Hi )+srcterm[3]*bfn(3, IP_Ryd_Hi ));
			//*/
		}

		return K;
	}
};

void add_equilibrium_levels( long nelem, long ipISO, long mmax )
{
	/* This routine adds the collisional-radiative matrix element of the asymptotic high n levels.
	 *
	 * These are not counted in the iso_sequence and are only used to stabilize the condensed matrix
	 * calculation
	 */
	t_iso_sp* sp=&iso_sp[ipISO][nelem];

	double factor = HION_LTE_POP*dense.AtomicWeight[nelem]/
		(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);

	double ConvLTEPOP = powpq(factor,3,2)/(2.*iso_ctrl.stat_ion[ipISO])/phycon.te32;

	double ghi =0.;
	/*This the matrix element */
	double K = 0.;

	/* number of extra discreet levels to be added to the condensed matrix */
	long top = 200;

	/* n of top discreet level */
	long ntop = sp->n_HighestResolved_local + sp->nCollapsed_local + top;

	/* Knot levels for extrapolation using eq. 4.11 of
	 * Brocklehurst MNRAS 148 (1970) 417-434
	 * see appendix E of Gordon book
	 */
	vector<double> srcterm, sink,enerj;
	vector<long> mj;
	multi_arr<double,2,C_TYPE> einji;

	srcterm.resize(4);
	enerj.resize(4);
	einji.alloc(4,4);
	einji.zero();

	for (long j=0; j<4; j++)
		srcterm[j] = 0.;

	get_mj(mmax, mj);

	for (long j=0; j<4; j++)
	{
		long n = sp->st[sp->fb[mj[j]].ipm].n();
		if( ipISO == ipH_LIKE )
		{
			enerj[j] = hydro_energy(nelem, n, 0,0,0)*WAVNRYD;
		}
		else if( ipISO == ipHE_LIKE )
		{
			enerj[j] = helike_energy(nelem, n, 0,0,0) * WAVNRYD;
		}
		else
		{
			/* Other iso sequences don't exist yet. */
			TotalInsanity();
		}

	}
//		enerj[j] = 1./pow2(sp->fb[mj[j]].ipm - sp->n_HighestResolved_local);

	/*Set the minimum level to which the extrapolation will apply */
	//const long MinCondensedLevel=0;

	/* First collapsed level to avoid oscillations */
	const long MinCondensedLevel = sp->numLevels_local - sp->nCollapsed_local;

	for (long level = MinCondensedLevel; level< sp->mmax; level ++ )
	{
		sp->fb[level].SourceBound = 0.;
		for (long j=0;j<4;j++)
			sp->fb[level].SinkBound[j] = 0.;
	}

	/* Obtain source terms */
	SourceTerms(enerj,srcterm,einji);

	/* loop on the extra discreet levels: first summation on the right hand side of
	 * eq. 4.12 of Brocklehurst MNRAS 148 (1970) 417-434
	 */
	for (long Nhi = sp->n_HighestResolved_local + sp->nCollapsed_local + 1;
			Nhi< ntop +1 ; Nhi ++)
	{


		/* energy of upper level in rydbergs
		 * because is Rydberg we consider that there are no differences between iso-sequences
		 */
		double Ehi;
		if( ipISO == ipH_LIKE )
		{
			Ehi = hydro_energy(nelem, Nhi, 0,0,0)*WAVNRYD;
		}
		else if( ipISO == ipHE_LIKE )
		{
			Ehi = helike_energy(nelem, Nhi, 0,0,0) * WAVNRYD;
		}
		else
		{
			/* Other iso sequences don't exist yet. */
			TotalInsanity();
		}
		fprintf(ioQQQ, "Nhi %li Ehi %e\n",Nhi,Ehi);
		//double Ehi = pow2((double)(nelem-ipISO+1)/(double)Nhi);
		ghi = 2.*pow2((double)Nhi);

		double srct1 = 0.;
		vector<double> snkt1(4);

		/* obtain the final source term without the matrix element */
		for (long j=0; j<4 ; j++)
			srct1 += srcterm[j] * bfn(j, Ehi );


		/* Print of the source factors for each Nhi*/
		/*
		 *
		fprintf(ioQQQ,"sfactors %li %g %g %g %g %g %g %g %g %g\n",Nhi,srct1,
				srcterm[0],srcterm[1],srcterm[2],srcterm[3],bfn(0, Ehi ),bfn(1, Ehi ),bfn(2, Ehi ),bfn(3, Ehi ));
		//*/
		srct1 = 1. -srct1;
		if (lgpop)
				srct1 += Ehi/phycon.te_ryd;

		/*if (srct1 < 0 )
			srct1 = 0.;*/
		/* ------ */

		/* sink terms without the matrix element */
		for (long num_bj=0;num_bj<4; num_bj++)
		{

			snkt1[num_bj] = 0.;
			for (long j=0 ; j<4 ; j++)
			{
				snkt1[num_bj] += bfn(j,Ehi)*einji[j][num_bj];
				if(lgpop)
					snkt1[num_bj]/=pow2(sp->fb[mj[num_bj]].ipm);
				/*print factors of each sink pivot factor
				 *
				fprintf(ioQQQ," j %li,Nhi %li num_bj %li, einji[j][num_bj] %g, Ehi %g bfn %g snkt1 %g\n",j,
						Nhi,num_bj, einji[j][num_bj],Ehi, bfn(j, Ehi ),snkt1[num_bj]);
				//*/
			}
		}

		/*print sink factor for extra levels
		 *
		fprintf(ioQQQ,"sinkterms Nhi %li Ehi %e s1 %e, s2 %e, s3 %e, s4 %e\n",Nhi,Ehi,snkt1[0]*pow2(Nhi),snkt1[1]*pow2(Nhi),snkt1[2]*pow2(Nhi),snkt1[3]*pow2(Nhi));
		//*/

		/* calculation of the matrix element for each of the pivotal levels */
		for (long level = MinCondensedLevel; level< mmax; level ++ )
		{

			long iplevel = sp->fb[level].ipm;

			/*if(Nhi % 50 == 0)
				fprintf(ioQQQ,"=====================ipISO %li, nelem %li, mmax %li, level %li,  %li %li\n",ipISO,
					nelem,sp->mmax, iplevel, sp->numLevels_local, sp->numLevels_max);
			*/

			/* calculate the matrix elements */
			double nhi = (double) Nhi;

			/* security check to make sure we are not overlapping the condensed levels */
			if(nhi == sp->st[iplevel].n())
			{
				fprintf(ioQQQ," HERE %li %li %g\n",iplevel, sp->fb[mmax].ipm,nhi);
				cdEXIT(EXIT_FAILURE);
			}

			/*Obtain the matrix elements*/
			K = Kbn( nelem, ipISO, nhi, ghi, iplevel, ntop, Ehi);


			/* sink terms */
			for (long num_bj=0;num_bj<4; num_bj++)
			{

				sp->fb[level].SinkBound[num_bj] += snkt1[num_bj]*K;

				/* second term in the right hand side of eq. 4.12 Brocklehurst 1970 */
				if (Nhi ==  ntop )
					sp->fb[level].SinkBound[num_bj] -=	0.5*snkt1[num_bj]*K;
			}

			/* source terms*/
			sp->fb[level].SourceBound += srct1*K;
			/*
			 *
			 if(sp->st[sp->fb[level].ipm].n() == 400)
				fprintf(ioQQQ,"source level %li level %li mmax %li Src %e %li %e %e\n",sp->st[sp->fb[level].ipm].n(),level,mmax,sp->fb[level].SourceBound, Nhi, K,srct1);
			//*/
			if (Nhi ==  ntop)
			{
				sp->fb[level].SourceBound -= 0.5*srct1*K;
				/*
				 *
				 if(sp->st[sp->fb[level].ipm].n() == 400)
					fprintf(ioQQQ,"ntop source level %li Src %e\n",sp->fb[level].ipm,sp->fb[level].SourceBound);
				//*/
			}
		}

	}

	/* finally calculate the integral term in the right hand side of eq. 4.12 Brocklehurst 1970 */
	for (long level = MinCondensedLevel; level< sp->mmax; level ++ )
	{
		long iplevel = sp->fb[level].ipm;
		double ELoRyd = sp->fb[iplevel].xIsoLevNIonRyd;
		double gLo = sp->st[iplevel].g();

		/* sink terms */
		bool lgsrc = false;
		bool lgsnk = true;
		long num_bj = 0;

		/* calculate the integral from nhi=ntop+1 to infinite */
		/* sink terms first */
		double lbound = (double) ntop+1;
		double epsabs = 0., epsrel = 1.e-5, res, abserr;
		long neval, ier, limit=25;
		//long lenw = 4*limit,
		long inf =1;
		for (long j=0;j<4;j++)
		{
			integrand_matrix_element func(nelem, ipISO, iplevel , ntop,j, gLo, ELoRyd, srcterm,mj,lgsnk,einji);
			/* use quadpack */

			double alist[limit],blist[limit],rlist[limit],elist[limit];
			long iord[limit];

			long last;

			/*
			long iwork[limit];
			double work[lenw];
			dqagi_(func, &lbound, &inf, &epsabs, &epsrel, &res, &abserr, &neval, &ier, &limit,
					&lenw, &last, iwork, work );*/
			dqagie_(func, &lbound, &inf,&epsabs, &epsrel, &limit, &res, &abserr, &neval, &ier,
							 alist, blist, rlist, elist,iord, &last);

			/*print integrated sink terms
			 *
			fprintf(ioQQQ, "sink results iplevel %li j %li %g ier %li\n",sp->st[iplevel].n(),sp->st[sp->fb[mj[j]].ipm].n(),res,ier);
			//*/

			sp->fb[level].SinkBound[j] += res;

		}

		/*now the source terms */
		lgsnk = false;
		lbound = (double) ntop+1;
		epsabs = 0., epsrel = 1.e-5;
		limit=25,
		//lenw = 4*limit,
		inf =1;
		double alist[limit],blist[limit],rlist[limit],elist[limit];
		long iord[limit];

		long last;

		integrand_matrix_element func(nelem, ipISO, iplevel , ntop,num_bj, gLo, ELoRyd,srcterm, mj,lgsrc, einji);
		/* use quadpack */
		/*
		long iwork[limit];
		double work[lenw];
		dqagi_(func, &lbound, &inf, &epsabs, &epsrel, &res, &abserr, &neval, &ier, &limit,
				&lenw, &last, iwork, work );
		*/
		dqagie_(func, &lbound, &inf,&epsabs, &epsrel, &limit, &res, &abserr, &neval, &ier,
						 alist, blist, rlist, elist,iord, &last);

		/*print source integrated terms
		*
		fprintf(ioQQQ, "source results iplevel %li %g %li\n",sp->st[iplevel].n(),res,ier);
		//*/

		sp->fb[level].SourceBound += res;

		sp->fb[level].SourceBound *= ConvLTEPOP*dense.eden*dense.xIonDense[nelem][nelem-ipISO+1];

		/* print final source and sink terms
		*
		fprintf(ioQQQ, "source %g level %li\n",sp->fb[level].SourceBound,sp->st[iplevel].n());
		for(long j=0;j<4;j++)
			fprintf(ioQQQ,"j %li,sink %g, level %li\n",sp->st[sp->fb[mj[j]].ipm].n(),sp->fb[level].SinkBound[j],sp->st[iplevel].n());
		//*/

	}

}

void get_mj(long mmax, vector<long>& mj)
{
	mj.resize(4);

	for (long j=0;j<4; j++)
	{
		if(mmax<4)
		{
			fprintf(ioQQQ," something is wrong number of condensed levels less than 4\n, "
					"mmax %li. I don't have enough condensed levels for the 4 pivotal levels\n",mmax);
			cdEXIT(EXIT_FAILURE);
		}
		else if (mmax < 14)
			mj[j] = mmax -1 - j;
		else
			mj[j] = mmax - 2 - 3*j;
	}
}

void SourceTerms(vector<double> enerj,vector<double>& srcterm,multi_arr<double,2,C_TYPE>& einji)
{
	double eps = 1e-7; //precision of the inverted matrix
	multi_arr<double,2,C_TYPE> eji,ident, oldeji;
	eji.alloc(4,4);
	eji.zero();
	oldeji.alloc(4,4);
	oldeji.zero();
	einji.alloc(4,4);
	einji.zero();
	ident.alloc(4,4);
	ident.zero();

	for (long int j=0; j<4 ; j++)
	{

		for (long int i =0; i<4 ; i++)
		{
			eji[j][i] = bfn(i, enerj[j] );
			oldeji[j][i] =eji[j][i];
		}
	}
	/// LAPACK INVERSION OF EJI HERE
		/* this will solve for the unknown Aj values in eq 4.11 of
		 * Brocklehurst MNRAS 148 (1970) 417-434
		 */
		static vector<int32> ipiv;
		int32 nerror = 0;
		ipiv.resize(4);

		getrf_wrapper(4,4,eji.data(),4,&ipiv[0],&nerror);

		/* we solve the equation system for the I matrix to get the inverse of eji*/
		for ( long int j=0; j< 4; j++)
			einji[j][j] = 1.;

		getrs_wrapper('N',4,4,eji.data(),4,&ipiv[0],einji.data(),4,&nerror);


		if (nerror < 0)
		{
			fprintf(ioQQQ,"SourceTerms: Error in lapack inversion. Element %i has an illegal value\n",-nerror);
			cdEXIT(EXIT_FAILURE);
		}
		/* now we use the inverted matrix to obtain he source and sink terms */

		// This is E-1 x 1 (1 being a 4 elements vector with every element =1)
		for (long i=0 ; i<4; i++)
		{
			srcterm[i] = 0.;
			for (long j =0 ; j<4 ; j++)
			{
				if(!lgpop)
					srcterm[i] += einji[i][j];
				else
					srcterm[i] += einji[i][j]*(1+enerj[j]/phycon.te_ryd);

				/* testing the inversion */

				for (long k = 0; k<4; k++)
					ident[i][j] += oldeji[i][k]*einji[k][j];

				double val = ident[i][j] -1;
				bool lgerror=false;
				if (i==j && abs(val) >eps  )
					lgerror = true;
				else if (i!= j && ident[i][j] > eps)
				{
					val = ident[i][j];
					lgerror = true;
				}

				if(lgerror)
				{
					fprintf(ioQQQ,"SourceTerms: error over the precision (%e) in the inverted matrix in i %li and j %li: %e\n",eps,i,j,val);
					cdEXIT(EXIT_FAILURE);
				}

			}
		}


}

double bfn(long j, double En)
{
	/* this routine accounts for the different terms of the asymptotic form of the
	 * departure coefficients chosen by Brocklehurst MNRAS 148 (1970) 417-434 eq. (4.11)
	 * NOTE: only the second term of the right hand side of the formula (4.11) is
	 * coded here as the 1 will be subtracted in add_equilibrium_levels()
	 *
	 * NOTE 2: In formula 4.11 j is starting at 1 while here is starting in 0
	 *
	 */
	if (En <= 0)
	{
		fprintf(ioQQQ,"bfn: Invalid value of En %e\n",En);
		cdEXIT(EXIT_FAILURE);
	}

	double bs = 0.;
	long mu = 5;
	/* this is like that in brocklehurst code Brocklehurst & Salem Comp. Phys. Comm. 13 (1977) 39
	 * a better formula not discontinuous in Te would be preferable
	 */
	if (phycon.te <= 1000)
		mu = 3;
	bs = powpq(En,mu+2*j,2)/log(En);

	if (lgpop)
			bs *=(1+En/phycon.te_ryd);

	/* test this routine printing values */
	if(lgprt)
	{
		fprintf(ioQQQ,"test bfn\n========\n");
		for (long n=10;n<1001;n++)
		{
			double Etest=hydro_energy(0, n, 0,0,0)*WAVNRYD;
			vector<double> btest;
			btest.resize(4);
			for (long jtest=0;jtest<4;jtest++)
			{
				btest[jtest] = powpq(Etest,mu+2*jtest,2)/log(Etest);
				if (lgpop )
					btest[jtest] *=(1+Etest/phycon.te_ryd);
			}


			fprintf(ioQQQ,"%li %e %e %e %e %e\n",n,btest[0],btest[1],btest[2],btest[3], Etest );
		}
		lgprt =false;
	}


	return bs;
}

double Kbn( long nelem, long ipISO, double Nhi, double gHi, long iplevel, long ntop, double IP_Ryd_Hi)
{
	/* This routine calculates matrix element for not integer quantum
	 * numbers to integers. To do that it is needed to use or reprogram h-like and
	 * he-like collisions to accept non-integers */
	t_iso_sp* sp=&iso_sp[ipISO][nelem];
	ColliderDensities colld(colliders);

	double K =0.,cs,rateCoef=0.;
	double Aul=0.;

	/*double factor = HION_LTE_POP*dense.AtomicWeight[nelem]/
		(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);
	double ConvLTEPOP = powpq(factor,3,2)/2./phycon.te32;

	double ltepop = dense.eden*dense.xIonDense[nelem][nelem-ipISO+1]*ConvLTEPOP*exp(IP_Ryd_Hi*TE1RYD/phycon.te);
	*/

	long nLo = sp->st[iplevel].n();
	long lLo = sp->st[iplevel].l();
	long sLo = sp->st[iplevel].S();
	long jLo = sp->st[iplevel].j();
	long gLo = sp->st[iplevel].g();
	double IP_Ryd_Lo = sp->fb[iplevel].xIsoLevNIonRyd;
	double tauLo = sp->st[iplevel].lifetime();
	//double DeltaE_Ryd = IP_Ryd_Hi - IP_Ryd_Lo;

	double EnerWN = RYD_INF* abs(IP_Ryd_Hi - IP_Ryd_Lo);
	double EnerErg = ERG1CM*EnerWN;
	/*
	double electronrate;
	*/
	/* Here the collisional rate coefficients for discreet n are obtained */
	if (Nhi <= ntop)
	{
		if (iplevel < sp->numLevels_local - sp->nCollapsed_local)
			Aul = hydro_transprob_collapsed_to_resolved( nelem, (long) Nhi, nLo, lLo );
		else if( iplevel >= sp->numLevels_local - sp->nCollapsed_local)
			Aul = hydro_transprob_collapsed_to_collapsed( nelem, (long) Nhi, nLo );


		/* obtain collision strength for all possible colliders*/
		for (long ipCollider = ipELECTRON ; ipCollider <= ipALPHA ; ipCollider++)
		{
		//long ipCollider = 0;
			cs = 0.;
			const char *where = "      ";

				if( ipISO == ipH_LIKE )
				{
					cs += GetHlikeCollisionStrength( nelem, ipCollider,
							(long) Nhi, -1, -1, gHi, IP_Ryd_Hi,
							nLo, lLo, sLo, gLo, IP_Ryd_Lo,
							Aul, tauLo, EnerWN, EnerErg, &where );
				}
				else if( ipISO == ipHE_LIKE )
				{
					cs += GetHelikeCollisionStrength( nelem, ipCollider,
							(long) Nhi, -1, -1, -1, gHi, IP_Ryd_Hi,
							nLo, lLo, sLo, jLo, gLo, IP_Ryd_Lo,
							Aul, tauLo, EnerWN, EnerErg, &where );
				}

			/* convert collision strength Y to rate coefficients */
			rateCoef = cs * ratefact(nelem,ipCollider)/gHi;

			/* fill the matrix element */
			K += rateCoef* colld.density(ipCollider);

			/*print rates from collider */
			/*if ((long)Nhi % 50 ==0 )
				fprintf(ioQQQ,"rates in K  nhi %g nlo %li ratecoef %g ipCollider %li dens %g\n",
				Nhi,nLo,rateCoef, ipCollider, colld.density(ipCollider));
			 */
			/*keep electron rates*/
			/*
			//double electronrate;
			if( ipCollider == ipELECTRON)
			 	 electronrate = ratecoef;
			*/
		}

		/*as the matrix elements are the contribution to the low levels radiative terms should be there*/
		K +=Aul;

		/* multiply by LTE populations as we will be working with departures */
		//K *=ltepop*gHi;//*exp(-IP_Ryd_Hi*TE1RYD/phycon.te);

	}
	/* in case that the Nhi is continuous double variable (for integration)*/
	else
	{

		Aul = HydroEinstA(Nhi,(double)nLo);

		/* in case the low state is resolved we split the collapsed Aul in each of the lLo assuming that only
		 * lLo-1 -> lLo and lLo+1-> lLo in the same spin system are relevant  */
		if (lLo >= 0)
			Aul /= (2.*nLo*nLo);

		cs =0.;
			cs = hydro_Lebedev_deexcit(nelem, ipH_LIKE, Nhi, nLo, gLo, IP_Ryd_Lo);
			//hydro_vs_deexcit( Nhi, gHi, IP_Ryd_Hi, nLo, gLo, IP_Ryd_Lo, Aul );


		/*effective rate coefficient */
		rateCoef = cs*ratefact(nelem,ipELECTRON)/gHi;
		/*matrix element multiplied by the LTE populations as we are working with b's*/
		K = (rateCoef*dense.eden+Aul);//*ltepop*gHi;//*exp(-IP_Ryd_Hi*TE1RYD/phycon.te);

		/*
		 electronrate=rateCoef;
		 */

	}
	if (lgpop)
		K *= pow2(Nhi);
//	K*= exp(DeltaE_Ryd*TE1RYD/phycon.te)*pow2(Nhi/(double)nLo);
/* test for comparison
 * PRINT THE K() TO CHECK WHAT COLLISIONS ARE IMPORTANT IN WHAT RANGES*/
/*
	if (Nhi > 600)
		fprintf(ioQQQ," in k() K(%g,%li) = %g, Aul %g %g %g %g %g again K %g \n",Nhi,nLo,K,Aul, rateCoef, ltepop, IP_Ryd_Hi,gHi,
				(electronrate*dense.eden+Aul)*ltepop*gHi*exp(DeltaE_Ryd*TE1RYD/phycon.te)*pow2(Nhi/(double)nLo));
	//*/
	return K;
}

