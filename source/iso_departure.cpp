/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_departure_solver solve for iso-sequence departure coefficients */
#include "cddefines.h"
#include "atmdat.h"
#include "continuum.h"
#include "conv.h"
#include "dense.h"
#include "dynamics.h"
#include "elementnames.h"
#include "freebound.h"
#include "grainvar.h"
#include "iso.h"
#include "ionbal.h"
#include "iterations.h"
#include "iso_condensed.h"
#include "mole.h"
#include "opacity.h"
#include "phycon.h"
#include "prt.h"
#include "rfield.h"
#include "thirdparty.h"
#include "trace.h"
#include "two_photon.h"
#include "vectorize.h"
#include "prt.h"

void get_condensed_matrix( long ipISO, long m[][NISO],vector < vector <double> >& L, multi_arr<double,2,C_TYPE> z, multi_arr<double,2,C_TYPE> &c,
		long mmax, long NumMaxLevels, long deg);
void get_interpolated_vector(long ipISO, vector <double>& xexpan, vector <double>& xcon,long m[][NISO],vector < vector <double> >& L,
		long mmax, long NumMaxLevels, long deg);
void set_L(long m[][NISO], vector <vector <double>>& L, long mmax, long ipISO, long NumMaxLevels, long deg);
int unitary_test_U(long ipISO, long NumMaxLevels, long mmax, long m[][NISO], long U);

multi_arr<double,2,C_TYPE> z, SaveZ,c;

/*solve for departure coefficients*/
void iso_departure_solver( const long int ipISO, const long int nelem, double &renorm, bool lgPrtMatrix )
{

	DEBUG_ENTRY( "iso_departure_solver()" );

	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	const long int numlevels_local = sp->numLevels_local;

	bool lgNegPop=false;
	long int level,
		ipHi,
		ipLo;

	/* this block of variables will be obtained and freed here */
	double source=0., sink=0.;

	/* obtain here the number and the index of pivotal levels */

	long m[sp->numLevels_max][NISO];
	long mmax= sp->mmax;

	if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
	{
		for (long ipHi =0 ; ipHi<mmax; ipHi++)
			m[ipHi][ipISO] =sp->fb[ipHi].ipm;
	}

	/* check that we were called with valid charge */
	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );

	/* these two collision rates must be the same or we are in big trouble,
	 * since used interchangeably */
	ASSERT( ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0]< SMALLFLOAT ||
		fabs( (sp->fb[0].ColIoniz* dense.EdenHCorr) /
		SDIV(ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0] ) - 1.) < 0.001 );

	if( (dense.IonHigh[nelem] >= nelem - ipISO) &&
		 (dense.IonLow[nelem] <= nelem - ipISO) &&
		 ionbal.RateRecomIso[nelem][ipISO] > 0. )
	{
		/* get simple estimate of atom to ion ratio */
		sp->xIonSimple = ionbal.RateIonizTot(nelem,nelem-ipISO)/ionbal.RateRecomIso[nelem][ipISO];
	}
	else
	{
		sp->xIonSimple = 0.;
	}

	/* which case atom to solve??? */
	if( dense.xIonDense[nelem][nelem+1-ipISO] < SMALLFLOAT )
	{
		/* don't bother if no ionizing radiation */
		strcpy( sp->chTypeAtomUsed, "zero " );
		if( trace.lgTrace && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "     iso_departure_solver iso=%2ld nelem=%2ld simple II/I=%10.2e so not doing equilibrium, doing %s.\n",
				ipISO, nelem, sp->xIonSimple, sp->chTypeAtomUsed );
		}

		/* total ionization will just be the ground state */
		lgNegPop = false;
		sp->st[0].Pop() = dense.xIonDense[nelem][nelem-ipISO];
		for( long n=1; n < numlevels_local; n++ )
		{
			sp->st[n].Pop() =  0.;
		}
		sp->qTot2S = 0.;
	}
	else
	{
		static vector<double> creation,Save_creation,CondCreation, error;

		creation.resize(numlevels_local);

		ASSERT( dense.xIonDense[nelem][nelem+1-ipISO] >= 0.f );

		/* Bolztmann and ionization potential factors */
		avx_ptr<double> arg(1,numlevels_local), bstep(1,numlevels_local);

		qList::const_iterator StElm = sp->st.begin();
		static vector<double> Boverg, IP_exp_overg;
		Boverg.resize(numlevels_local);
		IP_exp_overg.resize(numlevels_local);

		/*Boltzmann factors and energy exponentials */
		for (unsigned lev = 0; lev < Boverg.size(); ++lev)
		{
			Boverg[lev] = 1.0/(double)StElm[lev].g();
			IP_exp_overg[lev] = sp->st[lev].ConBoltz()/(double)StElm[lev].g();
		}

		for( level=0; level < numlevels_local; level++ )
		{
			if (level >0)
			{
				fixit("Wouldn't need to mask this out if levels were in order");
				/*WARNING: this argument id negative so the equations are adapted to
				 * it being <0
				 * */
				arg[level] = -max((StElm[level].energy().K()-StElm[level-1].energy().K()) / phycon.te, -38.);
			}
		}
		vexp( arg.ptr0(), bstep.ptr0(), 1, numlevels_local );

		double factor = HION_LTE_POP*dense.AtomicWeight[nelem]/
			(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);

		double ConvLTEPOP = powpq(factor,3,2)/(2.*iso_ctrl.stat_ion[ipISO])/phycon.te32;
		for( level=0; level < numlevels_local; level++ )
		{
			/* total recombination from once more ionized [cm-3 s-1]
			 * RateCont2Level units are s^-1 and contains electron density */
			creation[level] = sp->fb[level].RateCont2Level/ConvLTEPOP*IP_exp_overg[level]/dense.eden;
		}
		
		double ctsource=0.0, ctsink=0.0, ctrec=0.0;
		/* now charge transfer - all into/from ground, two cases, H and not H
		 * This is a source/sink term as this is coming from other ions
		 * (even if they are Hydrogen ions)
		 * */
		if( nelem==ipHYDROGEN )
		{
			/* charge transfer, hydrogen onto everything else */
			/* charge exchange recombination
			 * units s-1 including ions density*/
			ctrec += atmdat.CharExcRecTotal[nelem]/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden;
			ctsink += atmdat.CharExcIonTotal[nelem];
		}
		else if( nelem==ipHELIUM && ipISO==ipHE_LIKE )
		{
			/* this is recom of He due to ct with all other gas constituents */
			ctrec += atmdat.CharExcRecTotal[nelem]/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden;
			ctsink += atmdat.CharExcIonTotal[nelem];
		}
		else
		{
			for (long nelem1=ipHYDROGEN; nelem1 < t_atmdat::NCX; ++nelem1)
			{
				long ipISO=nelem1;
				ctrec += atmdat.CharExcRecTo[nelem1][nelem][nelem-ipISO]*iso_sp[ipISO][nelem1].st[0].Pop()/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden;
				ctsink += atmdat.CharExcIonOf[nelem1][nelem][nelem-ipISO]*dense.xIonDense[nelem1][1];
			}					
		}
		ctsource += ctrec;
		
		if ( nelem > ipISO )
		{
			double ction=0.0;
			if( nelem==ipHELIUM )
			{
				ctsink += atmdat.CharExcRecTotal[nelem];
				ction += atmdat.CharExcIonTotal[nelem]/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden;;
			}
			else
			{
				for (long nelem1=ipHYDROGEN; nelem1 < t_atmdat::NCX; ++nelem1)
				{
					long ipISO1=nelem1;
					ctsink += atmdat.CharExcRecTo[nelem1][nelem][nelem-(ipISO+1)]*iso_sp[ipISO1][nelem1].st[0].Pop();
					ction += atmdat.CharExcIonOf[nelem1][nelem][nelem-(ipISO+1)]*dense.xIonDense[nelem1][1]/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden;
				}
			}
			ctsource += ction;
		}

		/* now do the 2D array to fill the CR matrix*/

		z.alloc(numlevels_local,numlevels_local);
		z.zero();

		static vector<double> coll_down,RadDecay,pump;
		static vector<int32> ipiv;
		vector < vector <double> > L(sp->numLevels_max,vector<double> (sp->numLevels_max,0.));

		coll_down.resize(numlevels_local);
		RadDecay.resize(numlevels_local);
		pump.resize(numlevels_local);

		ColliderDensities colld(colliders);

		enum { DEBUG_RATES = false };

		if( DEBUG_RATES && iterations.lgLastIt )
			fprintf( stdout, "# ipISO\tnelem\tlevel\tcollDown\tcollIonz\tradRecom\tRadDecay\n" );

		/* this branch is main solver, full level populations
		 * assert since this code must change if NISO ever increased */
		ASSERT( NISO == 2 );

		strcpy( sp->chTypeAtomUsed, "Depart" );
		if( trace.lgTrace && (nelem == trace.ipIsoTrace[ipISO]) && iso_ctrl.lgMatCondOff[ipISO][nelem])
		{
			fprintf( ioQQQ, "     iso_departure_solver iso=%2ld nelem=%2ld doing regular matrix inversion, %s\n",
				ipISO, nelem, sp->chTypeAtomUsed );
		}

		/* master balance equation, use when significant population */
		for( level=0; level < numlevels_local; level++ )
		{
			double coll_down_total = 0.;
			if (level != 0)
			{
				for ( ipLo = 0; ipLo < level; ++ipLo )
					Boverg[ipLo] *= bstep[level];
			}

			/* all process depopulating level and placing into the continuum
			 * this does NOT include grain charge transfer ionization, added below */
			z[level][level] += sp->fb[level].RateLevel2Cont;

			/* all processes populating lower levels */
			for( ipLo=0; ipLo < level; ipLo++ )
			{
				coll_down[ipLo] = sp->trans(level,ipLo).Coll().ColUL( colld );
				if( DEBUG_RATES )
				{
					coll_down_total += coll_down[ipLo];
				}
				if ( rfield.lgPlasNu && sp->trans(level,ipLo).EnergyRyd()<rfield.plsfrq  )
				{
					RadDecay[ipLo] = iso_ctrl.SmallA;
					pump[ipLo] = iso_ctrl.SmallA;
				}
				else
				{
					RadDecay[ipLo] = MAX2( iso_ctrl.SmallA, sp->trans(level,ipLo).Emis().Aul()*
												  (sp->trans(level,ipLo).Emis().Ploss()) );
					/*pump is the UPWARD pumping rate */
					pump[ipLo] = MAX2( iso_ctrl.SmallA, sp->trans(level,ipLo).Emis().pump() );
				}
			}

			/*error generator*/
			if( iso_ctrl.lgRandErrGen[ipISO] )
			{
				for( ipLo=0; ipLo < level; ipLo++ )
				{
					coll_down[ipLo] *= sp->ex[level][ipLo].ErrorFactor[IPCOLLIS];
					RadDecay[ipLo] *= sp->ex[level][ipLo].ErrorFactor[IPRAD];
					pump[ipLo] *= sp->ex[level][ipLo].ErrorFactor[IPRAD];
				}
			}

			double glev = (double)StElm[level].g(), rglev = 1.0/glev;


			for( ipLo=0; ipLo < level; ipLo++ )
			{
				double coll_up = coll_down[ipLo] * glev *
					Boverg[ipLo];

				z[ipLo][ipLo] += coll_up + pump[ipLo] ;
				z[ipLo][level] -= coll_up + pump[ipLo] ;

				z[ipLo][level] *= rglev/Boverg[ipLo];

				double pump_down = pump[ipLo] *
					(double)StElm[ipLo].g() * rglev;


				z[level][level] += RadDecay[ipLo] + pump_down + coll_down[ipLo];
				z[level][ipLo] -= RadDecay[ipLo] + pump_down + coll_down[ipLo];

				z[level][ipLo] *=Boverg[ipLo]*glev;

				/* find total collisions out of 2^3S to singlets. */
				if( (level == 1) && (ipLo==0) )
				{
					sp->qTot2S += coll_down[ipLo]/dense.EdenHCorr;
				}
				if( (ipLo == 1) && (ipISO==ipH_LIKE || StElm[level].S()==1) )
				{
					sp->qTot2S += coll_up/dense.EdenHCorr;
				}
			}

			if( DEBUG_RATES && iterations.lgLastIt )
			{
				fprintf( stdout,
					"%2ld\t%2ld\t%2ld\t"
					"%.4e\t"
					"%.4e\t"
					"%.4e\t"
					"%.4e\n",
					ipISO, nelem, level,
					coll_down_total,
					iso_sp[ipISO][nelem].fb[level].ColIoniz * dense.eden,
					iso_sp[ipISO][nelem].fb[level].RadRecomb[0] * dense.eden,
					1. / iso_sp[ipISO][nelem].st[level].lifetime() );
			}
		}

		if( DEBUG_RATES && iterations.lgLastIt )
			fprintf( stdout, "\n\n" );

		ASSERT( ipISO <= ipHE_LIKE );
		for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
		{
			// induced two photon emission - special because upward and downward are
			// not related by ratio of statistical weights
			// iso.lgInd2nu_On is controlled with SET IND2 ON/OFF command

			fixit("need Pesc or the equivalent to multiply AulTotal?");
			// downward rate

			z[tnu->ipHi][tnu->ipLo] -= tnu->AulTotal;
			z[tnu->ipHi][tnu->ipHi] += tnu->AulTotal;

			z[tnu->ipHi][tnu->ipLo] *= Boverg[tnu->ipLo]*(double)StElm[tnu->ipHi].g();

		}

		if (iso_ctrl.lgInd2nu_On)
		{
			for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
			{
				z[tnu->ipHi][tnu->ipLo] -= tnu->induc_dn*Boverg[tnu->ipLo]*(double)StElm[tnu->ipHi].g();
;
				z[tnu->ipHi][tnu->ipHi] += tnu->induc_dn;

				// upward rate
				z[tnu->ipLo][tnu->ipHi] -= tnu->induc_up/Boverg[tnu->ipLo]/(double)StElm[tnu->ipHi].g();
				z[tnu->ipLo][tnu->ipLo] += tnu->induc_up;

			}
		}

		/*NOW  GRAINS */
		/* grain charge transfer recombination and ionization to ALL other stages */
		for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
		{
			if( ion!=nelem-ipISO )
			{
				source += gv.GrainChTrRate[nelem][ion][nelem-ipISO]*dense.xIonDense[nelem][ion]/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden/dense.xIonDense[nelem][nelem-ipISO+1];

				sink += gv.GrainChTrRate[nelem][nelem-ipISO][ion];

				source	+= mole.xMoleChTrRate[nelem][ion][nelem-ipISO]* dense.xIonDense[nelem][ion]* atmdat.lgCTOn/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden/dense.xIonDense[nelem][nelem-ipISO+1];

				sink += mole.xMoleChTrRate[nelem][nelem-ipISO][ion] * atmdat.lgCTOn;

			}
		}

		source += mole.source[nelem][nelem-ipISO]/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden/dense.xIonDense[nelem][nelem-ipISO+1];
		sink += mole.sink[nelem][nelem-ipISO];

		/* add in advection - these terms normally zero */
		if( iteration > dynamics.n_initial_relax+1 &&
			( dynamics.lgAdvection || dynamics.lgTimeDependentStatic ) &&
			dynamics.Rate != 0.0 &&
			!dynamics.lgEquilibrium && dynamics.lgISO[ipISO] )
		{
			source += dynamics.Source[nelem][nelem-ipISO]/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden/dense.xIonDense[nelem][nelem-ipISO+1];
			/* >>chng 02 Sep 06 rjrw -- advective term not recombination */
			sink += dynamics.Rate;
		}

		/* ionization from/recombination to next lower ionization stages */
		for(long ion_from=dense.IonLow[nelem]; ion_from < MIN2( dense.IonHigh[nelem], nelem-ipISO ) ; ion_from++ )
		{
			if( ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] >= 0. )
			{

				/* ionization from lower ionization stages, cm-3 s-1 */
				source += ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] * dense.xIonDense[nelem][ion_from]/ConvLTEPOP*2.*IP_exp_overg[0]/dense.eden/dense.xIonDense[nelem][nelem-ipISO+1];
			}

			/* recombination to next lower ionization stage, s-1 */
			if( ion_from == nelem-1-ipISO )
				sink += ionbal.RateRecomTot[nelem][ion_from];

		}

		ASSERT( source >= 0.f );
		if (0)
		{
			/*
			 * Collisional ionization and photoionization can only be to the ground state in H iso-sequences
			 * to conserve energy.
			 */
			creation[0] += source;
				/*recombination is only done from the ground state */
			z[0][0] += sink;
		}
		else
		{
			// Try Boltzmann weighting to capture LTE limit correctly
			t_iso_sp* sp = &iso_sp[ipISO][nelem];
			double partfun=0.0;
			for ( level = 0; level < numlevels_local; level++ )
			{
				partfun += sp->st[level].Boltzmann()*sp->st[level].g();
			}
			source /= partfun;
			for( level=0; level < numlevels_local; level++ )
			{

				creation[level] += source*
						sp->st[level].Boltzmann()*sp->st[level].g();
				z[level][level] += sink;
			}
		}
		creation[0] += ctsource;
		z[0][0] += ctsink;
		
		if( sp->trans(iso_ctrl.nLyaLevel[ipISO],0).Coll().rate_lu_nontherm() * iso_ctrl.lgColl_excite[ipISO] > 0. )
		{
			/* now add on supra thermal excitation */
			for( level=1; level < numlevels_local; level++ )
			{
				double RateUp , RateDown;

				RateUp = sp->trans(level,0).Coll().rate_lu_nontherm();
				RateDown = RateUp * (double)sp->st[ipH1s].g() /
					(double)sp->st[level].g();

				/* total rate out of lower level */
				z[ipH1s][ipH1s] += RateUp;

				/* rate from the upper level to ground */
				z[level][ipH1s] -= RateDown*Boverg[ipH1s]*(double)sp->st[level].g();

				/* rate from ground to upper level */
				z[ipH1s][level] -= RateUp/Boverg[ipH1s]/(double)sp->st[level].g();;

				z[level][level] += RateDown;
			}
		}
		/* ===================================================================
		 *
		 * at this point all matrix elements have been established
		 *
		 * ==================================================================== */
		long NumMaxLevels = numlevels_local;
		long deg =6;



		/* *******************************************************************
		 * Matrix condensation                                               *
		 * *******************************************************************  */
		if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
		{
			/*Unitary tests*/
			if (1)
			{
				for (long U=1;U<3;U++)
				{
					int Test;
					if ((Test = unitary_test_U(ipISO, NumMaxLevels, mmax, m, U))==0)
						cdEXIT(EXIT_FAILURE);
				}
			}

			set_L(m,L ,mmax, ipISO, NumMaxLevels,deg);

			c.alloc(mmax,mmax);
			c.zero();

			get_condensed_matrix(ipISO, m, L,z,c,mmax, NumMaxLevels, deg );

			SaveZ = c;

		}
		else
		{
			/* save matrix, this allocates SaveZ */
			SaveZ = z;
		}
		/********************************************************************** */

		Save_creation.resize(mmax);

		if(iso_ctrl.lgMatCondOff[ipISO][nelem] )
		{
			Save_creation.resize(NumMaxLevels);

			for( ipLo=0; ipLo < NumMaxLevels; ipLo++ )
				Save_creation[ipLo] = creation[ipLo];

			if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO])))
			{
				const long numlevels_print = numlevels_local;
				fprintf( ioQQQ, "  pop level     others => (iso_departure_solver)\n" );
				for( long n=0; n < numlevels_print; n++ )
				{
					fprintf( ioQQQ, "  %s %s %2ld", iso_ctrl.chISO[ipISO], elementnames.chElementNameShort[nelem], n );
					for( long j=0; j < numlevels_print; j++ )
					{
						fprintf( ioQQQ,"\t%.9e", z[j][n] );
					}
					fprintf( ioQQQ, "\n" );
				}
			}
		}

		if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO])))
		{
			for( long n=0; n < NumMaxLevels; n++ )
			{
				fprintf( ioQQQ,"\t%.9e", creation[n] );
			}

			fprintf( ioQQQ, "\n" );
		}

		if( lgPrtMatrix )
		{
			valarray<double> cr( get_ptr(creation), creation.size() );
			prt.matrix.prtRates( numlevels_local, z, cr );
		}

		int32 nerror;

		nerror = 0;
		ipiv.resize(numlevels_local);

		if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
		{
			CondCreation.resize(mmax);
			ipiv.resize(mmax);

			for( long j=0; j<mmax; j++)
				CondCreation[j]=creation[m[j][ipISO]];


			/*extrapolation source and sink terms */
			if (0)
			{
				vector<long> mj;

				get_mj(mmax, mj);

				for (long j=0; j<mmax ; j++)
				{
					CondCreation[j] += sp->fb[j].SourceBound;

						for (long i =0; i<4 ; i++)
					{

						c[j][mj[i]] -= sp->fb[j].SinkBound[i];
					}
				}
			}

			for( ipLo=0; ipLo < mmax; ipLo++ )
				Save_creation[ipLo] = CondCreation[ipLo];

			if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO])))
			{
				fprintf(ioQQQ,"mmax %li\n",mmax);

				for( long n=0; n < mmax; n++ )
				{
					fprintf( ioQQQ,"\t%.9e", CondCreation[n] );
				}
			}

			/**** INVERT THE CONDENSED MATRIX ****/
			getrf_wrapper(mmax,mmax,c.data(),mmax,&ipiv[0],&nerror);

			getrs_wrapper('N',mmax,1,c.data(),mmax,&ipiv[0],&CondCreation[0],mmax,&nerror);

		}
		else
		{

			/**** INVERT THE MATRIX ****/
			getrf_wrapper(numlevels_local,numlevels_local,
					z.data(),numlevels_local,&ipiv[0],&nerror);

			getrs_wrapper('N',numlevels_local,1,z.data(),numlevels_local,&ipiv[0],
					&creation[0],numlevels_local,&nerror);
		}

		if( nerror != 0 )
		{
			fprintf( ioQQQ, " iso_departure_solver: dgetrs finds singular or ill-conditioned matrix\n" );
			cdEXIT(EXIT_FAILURE);
		}
		//return to extended matrix
		if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
		{
			get_interpolated_vector(ipISO, creation,CondCreation,m,L,mmax, NumMaxLevels, deg );
		}

		long InvertMatrixLevels = NumMaxLevels;
		error.resize(numlevels_local);

		if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
		{
			InvertMatrixLevels = mmax;
			error.resize(mmax);
		}

		for( level=ipH1s; level < InvertMatrixLevels; level++ )
		{
			double qn = 0., qx = 0.;
			error[level] = 0.;
			for( long n=ipH1s; n < InvertMatrixLevels; n++ )
			{
				double q= SaveZ[n][level];

				if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
					q *=CondCreation[n];
				else
					q*=creation[n];

				/* remember the largest size of element in sum to div by below */
				if ( q > qx )
					qx = q;
				else if (q < qn)
					qn = q;

				error[level] += q;
			}

			if (-qn > qx)
				qx = -qn;

			if( qx > 0. )
			{
				error[level] = (error[level] - Save_creation[level])/qx;
			}
			else
			{
				error[level] = 0.;
			}
		}

		/* remember largest residual in matrix inversion */
		double BigError = -1.;
		long int level_error = -1;
		/* >>chng 06 aug 28, from numLevels_max to _local. */
		for( level=ipH1s; level < InvertMatrixLevels ; level++ )
		{
			double abserror;
			abserror = fabs( error[level]);
			/* this will be the largest residual in the matrix inversion */
			if( abserror > BigError )
			{
				BigError = abserror;
				level_error = level;
			}
		}
		if( BigError > FLT_EPSILON )
		{
			if( !conv.lgSearch )
				fprintf(ioQQQ," PROBLEM" );

			fprintf(ioQQQ,
				" iso_departure_solver zone %.2f - largest residual in iso=%li %s nelem=%li matrix inversion is %g "
				"level was %li Search?%c \n",
				fnzone,
				ipISO,
				elementnames.chElementName[nelem],
				nelem ,
				BigError ,
				level_error,
				TorF(conv.lgSearch) );
		}

		/*NOW EXTRA STUFF FOR LTE AND CHECKS */

		// Force level balance to LTE
		if ( iso_ctrl.lgLTE_levels[ipISO] )
		{
			for ( level = 0; level < numlevels_local; level++ )
			{
				creation[level] = 1.;
			}
		}

		for( level = numlevels_local-1; level > 0; --level )
		{
			/* check for negative populations */
			if( creation[level] < 0. )
				lgNegPop = true;
		}

		if( lgNegPop && dense.lgSetIoniz[nelem] )
		{
			//  simulation can become unphysical if ionization is fixed.
			//  in this case, just put everything in ground.
			//  It's really the best we can do.
			for( level = 1; level < numlevels_local; ++level )
				creation[level] = 0.;
			sp->st[0].Pop() = dense.xIonDense[nelem][nelem-ipISO];
			// flip flag back
			lgNegPop = false;
		}

		/*Fill departure coefficients */
		for( long level=0; level < numlevels_local; level++ )
		{
			sp->st[level].DepartCoef()=creation[level];
			fprintf(ioQQQ,"level %li dep %e\n",level,creation[level]);
		}
		for( long level=iso_sp[ipISO][nelem].numLevels_local; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
			sp->st[level].DepartCoef() = 1.;


		/* put level populations into master array */
		for( level=0; level < numlevels_local; level++ )
		{
			ASSERT( creation[level] >= 0. );
			sp->st[level].Pop() = sp->fb[level].PopLTE*dense.eden*dense.xIonDense[nelem][nelem-ipISO+1]*creation[level];

			if( sp->st[level].Pop() <= 0 && !conv.lgSearch )
			{
				fprintf(ioQQQ,
					"PROBLEM non-positive level pop for iso = %li, nelem = "
					"%li = %s, level=%li val=%.3e, dep=%.3e, nTotalIoniz %li\n",
					ipISO,
					nelem ,
					elementnames.chElementSym[nelem],
					level,
					sp->st[level].Pop() ,
					creation[level],
					conv.nTotalIoniz);
			}
		}

		/* zero populations of unused levels. */
		for( level=numlevels_local; level < sp->numLevels_max; level++ )
		{
			sp->st[level].Pop() = 0.;
			/* >>chng 06 jul 25, no need to zero this out, fix limit to 3-body heating elsewhere. */
			/* sp->st[level].PopLTE = 0.; */
		}

		/* TotalPopExcited is sum of excited level pops */
		/* renormalize the populations to agree with ion solver
		 * is this even needed? The populations are automatically renormalized when
		 * converted from departure coefficients by the ionized ion density
		 * */
		iso_renorm( nelem, ipISO, renorm );

		double TotalPopExcited = 0.;
		/* create sum of populations */
		for( level=1; level < numlevels_local; level++ )
			TotalPopExcited += sp->st[level].Pop();
		ASSERT( TotalPopExcited >= 0. );
		double TotalPop = TotalPopExcited + sp->st[0].Pop();

		/* option to force ionization */
		if( dense.lgSetIoniz[nelem] )
		{
			if( !fp_equal( TotalPop, (double)dense.xIonDense[nelem][nelem-ipISO] ) )
			{
				if( TotalPopExcited >= dense.xIonDense[nelem][nelem-ipISO] )
				{
					for( level=0; level < numlevels_local; level++ )
						sp->st[level].Pop() *=
							dense.xIonDense[nelem][nelem-ipISO] / TotalPop;
				}
				else
				{
					sp->st[0].Pop() =
						MAX2( 1e-30 * dense.xIonDense[nelem][nelem-ipISO],
						dense.xIonDense[nelem][nelem-ipISO] - TotalPopExcited );
				}
				sp->lgPopsRescaled = true;
			}
			ASSERT( sp->st[0].Pop() >= 0. );
		}

	}
	/* all solvers end up here */

	/* check on the sum of the populations */
	if( lgNegPop || dense.xIonDense[nelem][nelem-ipISO] < 0. )
	{
		fprintf( ioQQQ,
			" DISASTER iso_departure_solver finds negative ion fraction or population for iso=%2ld nelem=%2ld "
			"%s using solver %s, IonFrac=%.3e simple=%.3e TE=%.3e ZONE=%4ld\n",
			ipISO,
			nelem,
			elementnames.chElementSym[nelem],
			sp->chTypeAtomUsed,
			dense.xIonDense[nelem][nelem+1-ipISO]/SDIV(dense.xIonDense[nelem][nelem-ipISO]),
			sp->xIonSimple,
			phycon.te,
			nzone );

		fprintf( ioQQQ, " level pop are:\n" );
		for( long i=0; i < numlevels_local; i++ )
		{
			fprintf( ioQQQ,PrintEfmt("%8.1e", sp->st[i].Pop() ));
			if( (i!=0) && !(i%10) ) fprintf( ioQQQ,"\n" );
		}
		fprintf( ioQQQ, "\n" );
		ContNegative();
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	for( ipHi=1; ipHi<numlevels_local; ++ipHi )
	{
		for( ipLo=0; ipLo<ipHi; ++ipLo )
		{
			if( sp->trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
				continue;

			/* population of lower level, corrected for stimulated emission */
			sp->trans(ipHi,ipLo).Emis().PopOpc() =

				sp->st[ipLo].Pop() - sp->st[ipHi].Pop()*
				sp->st[ipLo].g()/sp->st[ipHi].g();

			// don't allow masers from collapsed levels or in Case B only when MASER OFF has been specified
			if( iso_ctrl.lgNoMaser[ipISO][nelem] && ( N_(ipHi) > sp->n_HighestResolved_local || opac.lgCaseB ))
				sp->trans(ipHi,ipLo).Emis().PopOpc() = MAX2( 0., sp->trans(ipHi,ipLo).Emis().PopOpc() );
		}
	}

	// Zero PopOpc of inactive transitions.
	for( ipHi=numlevels_local; ipHi < sp->numLevels_max; ++ipHi )
	{
		for( ipLo=0; ipLo<ipHi; ++ipLo )
		{
			sp->trans(ipHi,ipLo).Emis().PopOpc() = 0.;
		}
	}
	return;

}

int unitary_test_U(long ipISO, long NumMaxLevels, long mmax, long m[][NISO], long U)
{
	static multi_arr<double,2,C_TYPE> zTest,cTest;
	static vector<double> creationTest, yTest;
	static vector<int32> ipivTest;
	int32 nerror=0;
	long Nlev, mmTest, degTest;
	double fval =0.;

	/*
	 *  System: Zx=S
	 *
	 *
	 * Unitary test U1: Z is the identity matrix, S^t =[ 1,1,1]
	 *                  =====> x^t = [1,1,1]
	 * Unitary test U2:
	 *      | 2  1  1  |      | 1 |          |  3  -1  -1  |
	 * Z =  | 1  2  1  |  S = | 1 |  Z^-1 =  | -1   3  -1  |
	 * 		| 1  1  2  |      | 1 |          | -1  -1   3  |
	 *
	 *
	 * 		========>     x^t = [0.25, 0.25, 0.25]
	 *
	 * */


	 if (U == 1)
	 {
		 Nlev = NumMaxLevels;
	 	 mmTest =mmax;
	 	 degTest = 6;
	 	 fval = 1.;
	 }
	 else if (U==2)
	 {
		 Nlev =3;
		 mmTest = 2;
		 degTest = 2;
		 fval = 0.25;
	 }
	 else
	 {
		 fprintf(ioQQQ,"invalid value of U in unitary test for condensation U=%li\n",U);
		 TotalInsanity();
	 }

	vector < vector <double> > LTest(Nlev,vector<double> (Nlev));

	long mTest[mmTest][NISO];

	zTest.alloc(Nlev,Nlev);
	zTest.zero();
	cTest.alloc(mmTest,mmTest);
	cTest.zero();
	yTest.resize(mmTest);
	ipivTest.resize(mmTest);
	creationTest.resize(Nlev);

	if (U ==1)
	{
		for ( long t=0; t<mmTest; t++)
					mTest[t][ipISO] = m[t][ipISO];

		for (long p=0; p<Nlev;p++)
				{
					zTest[p][p]= 1.;
					creationTest[p] =1;
				}
	}
	else if(U==2)
	{
		mTest[0][ipISO] = 0;
		mTest[1][ipISO] = 2;

		for (long p= 0; p<3 ; p++)
		{
			for (long q=0;q<3;q++)
			{
				if (q == p)
					zTest[p][p] = 2.;
				else
					zTest[p][q] = 1.;
			}
			creationTest[p] =1;
		}

	}

	bool TestFail = false;



	set_L(mTest,LTest ,mmTest, ipISO, Nlev,degTest);

	get_condensed_matrix(ipISO, mTest, LTest,zTest,cTest,mmTest, Nlev, degTest );

	for (long t=0; t<mmTest; t++)
		yTest[t]=creationTest[m[t][ipISO]];

	getrf_wrapper(mmTest,mmTest,cTest.data(),mmTest,&ipivTest[0],&nerror);

	getrs_wrapper('N',mmTest,1,cTest.data(),mmTest,&ipivTest[0],&yTest[0],mmTest,&nerror);

	get_interpolated_vector(ipISO, creationTest,yTest,mTest,LTest,mmTest, Nlev,degTest );


	for (long j = 0; j< Nlev; j++)
	{

		if (abs(creationTest[j] - fval) >1e-5 )
		{
			fprintf(ioQQQ,"%li abs(creationTest[j] - fval) %g\n",j,abs(creationTest[j] - fval));
			TestFail = true;
		}
	}

	if (TestFail)
	{
		fprintf(ioQQQ,"U%li test failed\n",U);
		for (long j=0; j < Nlev; j++)
			fprintf(ioQQQ,"x[%li] %g fval %g\n",j,creationTest[j], fval);

		if (U==2)
		{
			fprintf(ioQQQ,"################test U2 - matrix Z #######\n");
			for (long p= 0; p<3 ; p++)
			{
				for (long q=0;q<3;q++)
				{
					if (q <2 )
						fprintf(ioQQQ," [%li][%li] = %g\t",p,q,zTest[p][q]);
					else
						fprintf(ioQQQ," [%li][%li] = %g\n",p,q,zTest[p][q]);
				}
			}
			fprintf(ioQQQ,"################test U2 - matrix C #######\n");

			for (long p= 0; p<2 ; p++)
			{
				for (long q=0;q<2;q++)
				{
					if (q <1 )
						fprintf(ioQQQ," [%li][%li] = %g\t",p,q,cTest[p][q]);
					else
						fprintf(ioQQQ," [%li][%li] = %g\n",p,q,cTest[p][q]);
				}
			}

		}

		return 0;
	}
	else
		return 1;

}

void set_L(long m[][NISO], vector< vector <double> >& L, long mmax, long ipISO, long NumMaxLevels, long deg)
{
	double num=1., denom=1.;

	long ir = (deg/2) -1;

	if(mmax ==1 )
	{
		L[0][0] =1;
		fprintf(ioQQQ,"Warning! Set_L() value for mmax is 1, is this a test? %li",mmax);
		return;
	}
	else if(mmax-1 < 1 || mmax < deg)
	{
		fprintf(ioQQQ,"Set_L() incorrect value for mmax %li or the degree of L %li",mmax, deg);
		cdEXIT(EXIT_FAILURE);
	}

	for ( long i=0; i < NumMaxLevels; i++)
		for (long t =0; t < mmax; t++)
			L[i][t] = 0.;

	for ( long i=0; i < NumMaxLevels; i++)
	{
		for (long t =0; t < mmax; t++)
		{
			long ttmax = (t > mmax-(deg-1)) ? mmax : t+deg-ir;
			long ttmin = (t > 1 ) ? t -ir : 0 ;
			if (t > mmax - deg +ir)
			{
				ttmax = mmax;
				ttmin = mmax-deg;
			}
			if(ttmin<0)
				ttmin = 0;

			if (m[t+1][ipISO] - m[t][ipISO] > 1
					&& i>= m[t][ipISO]+1 && i <= m[t+1][ipISO] && t+1 < mmax)
			{
				for (long intp = ttmin ; intp < ttmax ; intp ++ )
				{
					denom =1. ,num = 1.;

						for (long tau = ttmin; tau < ttmax; tau++)
						{

							if (tau != intp)
							{
								denom *= (double)(m[tau][ipISO] - m[intp][ipISO]);
								num *=(double)(m[tau][ipISO] - i);
							}

						}
						L[i][intp] = num/denom;
				}
			}
			else if(i == m[t][ipISO])
			{
				L[i][t] = 1.;
			}
		}

		/*security check */
		bool nonull = false;
		for (long t = 0 ; t<mmax ; t++)
		{
			if (L[i][t] != 0)
			{
				nonull = true;
				break;
			}
		}
		if (!nonull)
		{
			fprintf(ioQQQ,"Set_L() all zeros in Lagrangian for element %li \n",i);
			cdEXIT(EXIT_FAILURE);
		}


	}
	/* Debugging print option
	 *

	for ( long t=0; t < NumMaxLevels; t++)
	{
				for (long i =0; i < NumMaxLevels; i++)
		{
			if(L[i][t]!=0)
				fprintf(ioQQQ," Set_L L [%li][%li] = %g m[t] %li\n ",i,t,L[i][t],m[t][ipISO]);
		}
	}
	*
	*/

	return;
}

/* This condensed the matrix using an external interpolator L*/
void get_condensed_matrix( long ipISO, long m[][NISO],vector < vector <double> >& L, multi_arr<double,2,C_TYPE> z,
		multi_arr<double,2,C_TYPE>& c, long mmax , long NumMaxLevels, long deg)
{
	/* use Brocklehurst 4.4 */
	long ir = (deg/2)-1;

	c.zero();
	/* first fill the matrix with the condensed indices */
	for ( long s=0; s< mmax; s++)
	{
		for( long t=0; t< mmax ; t++)
			c[s][t] = z[m[s][ipISO]][m[t][ipISO]];

	}

	for ( long s=0; s< mmax; s++)
	{


		if (m[s+1][ipISO] - m[s][ipISO] > 1)
		{
			long ssmax = (s > mmax-(deg-1)) ? mmax : s+deg-ir;
			long ssmin = (s > 1 ) ? s -2 : 0 ;
			if (s > mmax -deg + ir)
			{
				ssmax = mmax;
				ssmin = mmax-deg;
			}
			for ( long j=0; j < NumMaxLevels; j++)
			{

				if (j>= m[s][ipISO]+1 && j < m[s+1][ipISO] )
				{
					for (long t = 0; t< mmax ; t++)
					{
						for (long intp = ssmin ; intp < ssmax ; intp ++ )
						{
							long k = m[intp][ipISO];
							c[intp][t] += z[k][j]*L[j][t];
						}

					}
				}
			}
		}
	}
	return;
}

/*This intepolates the input vector using an external interpolator L */
void get_interpolated_vector(long ipISO, vector <double>& xexpan,vector <double>& xcon, long m[][NISO],vector < vector <double> >& L,
		long mmax, long NumMaxLevels, long deg)
{
	long ir = (deg/2)-1;

	for ( long j = 0; j < NumMaxLevels; j++)
		xexpan[j]=0.;
	long p =0;
	for ( long t=0; t< mmax ; t++)
	{
		xexpan[m[t][ipISO]] = xcon[t];

		if (m[t+1][ipISO] - m[t][ipISO] >1)
		{
			p=t;
			break;
		}
	}
	for ( long t=p; t< mmax ; t++)
	{
		if (m[t+1][ipISO] - m[t][ipISO] >1 || t >= mmax-2)
		{

			long ttmax = (t > mmax-(deg-1)) ? mmax : t+deg-ir;
			long ttmin = (t > 1 ) ? t -2 : 0 ;
			if (t > mmax -deg+ir)
			{
				ttmax = mmax;
				ttmin = mmax-deg;
			}

			for (long intp = ttmin ; intp < ttmax ; intp ++ )
			{

				for ( long j = m[t][ipISO]+1; j <= m[t+1][ipISO]; j++)
				{
					if (j< NumMaxLevels)
						xexpan[j] += L[j][intp]*xcon[intp];
					else
						break;
				}

			}

		}
		else
			xexpan[m[t+1][ipISO]] = xcon[t+1];
	}

	return;

}
