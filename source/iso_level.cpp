/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_level solve for iso-sequence level populations */
#include "cddefines.h"
#include "atmdat.h"
#include "continuum.h"
#include "conv.h"
#include "doppvel.h"
#include "dynamics.h"
#include "elementnames.h"
#include "grainvar.h"
#include "he.h"
#include "hydroeinsta.h"
#include "ionbal.h"
#include "iso.h"
#include "iso_condensed.h"
#include "opacity.h"
#include "phycon.h"
#include "rfield.h"
#include "trace.h"
#include "mole.h"
#include "freebound.h"
#include "two_photon.h"
#include "dense.h"
#include "vectorize.h"
#include "prt.h"
#include "iterations.h"


STATIC void iso_multiplet_opacities_one(
	const long int ipISO, const long int nelem);

void set_L(long m[][NISO], vector <vector <double>>& L, long mmax, long ipISO, long NumMaxLevels, long deg);

void get_condensed_matrix( long ipISO, long m[][NISO],vector < vector <double> >& L, multi_arr<double,2,C_TYPE> z, multi_arr<double,2,C_TYPE> &c,
		long mmax, long NumMaxLevels, long deg);
void get_interpolated_vector(long ipISO, vector <double>& xexpan, vector <double>& xcon,long m[][NISO],vector < vector <double> >& L,
		long mmax, long NumMaxLevels, long deg);
int unitary_test_U(long ipISO, long NumMaxLevels, long mmax, long m[][NISO], long U);

multi_arr<double,2,C_TYPE> z, SaveZ,c;

// vector <vector <double> >&
/* solve for level populations  */
void iso_level( const long int ipISO, const long int nelem, double &renorm,
		bool lgPrtMatrix )
{
	long int ipHi,
		ipLo,
		i,
		level,
		level_error;

 	t_iso_sp* sp = &iso_sp[ipISO][nelem]; 
	const long int numlevels_local = sp->numLevels_local;

	double BigError;

	int32 nerror;
	bool lgNegPop=false;
	long mmax=1;// = 8;
	static vector<int32> ipiv;
	/* this block of variables will be obtained and freed here */
	double source=0., sink=0.;
	static vector<double> PopPerN;
	PopPerN.resize(sp->n_HighestResolved_local+1);

	DEBUG_ENTRY( "iso_level()" );
	renorm = 1.;

	/* check that we were called with valid charge */
	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );

	/* these two collision rates must be the same or we are in big trouble,
	 * since used interchangeably */
	ASSERT( ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0]< SMALLFLOAT ||
		fabs( (sp->fb[0].ColIoniz* dense.EdenHCorr) /
		SDIV(ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0] ) - 1.) < 0.001 );


	/* obtain here for the moment the number and the index of pivotal levels */

	long m[sp->numLevels_max][NISO];

	if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
	{
		for (long ipHi =0 ; ipHi < sp->mmax; ipHi++)
			m[ipHi][ipISO] =sp->fb[ipHi].ipm;
		mmax = sp->mmax;

		ipiv.resize(sp->mmax);
	}
	else
		ipiv.resize(numlevels_local);

	/*for ( long j = 0; j <mmax; j++)
		fprintf(ioQQQ," nelem %li m ( j %li , ipISO %li ) = %li mmax %li\n", nelem, j, ipISO, m[j][ipISO], mmax); // debug!!!*/




	vector < vector <double> > L(sp->numLevels_max,vector<double> (sp->numLevels_max,0.));

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
			fprintf( ioQQQ, "     iso_level iso=%2ld nelem=%2ld simple II/I=%10.2e so not doing equilibrium, doing %s.\n", 
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
		/* fill in recombination vector - values were set in iso_ionize_recombine.cpp */
		static vector<double> creation, error, work, Save_creation, y;//, qterm, zterm;
		creation.resize(numlevels_local);
		error.resize(numlevels_local);
		work.resize(numlevels_local);
		//Save_creation.resize(numlevels_local);
		y.resize(mmax);

		ASSERT( dense.xIonDense[nelem][nelem+1-ipISO] >= 0.f );
		for( level=0; level < numlevels_local; level++ )
		{
			/* total recombination from once more ionized [cm-3 s-1] */
			creation[level] = sp->fb[level].RateCont2Level * dense.xIonDense[nelem][nelem+1-ipISO];
		}
		
		double ctsource=0.0, ctsink=0.0, ctrec=0.0;
		/* now charge transfer - all into/from ground, two cases, H and not H */
		if( nelem==ipHYDROGEN )
		{
			/* charge transfer, hydrogen onto everything else */
			/* charge exchange recombination */
			ctrec += atmdat.CharExcRecTotal[nelem];
			ctsink += atmdat.CharExcIonTotal[nelem];
		}
		else if( nelem==ipHELIUM && ipISO==ipHE_LIKE )
		{
			/* this is recom of He due to ct with all other gas constituents */
			ctrec += atmdat.CharExcRecTotal[nelem];
			ctsink += atmdat.CharExcIonTotal[nelem];
		}
		else
		{
			for (long nelem1=ipHYDROGEN; nelem1 < t_atmdat::NCX; ++nelem1)
			{
				long ipISO=nelem1;
				ctrec += atmdat.CharExcRecTo[nelem1][nelem][nelem-ipISO]*iso_sp[ipISO][nelem1].st[0].Pop();
				ctsink += atmdat.CharExcIonOf[nelem1][nelem][nelem-ipISO]*dense.xIonDense[nelem1][1];
			}					
		}
		ctsource += ctrec*dense.xIonDense[nelem][nelem+1-ipISO];
		
		if ( nelem > ipISO )
		{
			double ction=0.0;
			if( nelem==ipHELIUM )
			{
				ctsink += atmdat.CharExcRecTotal[nelem];
				ction += atmdat.CharExcIonTotal[nelem];
			}
			else
			{
				for (long nelem1=ipHYDROGEN; nelem1 < t_atmdat::NCX; ++nelem1)
				{
					long ipISO1=nelem1;
					ctsink += atmdat.CharExcRecTo[nelem1][nelem][nelem-(ipISO+1)]*iso_sp[ipISO1][nelem1].st[0].Pop();
					ction += atmdat.CharExcIonOf[nelem1][nelem][nelem-(ipISO+1)]*dense.xIonDense[nelem1][1];
				}
			}
			ctsource += ction * dense.xIonDense[nelem][nelem-(ipISO+1)];
		}
		
		/* now do the 2D array */
//		multi_arr<double,2,C_TYPE> z, SaveZ,c;
		z.alloc(numlevels_local,numlevels_local);
		z.zero();

		c.alloc(mmax,mmax);
		c.zero();
		vector <double> xin,xout;
		xin.resize(mmax);
		xout.resize(numlevels_local);

		for (i=0;i<numlevels_local;i++)
			xout[i] = 0.;
		for (i=0; i < mmax ; i++)
			xin[i] = 0.;

		/* this branch is main solver, full level populations 
		 * assert since this code must change if NISO ever increased */
		ASSERT( NISO == 2 );

		strcpy( sp->chTypeAtomUsed, "Popul" );
		if( trace.lgTrace && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "     iso_level iso=%2ld nelem=%2ld doing regular matrix inversion, %s\n", 
				ipISO, nelem, sp->chTypeAtomUsed );
		}

		//qList::const_iterator StElm = StatesElemNEW[nelem][nelem-ipISO].begin();
		qList::const_iterator StElm = sp->st.begin();
		static vector<double> Boltzmann_overg;
		Boltzmann_overg.resize(numlevels_local-1);
		for (unsigned lev = 0; lev < Boltzmann_overg.size(); ++lev)
			Boltzmann_overg[lev] = 1.0/(double)StElm[lev].g();
		
		/* >>chng 05 dec 21, rm eden to make into rate coefficient */
		sp->qTot2S = sp->fb[1].ColIoniz;

		static vector<double> coll_down, RadDecay, pump;
		coll_down.resize(numlevels_local);
		RadDecay.resize(numlevels_local);
		pump.resize(numlevels_local);
		avx_ptr<double> arg(1,numlevels_local), bstep(1,numlevels_local);

		for( level=1; level < numlevels_local; level++ )
		{
			fixit("Wouldn't need to mask this out if levels were in order");
			arg[level] = -max((StElm[level].energy().K()-StElm[level-1].energy().K()) / phycon.te, -38.);
		}
		vexp( arg.ptr0(), bstep.ptr0(), 1, numlevels_local );

		ColliderDensities colld(colliders);

		enum { DEBUG_RATES = false };

		if( DEBUG_RATES && iterations.lgLastIt )
			fprintf( stdout, "# ipISO\tnelem\tlevel\tcollDown\tcollIonz\tradRecom\tRadDecay\n" );

		/* master balance equation, use when significant population */
		for( level=0; level < numlevels_local; level++ )
		{
			double coll_down_total = 0.;

			/* all process depopulating level and placing into the continuum
			 * this does NOT include grain charge transfer ionization, added below */
			z[level][level] += sp->fb[level].RateLevel2Cont;

			if (level != 0)
			{
				for ( ipLo = 0; ipLo < level; ++ipLo )
					Boltzmann_overg[ipLo] *= bstep[level];
			}

			/* all processes populating level from below */
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
					pump[ipLo] = MAX2( iso_ctrl.SmallA, sp->trans(level,ipLo).Emis().pump() );
				}
			}

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
					Boltzmann_overg[ipLo];

				z[ipLo][ipLo] += coll_up + pump[ipLo] ;
				z[ipLo][level] -= coll_up + pump[ipLo] ;

				double pump_down = pump[ipLo] *
					(double)StElm[ipLo].g() * rglev;

/*				if(!iso_ctrl.lgMatCondOff[ipISO][nelem])
					z[ipLo][level] -= (RadDecay[ipLo] + pump_down)*Boltzmann_overg[ipLo]*glev;
*/

				z[level][level] += RadDecay[ipLo] + pump_down + coll_down[ipLo];
				z[level][ipLo] -= RadDecay[ipLo] + pump_down + coll_down[ipLo];
/*				if(!iso_ctrl.lgMatCondOff[ipISO][nelem])
				{
					z[level][ipLo] += RadDecay[ipLo];
				}
*/

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

		/** \todo 2 the indices for the two-photon rates must be changed for further iso sequences. */  
		ASSERT( ipISO <= ipHE_LIKE );
		for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
		{
			// induced two photon emission - special because upward and downward are
			// not related by ratio of statistical weights 
			// iso.lgInd2nu_On is controlled with SET IND2 ON/OFF command 

			fixit("need Pesc or the equivalent to multiply AulTotal?");
			// downward rate
/*			if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
				z[tnu->ipLo][tnu->ipHi] -= tnu->AulTotal*Boltzmann_overg[tnu->ipLo]*(double)StElm[tnu->ipHi].g();
			else
*/
			z[tnu->ipHi][tnu->ipLo] -= tnu->AulTotal;
			z[tnu->ipHi][tnu->ipHi] += tnu->AulTotal; 

		}
		
		if (iso_ctrl.lgInd2nu_On)
		{
			for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
			{
				// induced two photon emission - special because upward and downward are
				// not related by ratio of statistical weights 
				// iso.lgInd2nu_On is controlled with SET IND2 ON/OFF command 


/*				if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
					z[tnu->ipLo][tnu->ipHi] -= tnu->induc_dn*Boltzmann_overg[tnu->ipLo]*(double)StElm[tnu->ipHi].g();
				else
*/
				z[tnu->ipHi][tnu->ipLo] -= tnu->induc_dn;
				z[tnu->ipHi][tnu->ipHi] += tnu->induc_dn;
				
				// upward rate
				z[tnu->ipLo][tnu->ipHi] -= tnu->induc_up;
				z[tnu->ipLo][tnu->ipLo] += tnu->induc_up;
			}
		}
		/* HION_LTE_POP	is planck^2 / (2 pi m_e k ), raised to 3/2 next */
		//double factor = HION_LTE_POP*dense.AtomicWeight[nelem]/
		//		(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);


				/* term in () is stat weight of electron * ion */
		//double ConvLTEPOP = powpq(factor,3,2)/(2.*iso_ctrl.stat_ion[ipISO])/phycon.te32;
		/* grain charge transfer recombination and ionization to ALL other stages */
		for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
		{
			if( ion!=nelem-ipISO )
			{


				/*AQUI,AQUI*/ /*METER ESTOS TERMINOS DESPUES YA QUE DEPENDEN DEL NIVEL*/
/*				if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
				{
					for( level=0; level < numlevels_local; level++ )
						creation[level] += gv.GrainChTrRate[nelem][ion][nelem-ipISO] / ConvLTEPOP/sp->st[level].ConBoltz()
							+ gv.GrainChTrRate[nelem][nelem-ipISO][ion];
				}
*/

				source += gv.GrainChTrRate[nelem][ion][nelem-ipISO]*dense.xIonDense[nelem][ion];

				sink += gv.GrainChTrRate[nelem][nelem-ipISO][ion];
/*
				double xmolechtrrate = mole.xMoleChTrRate[nelem][ion][nelem-ipISO];
*/
				/*AQUI AQUI*/
				/*double xmolesink = mole.xMoleChTrRate[nelem][nelem-ipISO][ion];
*/
/*				if (!iso_ctrl.lgMatCondOff[ipISO][nelem] )
				{
					for( level=0; level < numlevels_local; level++ )
					creation[level] += xmolechtrrate / ConvLTEPOP/ sp->st[level].ConBoltz()*atmdat.lgCTOn
					+  xmolesink* atmdat.lgCTOn;
				}
*/

				source	+= mole.xMoleChTrRate[nelem][ion][nelem-ipISO]* dense.xIonDense[nelem][ion]* atmdat.lgCTOn;

				sink += mole.xMoleChTrRate[nelem][nelem-ipISO][ion] * atmdat.lgCTOn;

			}
		}
		
		/* add in source and sink terms from molecular network. */
		/*if (!iso_ctrl.lgMatCondOff[ipISO][nelem] && 0)
		{
			/for( level=0; level < numlevels_local; level++ )
				creation[level] +=mole.source[nelem][nelem-ipISO]/ ConvLTEPOP/ sp->st[level].ConBoltz()/ dense.xIonDense[nelem][nelem-ipISO]
			+mole.sink[nelem][nelem-ipISO];
		}
*/

		source += mole.source[nelem][nelem-ipISO];
		sink += mole.sink[nelem][nelem-ipISO];

#if	1
		/* >>chng 02 Sep 06 rjrw -- all elements have these terms */
		/*>>>chng 02 oct 01, only include if lgAdvection or lgTimeDependentStatic is set */
		if( iteration > dynamics.n_initial_relax+1 &&
			( dynamics.lgAdvection || dynamics.lgTimeDependentStatic ) &&
			dynamics.Rate != 0.0 &&
			!dynamics.lgEquilibrium && dynamics.lgISO[ipISO] )
		{
			/* add in advection - these terms normally zero */
			/*if (!iso_ctrl.lgMatCondOff[ipISO][nelem] )
			{
				for( level=0; level < numlevels_local; level++ )
				creation[level] += dynamics.Source[nelem][nelem-ipISO]/ ConvLTEPOP/ sp->st[level].ConBoltz()/ dense.xIonDense[nelem][nelem-ipISO] +
				dynamics.Rate;
			}
*/

			source += dynamics.Source[nelem][nelem-ipISO];
			/* >>chng 02 Sep 06 rjrw -- advective term not recombination */
			sink += dynamics.Rate;
		}
#else
		/*>>>chng 02 oct 01, only include if lgAdvection or lgTimeDependentStatic is set */
		if( iteration > dynamics.n_initial_relax+1 &&
			( dynamics.lgAdvection || dynamics.lgTimeDependentStatic ) &&
			dynamics.Rate != 0.0 &&
			!dynamics.lgEquilibrium && dynamics.lgISO[ipISO])
		{
			for( level=0; level < numlevels_local; level++ )
			{
				creation[level] += dynamics.StatesElem[nelem][nelem-ipISO][level];
				z[level][level] += dynamics.Rate;
			}
		}
#endif

		/* ionization from/recombination from lower ionization stages */
		for(long ion_from=dense.IonLow[nelem]; ion_from < MIN2( dense.IonHigh[nelem], nelem-ipISO ) ; ion_from++ )
		{
			if( ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] >= 0. )
			{
				/*if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
				{
					for( level=0; level < numlevels_local; level++ )
					{
						creation[level]+=ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] / ConvLTEPOP/ sp->st[level].ConBoltz();
						if( ion_from == nelem-1-ipISO )
							creation[level] += ionbal.RateRecomTot[nelem][ion_from];
					}
				}
*/

				/* ionization from lower ionization stages, cm-3 s-1 */
				source += ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] * dense.xIonDense[nelem][ion_from];
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
			/*for( level=0; level < numlevels_local; level++ )
			{
				z[level][level] += sink;
			}*/
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
				/*if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
					creation[level] += source/ ConvLTEPOP/sp->st[level].Boltzmann()/dense.xIonDense[0][1] + sink;

							 //  *= sp->st[level].Boltzmann()*sp->st[level].g()/partfun;
				else
*/
				creation[level] += source*
						sp->st[level].Boltzmann()*sp->st[level].g();
				z[level][level] += sink;
				}
		}
		/*if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
			creation[0] += ctsource/ ConvLTEPOP/ sp->st[0].Boltzmann()/dense.xIonDense[0][1]+ctsink;
		else
*/
		creation[0] += ctsource;
		z[0][0] += ctsink;

		/* >>chng 04 nov 30, atom XX-like collisions off turns this off too */
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
				z[level][ipH1s] -= RateDown;

				/* rate from ground to upper level */
				z[ipH1s][level] -= RateUp;

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
				fprintf( ioQQQ, "  pop level     others => (iso_level)\n" );
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
		/*
			fprintf(ioQQQ," recomb ct %.2e co %.2e hectr %.2e hctr %.2e\n",
					  atmdat.CharExcRecTotal[ipHELIUM],
					  findspecieslocal("CO")->den ,
					  atmdat.CharExcRecTo[ipHELIUM][nelem][nelem-ipISO]*iso_sp[ipHE_LIKE][ipHELIUM].st[0].Pop() ,
					  atmdat.CharExcRecTo[ipHYDROGEN][nelem][nelem-ipISO]*iso_sp[ipH_LIKE][ipHYDROGEN].st[0].Pop() );
			fprintf(ioQQQ," recomb          ");
			*/
		if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO])) && iso_ctrl.lgMatCondOff[ipISO][nelem] )
		{
			for( long n=0; n < NumMaxLevels; n++ )
			{
				fprintf( ioQQQ,"\t%.9e", creation[n] );
			}

			fprintf( ioQQQ, "\n" );
		}




		if( lgPrtMatrix || 1 )
		{
			valarray<double> cr( get_ptr(creation), creation.size() );
			prt.matrix.prtRates( numlevels_local, z, cr );
		}

		nerror = 0;

		if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
		{

			/*qterm.resize(mmax);
			zterm.resize(mmax);*/

			for( long j=0; j<mmax; j++)
			{

				y[j]=creation[m[j][ipISO]];
				//fprintf(ioQQQ,"ipISO %li, %li j%li\n",nelem,m[j][ipISO],j);
				//fprintf(ioQQQ," qterm %g\n",sp->fb[m[j][ipISO]].SourceBound);
				//for (long i=0;i<4;i++)
					//fprintf(ioQQQ," zterm %g i %li\n",sp->fb[sp->fb[j].ipm].SinkBound[i],i);
				//fprintf(ioQQQ," j %li vector %g\n",j,y[j]);
				if (1)
				{
					vector<long> mj;

					get_mj(mmax, mj);

					y[j] += sp->fb[j].SourceBound;
					fprintf(ioQQQ,"sources j %li y %e sbound %e %e %e %e %e\n",j,y[j],sp->fb[j].SourceBound,
							sp->fb[j].SinkBound[0],sp->fb[j].SinkBound[1],sp->fb[j].SinkBound[2],sp->fb[j].SinkBound[3]);

					for (long i =0; i<4 ; i++)
					{
						//if (j != mj[i] )
							c[j][mj[i]] -= sp->fb[j].SinkBound[i];
						//else
							//c[j][j] += sp->fb[j].SinkBound[i];
					}
				}
			}

			for( ipLo=0; ipLo < mmax; ipLo++ )
				Save_creation[ipLo] = y[ipLo];

			if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO])))
			{
				fprintf(ioQQQ,"mmax %li\n",mmax);

				for( long n=0; n < mmax; n++ )
				{
					fprintf( ioQQQ,"\t%.9e", y[n] );
				}
			}

			getrf_wrapper(mmax,mmax,c.data(),mmax,&ipiv[0],&nerror);

			getrs_wrapper('N',mmax,1,c.data(),mmax,&ipiv[0],&y[0],mmax,&nerror);


		}
		else
		{
			getrf_wrapper(numlevels_local,numlevels_local,
					z.data(),numlevels_local,&ipiv[0],&nerror);

			getrs_wrapper('N',numlevels_local,1,z.data(),numlevels_local,&ipiv[0],
					&creation[0],numlevels_local,&nerror);
		}
		
		if( nerror != 0 )
		{
			fprintf( ioQQQ, " iso_level: dgetrs finds singular or ill-conditioned matrix\n" );
			cdEXIT(EXIT_FAILURE);
		}

		//output, solution condensed vector
		if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
		{
			fprintf( ioQQQ,"mmax %li\n",mmax);
			for( long n=0; n < mmax; n++ )
			{
				fprintf( ioQQQ,"x[%li] m %li %.9e\n",n, m[n][ipISO], y[n] ); // debug!!!
			}
		}

		//return to extended matrix
		if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
		{
			get_interpolated_vector(ipISO, creation,y,m,L,mmax, NumMaxLevels, deg );
		}

		long InvertMatrixLevels = NumMaxLevels;
		if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
			InvertMatrixLevels = mmax;
/*		for (long j = 0; j< sp->numLevels_local; j++)
			fprintf(ioQQQ,"j %li creation[j] %g\n",j,creation[j]); //debug!!! */
		/* check whether solution is valid */
		/* >>chng 06 aug 28, both of these from numLevels_max to _local. */
		for( level=ipH1s; level < InvertMatrixLevels; level++ )
		{
			double qn = 0., qx = 0.;
			error[level] = 0.;
			for( long n=ipH1s; n < InvertMatrixLevels; n++ )
			{
				double q= SaveZ[n][level];

				if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
					q *=y[n];
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
		BigError = -1.;
		level_error = -1;
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

		/* matrix inversion should be nearly as good as the accuracy of a double,
		 * but demand that it is better than epsilon for a float */
		if( BigError > FLT_EPSILON ) 
		{
			if( !conv.lgSearch )
				fprintf(ioQQQ," PROBLEM" );

			fprintf(ioQQQ,
				" iso_level zone %.2f - largest residual in iso=%li %s nelem=%li matrix inversion is %g "
				"level was %li Search?%c \n", 
				fnzone,
				ipISO,
				elementnames.chElementName[nelem],
				nelem , 
				BigError , 
				level_error,
				TorF(conv.lgSearch) );
		}
		
		// Force level balance to LTE
		if ( iso_ctrl.lgLTE_levels[ipISO] )
		{
			t_iso_sp* sp = &iso_sp[ipISO][nelem];
			double partfun=0.0;
			for ( level = 0; level < numlevels_local; level++ )
			{
				partfun += sp->st[level].Boltzmann()*sp->st[level].g();
			}
			double scale = dense.xIonDense[nelem][nelem-ipISO]/partfun;
			for ( level = 0; level < numlevels_local; level++ )
			{
				creation[level] = sp->st[level].Boltzmann()*sp->st[level].g()*scale;
			}
		}

		double creationzero = creation[0];
		for( level = numlevels_local-1; level >= 0; --level )
		{
			/* check for negative populations */
			fprintf(ioQQQ,"creation[%li] = %g \n",level,creation[level]);

			if(creationzero < 0.)
			{
				creation[level] *= -1.;
				if (creation[level] < 0.)
					fprintf(ioQQQ,"no cambio signo creation %g level %li\n",creation[level],level);
			}
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
			creation[0] = dense.xIonDense[nelem][nelem-ipISO];
			// flip flag back
			lgNegPop = false;
		}

		/* put level populations into master array */
		for( level=0; level < numlevels_local; level++ )
		{
			if (creation[level] < 0)
				fprintf(ioQQQ,"creation %g level %li\n",creation[level],level);
			ASSERT( creation[level] >= 0. );
/*			if (!iso_ctrl.lgMatCondOff[ipISO][nelem])
 */
			sp->st[level].Pop() = creation[level];
/*
			else
				sp->st[level].Pop() = creation[level]*dense.eden*dense.xIonDense[0][1]*ConvLTEPOP*sp->st[level].Boltzmann()*sp->st[level].g()/2.;
*/
			if( sp->st[level].Pop() <= 0 && !conv.lgSearch )
			{
				fprintf(ioQQQ,
					"PROBLEM non-positive level pop for iso = %li, nelem = "
					"%li = %s, level=%li val=%.3e nTotalIoniz %li\n", 
					ipISO,
					nelem , 
					elementnames.chElementSym[nelem],
					level,
					sp->st[level].Pop() ,
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
		/* renormalize the populations to agree with ion solver */
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
			" DISASTER iso_level finds negative ion fraction for iso=%2ld nelem=%2ld "
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
		for( i=0; i < numlevels_local; i++ )
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

/** update multiplet opacities */
void iso_multiplet_opacities( void )
{
	for (long nelem = ipHYDROGEN; nelem < LIMELM; ++nelem)
	{
		if( ! dense.lgElmtOn[nelem] )
			continue;
		for ( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			if( (dense.IonHigh[nelem] >= nelem - ipISO) &&
				 (dense.IonLow[nelem] <= nelem - ipISO) )
			{
				iso_multiplet_opacities_one(ipISO, nelem);
			}
		}
	}
}

STATIC void iso_multiplet_opacities_one(
	const long int ipISO, const long int nelem)
{
	DEBUG_ENTRY( "iso_multiplet_opacities_one()" );
 	t_iso_sp* sp = &iso_sp[ipISO][nelem]; 
	// Allocate array for calculating n-n' multiplet opacities.
	multi_arr<double,2> MultOpacs;
	long nMax = sp->n_HighestResolved_max + sp->nCollapsed_max;
	MultOpacs.reserve( nMax+1 );
	for( long n=2; n <= nMax; ++n )
		MultOpacs.reserve( n, n+1 );
	MultOpacs.alloc();
	MultOpacs.zero();
	
	double rDopplerWidth = 1.0/GetDopplerWidth(dense.AtomicWeight[nelem]);
	
	// Sum n-n' multiplet opacities.	
	for( long ipHi=1; ipHi < sp->numLevels_max; ++ipHi )
	{
		for( long ipLo=0; ipLo < ipHi; ++ipLo )
		{
			const TransitionProxy& tr = sp->trans(ipHi,ipLo);
			MultOpacs[ sp->st[ipHi].n() ][ sp->st[ipLo].n() ] += tr.Emis().PopOpc() *
				tr.Emis().opacity() * rDopplerWidth;
		}
	}

	// Now store n-n' multiplet opacities.
	for( long ipHi=1; ipHi < sp->numLevels_max; ++ipHi )
	{
		for( long ipLo=0; ipLo < ipHi; ++ipLo )
		{
			const TransitionProxy& tr = sp->trans(ipHi,ipLo);
			tr.Emis().mult_opac() = MultOpacs[ sp->st[ipHi].n() ][ sp->st[ipLo].n() ];
		}
	}
}

void iso_set_ion_rates( long ipISO, long nelem)
{
	DEBUG_ENTRY( "iso_set_ion_rates()" );
 	t_iso_sp* sp = &iso_sp[ipISO][nelem]; 
	const long int numlevels_local = sp->numLevels_local;
	/* this is total ionization rate, s-1, of this species referenced to
	 * the total abundance */
	double TotalPop = 0.;
	ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] = 0.;
	for( long level=0; level < numlevels_local; level++ )
	{
		/* sum of all ionization processes from this atom to ion, cm-3 s-1 now,
		 * but is divided by TotalPop below to become s-1 */
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] += 
			sp->st[level].Pop() * sp->fb[level].RateLevel2Cont;
		TotalPop += sp->st[level].Pop();
	}
	
	if( ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] > BIGDOUBLE )
	{
		fprintf( ioQQQ, "DISASTER RateIonizTot for Z=%li, ion %li is larger than BIGDOUBLE.  This is a big problem.",
					nelem+1, nelem-ipISO);
		cdEXIT(EXIT_FAILURE);
	}

	if (TotalPop <= SMALLFLOAT)
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] = sp->fb[0].RateLevel2Cont;
	else
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] /= TotalPop;

	if( ionbal.RateRecomIso[nelem][ipISO] > 0. )
		sp->xIonSimple = ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1]/ionbal.RateRecomIso[nelem][ipISO];
	else
		sp->xIonSimple = 0.;

	ASSERT( ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] >= 0. );

	if( ipISO == ipHE_LIKE && nelem==ipHELIUM && nzone>0 )
	{
		/* find fraction of He0 destructions due to photoionization of 2^3S */
		double ratio;
		double rateOutOf2TripS = sp->st[ipHe2s3S].Pop() * sp->fb[ipHe2s3S].RateLevel2Cont;
		if( rateOutOf2TripS > SMALLFLOAT )
		{
			ratio = rateOutOf2TripS /
				( rateOutOf2TripS + sp->st[ipHe1s1S].Pop()*ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] );
		}
		else
		{
			ratio = 0.;
		}
		if( ratio > he.frac_he0dest_23S )
		{
			/* remember zone where this happended and fraction, and frac due to photoionization */
			he.nzone = nzone;
			he.frac_he0dest_23S = ratio;
			rateOutOf2TripS = sp->st[ipHe2s3S].Pop() *sp->fb[ipHe2s3S].gamnc;
			if( rateOutOf2TripS > SMALLFLOAT )
			{
				he.frac_he0dest_23S_photo = rateOutOf2TripS /
					( rateOutOf2TripS + sp->st[ipHe1s1S].Pop()*ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] );
			}
			else
			{
				he.frac_he0dest_23S_photo = 0.;
			}
		}
	}
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
	/*
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
		yTest[t]=creationTest[mTest[t][ipISO]];

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
