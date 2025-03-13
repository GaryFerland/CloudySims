/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_create create data for hydrogen and helium, 1 per coreload, called by ContCreatePointers 
 * in turn called after commands parsed */
/*iso_zero zero data for hydrogen and helium */
#include "cddefines.h"
#include "atmdat_adfa.h"
#include "dense.h"
#include "helike.h"
#include "helike_einsta.h"
#include "hydro_bauman.h"
#include "hydrogenic.h"
#include "hydroeinsta.h"
#include "hydro_tbl.h"
#include "iso.h"
#include "opacity.h"
#include "phycon.h"
#include "taulines.h"
#include "mole.h"
#include "freebound.h"
#include "lines_service.h"
#include "prt.h"
#include "save.h"
#include "rfield.h"

/*iso_zero zero data for hydrogen and helium */
STATIC void iso_zero(void);

/* allocate memory for iso sequence structures */
STATIC void iso_allocate(void);

/* define levels of iso sequences and assign quantum numbers to those levels */
STATIC void iso_assign_quantum_numbers(void);

STATIC void iso_assign_extralyman_levels();

STATIC void FillExtraLymanLine( const TransitionList::iterator& t, long ipISO, long nelem, long nHi, double j );

STATIC void iso_satellite( void );

static const char chL[21] = {'S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V','W','X','Y','Z'};

static vector<species> isoSpecies( NISO );

int getL(char l)
{
	char ll = toupper(l);
	for( int i=0; i < 21; ++i )
		if( ll == chL[i] )
			return i;
	return -1;
}

/** iso_setRedisFun assign the line redistribution function type
 * \param ipISO isoelectronic sequence
 * \param nelem element index
 * \param ipLo  index to lower state
 * \param ipHi  index to upper state
 */
void iso_setRedisFun (long ipISO, long nelem, long ipLo, long ipHi)
{
	if( ipLo == 0 && ipHi == iso_ctrl.nLyaLevel[ipISO] )
	{
		long redis = iso_ctrl.ipLyaRedist[ipISO];
		// H LyA has a special redistribution function 
		if( ipISO==ipH_LIKE && nelem==ipHYDROGEN )
			redis = ipLY_A;
		iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().iRedisFun() = redis;
	}
	else if( ipLo == 0 )
	{
		/* these are rest of Lyman lines, 
		 * complete redistribution, doppler core only, K2 core, default ipCRD */
		iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().iRedisFun() = iso_ctrl.ipResoRedist[ipISO];
	}
	else
	{
		/* all lines coming from excited states, default is complete
		 * redis with wings, ipCRDW*/
		iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().iRedisFun() = iso_ctrl.ipSubRedist[ipISO];
	}

	return;
}



/** iso_setOpacity compute line opacity
 * \param ipISO isoelectronic sequence
 * \param nelem element index
 * \param ipLo  index to lower state
 * \param ipHi  index to upper state
 */
void iso_setOpacity (long ipISO, long nelem, long ipLo, long ipHi)
{
	if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA ||
		iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() <= 0.)
	{
		iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().gf() = 0.;
		iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().opacity() = 0.;
	}
	else
	{
		iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().gf() = 
			(realnum)(GetGF(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul(),
			iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN(),
			iso_sp[ipISO][nelem].st[ipHi].g()));
		ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().gf() > 0.);

		/* derive the abs coef, call to function is gf, wl (A), g_low */
		iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().opacity() = 
			(realnum)(abscf(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().gf(),
			iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN(),
			iso_sp[ipISO][nelem].st[ipLo].g()));
		ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().opacity() > 0.);
	}

	return;
}

#if 0
// hydro_energy: get the energy needed to ionize level (in cm^-1)
// the parameters l, s, and g are currently not used, but will be in the future
double hydro_energy(long nelem, long n, long /* l */, long /* s  */, long /* g */)
{
	DEBUG_ENTRY( "hydro_energy()" );

	double HIonPoten;
	
	/* the level energies and line wavelengths are derived from these
	 * Rydberg constants.  Changing the constants will change the wavelength,
	 * but beware that many H I, He I, and He II wavelengths are hardwired into
	 * the code.  These include lists of Case B line predictions and places
	 * where lines are retrieved.  If the Rydberg constants are changed it would
	 * be wise to invest the time to derive all wavelengths self consistently
	 * with global variables and remove the explicit lists.  See the source
	 * files changed in r11678 for an idea of where these lines occur.
	 */
	if( nelem==ipHYDROGEN )
		/* this is the NIST value for H itself */
		HIonPoten = 0.99946650834637;
	else if( nelem==ipHELIUM )
		HIonPoten = 3.9996319917;
	else
		/* Dima's data in phfit.dat have ionization potentials in eV
		 * with four significant figures*/
		HIonPoten = atmdat.getIonPot(nelem, nelem)/EVRYD;
	ASSERT(HIonPoten > 0.);

	return HIonPoten/POW2((double)n)*RYD_INF;
}
#endif


/* hydro_energy calculates energy of a given level in cm^-1. */
double hydro_energy(long nelem, long n, long l, long s, long j)
{
        DEBUG_ENTRY( "hydro_energy()" );

	enum{ DEBUG_LOC = false };


	/* NB NB
	 * 2021-06-02, M.Chatzikos
	 *
	 * Do not use logic based on the number of resolved levels in the model,
	 * as it leads to an ASSERT() being thrown when the command
	 *    'compile recombination coefficients H-like'
	 * is used.  That's b/c that command expects more resolved levels than
	 * typically available in the Stout .nrg.dat files.
	 */

	double E_H = iso_sp[ipH_LIKE][nelem].energy_ioniz(n, l, s, 2*j+1);

	if( E_H < 0. )
	{
		E_H = H_RYD_FACTOR * RYD_INF * POW2((double)(nelem+1)/(double)n);
	}

	if( DEBUG_LOC )
	{
        	fprintf( ioQQQ, "debug E_H n l s j QN2ind IonPot"
				" %.2e %ld %ld %ld %ld %ld %.2e\n",
				E_H, n, l, s, j, QN2ind(n, l, s, 2*j+1),
				iso_sp[ipH_LIKE][nelem].IonPot );
	}

        ASSERT(E_H > 0.);

        return E_H;
}


void iso_create()
{
	long int ipHi, 
		ipLo;

	static int nCalled = 0;

	DEBUG_ENTRY( "iso_create()" );

	/* > 1 if not first call, then just zero arrays out */
	if( nCalled > 0 )
	{
		iso_zero();
		return;
	}

	/* this is first call, increment the nCalled counterso never do this again */
	++nCalled;

	/* these are the statistical weights of the ions */
	iso_ctrl.stat_ion[ipH_LIKE] = 1.f;
	iso_ctrl.stat_ion[ipHE_LIKE] = 2.f;

	/* this routine allocates all the memory
	 * needed for the iso sequence structures */
	iso_allocate();

	/* loop over iso sequences and assign quantum numbers to all levels */
	iso_assign_quantum_numbers();

	/* initialize extra lyman lines */
	iso_assign_extralyman_levels();

	/* this is a dummy line, junk it too. */
	(*TauDummy).Junk();
	(*TauDummy).AddHiState();
	(*TauDummy).AddLoState();
	(*TauDummy).AddLine2Stack();

	/********************************************/
	/**********  Line and level energies ********/
	/********************************************/
	iso_init_energies();
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* main hydrogenic arrays, fill with sane values */
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				double EnergyRydGround = 0.;
				/* go from ground to the highest level */
				for( ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					double EnergyWN, EnergyRyd;

					if( ipISO == ipH_LIKE )
					{
						/* NB NB
						 *
						 * 2020-08-16
						 *
						 * The following assumes that the data are resolved in the
						 * input data file, but the stored values of l, s, j are -1
						 * for each, and the following call crashes.
						 *
						 * The hack below is temporary until we resolve the data in
						 * the data file.
						 */
						// EnergyRyd = hydro_energy(nelem, N_(ipHi), L_(ipHi),
						//					S_(ipHi), J_(ipHi)) * WAVNRYD;
						EnergyRyd = hydro_energy(nelem, N_(ipHi), -1, -1, -1) * WAVNRYD;
					}
					else if( ipISO == ipHE_LIKE )
					{
						EnergyRyd = helike_energy(nelem, N_(ipHi), L_(ipHi),
									  S_(ipHi), J_(ipHi)) * WAVNRYD;
					}
					else
					{
						/* Other iso sequences don't exist yet. */
						TotalInsanity();
					}

					/* >>chng 02 feb 09, change test to >= 0 since we now use 0 for 2s-2p */
					ASSERT( EnergyRyd >= 0. );

					iso_sp[ipISO][nelem].fb[ipHi].xIsoLevNIonRyd = EnergyRyd;
					if (ipHi == 0)
						EnergyRydGround = EnergyRyd;
					iso_sp[ipISO][nelem].st[ipHi].energy().set(EnergyRydGround-EnergyRyd);

					/* now loop from ground to level below ipHi */
					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						if( false && N_(ipLo) == N_(ipHi) && ipISO == ipH_LIKE &&
							L_(ipHi) == L_(ipLo)+1 && size_t(N_(ipHi)) <= t_hydro_tbl::Inst().nmaxnn() )
							/** \todo 2 this should support j-resolved transitions as wavl depends strongly on that */
							EnergyWN = t_hydro_tbl::Inst().wn( N_(ipLo), L_(ipLo), L_(ipHi), nelem+1 );
						else
							EnergyWN = RYD_INF * fabs(iso_sp[ipISO][nelem].fb[ipLo].xIsoLevNIonRyd -
													  iso_sp[ipISO][nelem].fb[ipHi].xIsoLevNIonRyd);

						if( ipISO == ipH_LIKE && ipHi == 1 )
						{
							EnergyWN = t_hydro_tbl::Inst().m1wn(N_(ipHi), L_(ipHi), nelem+1);
						}

						/* transition energy in various units: */
						iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() = (realnum)EnergyWN;

						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() >= 0.);
						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyErg() >= 0.);
						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyK() >= 0.);

						if( rfield.isEnergyBound(Energy(EnergyWN, "cm^-1")) )
						{
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLangVac() = 
								wn2angVac( double( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() ) );
							ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLangVac() > 0.);
						}
						else
						{
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLangVac() = 1.e30_r;
						}
					}
				}

				if( ipISO == ipH_LIKE)
				{
					/* fill the extra Lyman lines */
					for( long nHi=2; nHi < iso_ctrl.nLymanHLike[nelem]; nHi++ )
					{
						FillExtraLymanLine( ExtraLymanLinesJ05[nelem].begin()+ipExtraLymanLinesJ05[nelem][nHi], ipISO, nelem, nHi, 0.5 );
						FillExtraLymanLine( ExtraLymanLinesJ15[nelem].begin()+ipExtraLymanLinesJ15[nelem][nHi], ipISO, nelem, nHi, 1.5 );

						enum {DEBUG_LOC=false};
						if(DEBUG_LOC && ipISO == ipH_LIKE && nelem == ipIRON && nHi == 2)
						{
							fprintf( ioQQQ, "%li\t%li\t%f\n",
										nelem,
										nHi,
										ExtraLymanLinesJ05[nelem][ipExtraLymanLinesJ05[nelem][nHi]].EnergyWN()
										);

							cdEXIT(EXIT_FAILURE);
						}
					}
				}
				else if( ipISO == ipHE_LIKE )
				{
					/* fill the extra Lyman lines */
					for( long ipHi=2; ipHi < iso_ctrl.nLyman_alloc[ipISO]; ipHi++ )
					{
						FillExtraLymanLine( ExtraLymanLinesHeLike[nelem].begin()+ipExtraLymanLinesHeLike[nelem][ipHi], ipISO, nelem, ipHi, -1. );
					}
				}
				else
					TotalInsanity();
			}
		}
	}

	/***************************************************************/
	/***** Set up recombination tables for later interpolation *****/
	/***************************************************************/
	/* NB - the above is all we need if we are compiling recombination tables. */
	iso_recomb_alloc();
	iso_recomb_setup( ipH_LIKE );
	iso_recomb_setup( ipHE_LIKE );
	iso_recomb_auxiliary_free();

	/* set up helium collision tables */
	HeCollidSetup();

	/***********************************************************************************/
	/**********  Transition Probabilities, Redistribution Functions, Opacitites ********/
	/***********************************************************************************/
	enum { DEBUG_EINA_A_NNP = false };

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		if( ipISO == ipH_LIKE )
		{
			/* do nothing here */
		}
		else if( ipISO == ipHE_LIKE )
		{
			/* This routine reads in transition probabilities from a file. */ 
			HelikeTransProbSetup();
		}
		else
		{
			TotalInsanity();
		}

		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
					{
						realnum Aul;

						/* transition prob, EinstA uses current H atom indices */
						if( ipISO == ipH_LIKE )
						{
							Aul = hydro_transprob( nelem, ipHi, ipLo );
						}
						else if( ipISO == ipHE_LIKE )
						{
							Aul = helike_transprob(nelem, ipHi, ipLo);
						}
						else
						{
							TotalInsanity();
						}

						if( Aul <= iso_ctrl.SmallA )
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipEmis() = -1;
						else
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).AddLine2Stack();

						iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() = Aul;

						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() > 0.);

						iso_setRedisFun(ipISO, nelem, ipLo, ipHi);

						iso_setOpacity(ipISO, nelem, ipLo, ipHi);

						if( DEBUG_EINA_A_NNP )
						{
							printf( "%ld\t %ld\t %ld\t %.4e\n",
								nelem+1,
								N_( ipHi ), N_( ipLo ),
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() );
						}
					}
				}
			}
		}
	}

	/****************************************************/
	/**********  lifetimes and damping constants ********/
	/****************************************************/

	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				/* these are not defined and must never be used */
				iso_sp[ipISO][nelem].st[0].lifetime() = -FLT_MAX;

				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					iso_sp[ipISO][nelem].st[ipHi].lifetime() = SMALLFLOAT;
					
					for( ipLo=0; ipLo < ipHi; ipLo++ )  
					{
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
							continue;

						iso_sp[ipISO][nelem].st[ipHi].lifetime() += iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul();
					}

					/* sum of A's was just stuffed, now invert for lifetime. */
					iso_sp[ipISO][nelem].st[ipHi].lifetime() = 1./iso_sp[ipISO][nelem].st[ipHi].lifetime();

					enum { DEBUG_LIFETIMES = false };
					if( DEBUG_LIFETIMES )
					{
						fprintf(ioQQQ , "DEBUG_LIFETIMES\t%ld\t%ld\t %ld\t %.4e\n",
							ipISO , nelem+1,
							N_( ipHi ),
							iso_sp[ipISO][nelem].st[ipHi].lifetime() );
					}

					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() <= 0. )
							continue;

						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
							continue;

						iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().dampXvel() = (realnum)( 
							(1.f/iso_sp[ipISO][nelem].st[ipHi].lifetime())/
							PI4/iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN());

						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().dampXvel()> 0.);
					}
				}
			}
		}
	}

	/* zero out some line information */
	iso_zero();

	/* loop over iso sequences */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( long nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] ) 
			{
				/* calculate cascade probabilities, branching ratios, and associated errors. */
				iso_cascade( ipISO, nelem);
			}
		}
	}

	iso_satellite();

	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
		iso_satellite_update( nelem );

	/***************************************/
	/**********  Stark Broadening **********/
	/***************************************/
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long ipHi= 1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ++ipHi )
				{
					for( long ipLo=0; ipLo < ipHi; ++ipLo )
					{
						iso_sp[ipISO][nelem].ex[ipHi][ipLo].pestrk = 0.;
						iso_sp[ipISO][nelem].ex[ipHi][ipLo].pestrk_up = 0.;
					}
				}

				// signal that model atom has been set up
				long ion = nelem - ipISO;
				strcpy(atmdat.chdBaseSources[nelem][ion], "IsoSeq");
				atmdat.lgdBaseSourceExists[nelem][ion] = true;
			}
		}
	}

	// arrays set up, do not allow number of levels to change in later sims
	lgHydroAlloc = true;

	return;
}


/* ============================================================================== */
STATIC void iso_zero(void)
{
	DEBUG_ENTRY( "iso_zero()" );

	hydro.HLineWidth = 0.;

	/****************************************************/
	/**********  initialize some variables **********/
	/****************************************************/
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					iso_sp[ipISO][nelem].st[ipHi].Pop() = 0.;
					iso_sp[ipISO][nelem].fb[ipHi].Reset();
				}
				if (ipISO <= nelem) 
					iso_sp[ipISO][nelem].st[0].Pop() =  
						dense.xIonDense[nelem][nelem-ipISO]; 
			}
		}
	}

	/* ground state of H and He is different since totally determine
	 * their own opacities */
	iso_sp[ipH_LIKE][ipHYDROGEN].fb[0].ConOpacRatio = 1e-5;
	if( dense.lgElmtOn[ipHELIUM] )
	{
		iso_sp[ipH_LIKE][ipHELIUM].fb[0].ConOpacRatio = 1e-5;
		iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ConOpacRatio = 1e-5;
	}
	return;
}

STATIC void iso_allocate(void)
{

	DEBUG_ENTRY( "iso_allocate()" );

	/* the hydrogen and helium like iso-sequences */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		isoSpecies[ipISO].database = iso_ctrl.chISO[ipISO];

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( dense.lgElmtOn[nelem] )
			{
				t_iso_sp *sp = &iso_sp[ipISO][nelem];
				sp->numLevels_alloc = sp->numLevels_max;

				{
					string chemicalLabel = makeChemical( nelem, nelem-ipISO );

					sp->lgPrtMatrix = false;
					if( chemicalLabel == prt.matrix.species )
						sp->lgPrtMatrix = true;

					sp->lgImgMatrix = false;
					if( chemicalLabel == save.img_matrix.species )
						sp->lgImgMatrix = true;
				}

				ASSERT( sp->numLevels_max > 0 );
				ASSERT( iso_ctrl.nLyman_alloc[ipISO] == iso_ctrl.nLyman[ipISO] );

				sp->ipTrans.reserve( sp->numLevels_max );
				sp->ex.reserve( sp->numLevels_max );
				sp->CascadeProb.reserve( sp->numLevels_max );
				sp->BranchRatio.reserve( sp->numLevels_max );
				//sp->st.resize( sp->numLevels_max );
				sp->fb.resize( sp->numLevels_max );

				for( long n=1; n < sp->numLevels_max; ++n )
				{
					sp->ipTrans.reserve( n, n );
				}

				for( long n=0; n < sp->numLevels_max; ++n )
				{
					sp->ex.reserve( n, sp->numLevels_max );
					sp->CascadeProb.reserve( n, sp->numLevels_max );
					sp->BranchRatio.reserve( n, sp->numLevels_max );
				}

				sp->ipTrans.alloc();
				sp->ex.alloc();
				sp->CascadeProb.alloc();
				sp->BranchRatio.alloc();
			}
		}
	}

	{
		long ipISO = ipHE_LIKE;

		ipExtraLymanLinesHeLike.reserve( LIMELM );

		ExtraLymanLinesHeLike.reserve(LIMELM);

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( dense.lgElmtOn[nelem] )
			{
				ASSERT( iso_sp[ipISO][nelem].numLevels_max > 0 );

				ipExtraLymanLinesHeLike.reserve( nelem, iso_ctrl.nLyman_alloc[ipISO] );
			}
		}

		ipExtraLymanLinesHeLike.alloc();

		for( long nelem=0; nelem < ipISO; ++nelem )
		{
			ExtraLymanLinesHeLike.push_back(
				TransitionList("Insanity",&AnonStates));
		}

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				ExtraLymanLinesHeLike.push_back(
					TransitionList("ExtraLymanLinesHeLike",&iso_sp[ipISO][nelem].st));
			}
			else
			{
				ExtraLymanLinesHeLike.push_back(
					TransitionList("Insanity",&AnonStates));
			}
		}

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( dense.lgElmtOn[nelem] )
			{
				AllTransitions.push_back(ExtraLymanLinesHeLike[nelem]);
				ExtraLymanLinesHeLike[nelem].resize(iso_ctrl.nLyman_alloc[ipISO]-2);
			}
		}
	}

	/* Set up H-like extra lyman lines */
	{
		long ipISO = ipH_LIKE;

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				iso_ctrl.nLymanHLike[nelem] =
					iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max + iso_ctrl.nLyman_alloc[ipISO] + 1;
			}
		}

		ipExtraLymanLinesJ05.reserve( LIMELM );
		ipExtraLymanLinesJ15.reserve( LIMELM );

		ExtraLymanLinesJ05.reserve(LIMELM);
		ExtraLymanLinesJ15.reserve(LIMELM);
	
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				ipExtraLymanLinesJ05.reserve( nelem, iso_ctrl.nLymanHLike[nelem] );
				ipExtraLymanLinesJ15.reserve( nelem, iso_ctrl.nLymanHLike[nelem] );
			}
		}

		ipExtraLymanLinesJ05.alloc();
		ipExtraLymanLinesJ15.alloc();

		for( long nelem=0; nelem < ipISO; ++nelem )
		{
			ExtraLymanLinesJ05.push_back(
				TransitionList("Insanity",&AnonStates));
			ExtraLymanLinesJ15.push_back(
				TransitionList("Insanity",&AnonStates));
		}

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				ExtraLymanLinesJ05.push_back(
					TransitionList("ExtraLymanLinesJ05",&iso_sp[ipISO][nelem].stJ05));
				ExtraLymanLinesJ15.push_back(
					TransitionList("ExtraLymanLinesJ15",&iso_sp[ipISO][nelem].stJ15));
			}
			else
			{
				ExtraLymanLinesJ05.push_back(
					TransitionList("Insanity",&AnonStates));
				ExtraLymanLinesJ15.push_back(
					TransitionList("Insanity",&AnonStates));
			}
		}

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( dense.lgElmtOn[nelem] )
			{
				AllTransitions.push_back(ExtraLymanLinesJ05[nelem]);
				ExtraLymanLinesJ05[nelem].resize(iso_ctrl.nLymanHLike[nelem]);

				AllTransitions.push_back(ExtraLymanLinesJ15[nelem]);
				ExtraLymanLinesJ15[nelem].resize(iso_ctrl.nLymanHLike[nelem]);
			}
		}
	}

	ipSatelliteLines.reserve( NISO );

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		ipSatelliteLines.reserve( ipISO, LIMELM );

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( dense.lgElmtOn[nelem] )
			{
				ASSERT( iso_sp[ipISO][nelem].numLevels_max > 0 );

				ipSatelliteLines.reserve( ipISO, nelem, iso_sp[ipISO][nelem].numLevels_max );
			}
		}
	}

	ipSatelliteLines.alloc();

	Transitions.resize(NISO);
	SatelliteLines.resize(NISO);

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		Transitions[ipISO].reserve(LIMELM);
		SatelliteLines[ipISO].reserve(LIMELM);
		
		for( long nelem=0; nelem < ipISO; ++nelem )
		{
			Transitions[ipISO].push_back(
				TransitionList("Insanity",&AnonStates));
			SatelliteLines[ipISO].push_back(
				TransitionList("Insanity",&AnonStates));
		}
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				Transitions[ipISO].push_back(
					TransitionList("Isosequence",&iso_sp[ipISO][nelem].st));
				SatelliteLines[ipISO].push_back(
					TransitionList("SatelliteLines",&iso_sp[ipISO][nelem].st));
			}
			else
			{
				Transitions[ipISO].push_back(
					TransitionList("Insanity",&AnonStates));
				SatelliteLines[ipISO].push_back(
					TransitionList("Insanity",&AnonStates));
			}
		}

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( dense.lgElmtOn[nelem] )
			{
				if( iso_ctrl.lgDielRecom[ipISO] )
				{
					SatelliteLines[ipISO][nelem].resize( iso_sp[ipISO][nelem].numLevels_max );
					AllTransitions.push_back(SatelliteLines[ipISO][nelem]);
					unsigned int nLine = 0;
					for( long ipLo=0; ipLo<iso_sp[ipISO][nelem].numLevels_max; ipLo++ )
					{
						/* Upper level is continuum, use a generic state
						 * lower level is the same as the index. */
						ipSatelliteLines[ipISO][nelem][ipLo] = nLine;
						SatelliteLines[ipISO][nelem][nLine].Junk();
						long ipHi = iso_sp[ipISO][nelem].numLevels_max;
						SatelliteLines[ipISO][nelem][nLine].setHi(ipHi);
						SatelliteLines[ipISO][nelem][nLine].setLo(ipLo);
						SatelliteLines[ipISO][nelem][nLine].AddLine2Stack();
						++nLine;
					}
					ASSERT(SatelliteLines[ipISO][nelem].size() == nLine);
				}

				//iso_sp[ipISO][nelem].tr.resize( iso_sp[ipISO][nelem].ipTrans.size() );
				//iso_sp[ipISO][nelem].tr.states() = &iso_sp[ipISO][nelem].st;
				Transitions[ipISO][nelem].resize( iso_sp[ipISO][nelem].ipTrans.size() );
				AllTransitions.push_back(Transitions[ipISO][nelem]);
				unsigned int nTransition=0;
				for( long ipHi=1; ipHi<iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					for( long ipLo=0; ipLo < ipHi; ipLo++ )
					{
						/* set ENTIRE array to impossible values, in case of bad pointer */
						iso_sp[ipISO][nelem].ipTrans[ipHi][ipLo] = nTransition;
						Transitions[ipISO][nelem][nTransition].Junk();
						Transitions[ipISO][nelem][nTransition].setHi(ipHi);
						Transitions[ipISO][nelem][nTransition].setLo(ipLo);
						++nTransition;
					}
				}
				ASSERT(Transitions[ipISO][nelem].size() == nTransition);
				iso_sp[ipISO][nelem].tr = &Transitions[ipISO][nelem];
			}
		}
	}

	// associate line and level stacks with species
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				long ion = nelem - ipISO;
				ASSERT( ion >= 0 && ion <= nelem );

				molecule *spmole = findspecies( makeChemical( nelem, ion ).c_str() );
				ASSERT( spmole != null_mole );

				mole.species[ spmole->index ].dbase = &isoSpecies[ipISO];
				mole.species[ spmole->index ].dbase->numLevels_local = iso_sp[ipISO][nelem].numLevels_local;
				mole.species[ spmole->index ].dbase->numLevels_max = iso_sp[ipISO][nelem].numLevels_max;
				mole.species[ spmole->index ].levels = &iso_sp[ipISO][nelem].st;
				mole.species[ spmole->index ].lines = &Transitions[ipISO][nelem];
			}
		}
	}

	return;
}

STATIC void iso_assign_quantum_numbers(void)
{
	long int
	  ipLo,
	  level,
	  i,
	  in,
	  il,
	  is,
	  ij;

	DEBUG_ENTRY( "iso_assign_quantum_numbers()" );

	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		long ipISO = ipH_LIKE;
		/* only check elements that are turned on */
		if( dense.lgElmtOn[nelem] )
		{
			i = 0;

			/* 2 for doublet */
			is = ipDOUBLET;

			/* this loop is over quantum number n */
			for( in = 1L; in <= iso_sp[ipISO][nelem].n_HighestResolved_max; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					iso_sp[ipISO][nelem].st[i].n() = in;
					iso_sp[ipISO][nelem].st[i].S() = is;
					iso_sp[ipISO][nelem].st[i].l() = il;
					iso_sp[ipISO][nelem].st[i].j() = -1._r;
					iso_sp[ipISO][nelem].stJ05[i].n() = in;
					iso_sp[ipISO][nelem].stJ05[i].S() = is;
					iso_sp[ipISO][nelem].stJ05[i].l() = il;
					iso_sp[ipISO][nelem].stJ15[i].n() = in;
					iso_sp[ipISO][nelem].stJ15[i].S() = is;
					iso_sp[ipISO][nelem].stJ15[i].l() = il;
					if(in == 1)
					{
						iso_sp[ipISO][nelem].stJ05[i].j() = 0.5_r;
						/* This is the ground state for np3/2 -> 1s1/2 transition */
						iso_sp[ipISO][nelem].stJ15[i].j() = 0.5_r;
					}
					else if(il == 1)
					{
						iso_sp[ipISO][nelem].stJ05[i].j() = 0.5_r;
						iso_sp[ipISO][nelem].stJ15[i].j() = 1.5_r;
					}
					else
					{
						iso_sp[ipISO][nelem].stJ05[i].j() = -1._r;
						iso_sp[ipISO][nelem].stJ15[i].j() = -1._r;
					}
					++i;
				}
			}
			/* now do the collapsed levels */
			in = iso_sp[ipISO][nelem].n_HighestResolved_max + 1;
			for( level = i; level< iso_sp[ipISO][nelem].numLevels_max; ++level)
			{
				iso_sp[ipISO][nelem].st[level].n() = in;
				iso_sp[ipISO][nelem].st[level].S() = -LONG_MAX;
				iso_sp[ipISO][nelem].st[level].l() = -LONG_MAX;
				iso_sp[ipISO][nelem].st[level].j() = -1._r;
				iso_sp[ipISO][nelem].stJ05[level].n() = in;
				iso_sp[ipISO][nelem].stJ05[level].S() = -LONG_MAX;
				iso_sp[ipISO][nelem].stJ05[level].l() = -LONG_MAX;
				iso_sp[ipISO][nelem].stJ05[level].j() = 0.5_r;
				iso_sp[ipISO][nelem].stJ15[level].n() = in;
				iso_sp[ipISO][nelem].stJ15[level].S() = -LONG_MAX;
				iso_sp[ipISO][nelem].stJ15[level].l() = -LONG_MAX;
				iso_sp[ipISO][nelem].stJ15[level].j() = 1.5_r;
				++in;
			}
			--in;

			/* confirm that we did not overrun the array */
			ASSERT( i <= iso_sp[ipISO][nelem].numLevels_max );

			/* confirm that n is positive and not greater than the max n. */
			ASSERT( (in > 0) && (in < (iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max + 1) ) );
		}
	}

	/* then do he-like */
	for( long nelem=ipHELIUM; nelem < LIMELM; nelem++ )
	{
		long ipISO = ipHE_LIKE;
		/* only check elements that are turned on */
		if( dense.lgElmtOn[nelem] )
		{
			i = 0;

			/* this loop is over quantum number n */
			for( in = 1L; in <= iso_sp[ipISO][nelem].n_HighestResolved_max; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					for( is = 3L; is >= 1L; is -= 2 )
					{
						/* All levels except singlet P follow the ordering scheme:	*/
						/*	lower l's have lower energy	*/
						/* 	triplets have lower energy	*/
						if( (il == 1L) && (is == 1L) )
							continue;
						/* n = 1 has no triplet, of course.	*/
						if( (in == 1L) && (is == 3L) )
							continue;

						/* singlets */
						if( is == 1 )
						{
							iso_sp[ipISO][nelem].st[i].n() = in;
							iso_sp[ipISO][nelem].st[i].S() = is;
							iso_sp[ipISO][nelem].st[i].l() = il;
							/* this is not a typo, J=L for singlets.  */
							iso_sp[ipISO][nelem].st[i].j() = (realnum)il;
							++i;
						}
						/* 2 triplet P is j-resolved */
						else if( (in == 2) && (il == 1) && (is == 3) )
						{
							ij = 0;
							do 
							{
								iso_sp[ipISO][nelem].st[i].n() = in;
								iso_sp[ipISO][nelem].st[i].S() = is;
								iso_sp[ipISO][nelem].st[i].l() = il;
								iso_sp[ipISO][nelem].st[i].j() = (realnum)ij;
								++i;
								++ij;
								/* repeat this for the separate j-levels within 2^3P. */
							}	while ( ij < 3 );
						}
						else
						{
							iso_sp[ipISO][nelem].st[i].n() = in;
							iso_sp[ipISO][nelem].st[i].S() = is;
							iso_sp[ipISO][nelem].st[i].l() = il;
							iso_sp[ipISO][nelem].st[i].j() = -1._r;
							++i;
						}
					}
				}
				/*	Insert singlet P at the end of every sequence for a given n.	*/
				if( in > 1L )
				{
					iso_sp[ipISO][nelem].st[i].n() = in;
					iso_sp[ipISO][nelem].st[i].S() = 1L;
					iso_sp[ipISO][nelem].st[i].l() = 1L;
					iso_sp[ipISO][nelem].st[i].j() = 1._r;
					++i;
				}
			}
			/* now do the collapsed levels */
			in = iso_sp[ipISO][nelem].n_HighestResolved_max + 1;
			for( level = i; level< iso_sp[ipISO][nelem].numLevels_max; ++level)
			{
				iso_sp[ipISO][nelem].st[level].n() = in;
				iso_sp[ipISO][nelem].st[level].S() = -LONG_MAX;
				iso_sp[ipISO][nelem].st[level].l() = -LONG_MAX;
				iso_sp[ipISO][nelem].st[level].j() = -1._r;
				++in;
			}
			--in;

			/* confirm that we did not overrun the array */
			ASSERT( i <= iso_sp[ipISO][nelem].numLevels_max );

			/* confirm that n is positive and not greater than the max n. */
			ASSERT( (in > 0) && (in < (iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max + 1) ) );
		}
	}

	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( ipLo=ipH1s; ipLo < iso_sp[ipISO][nelem].numLevels_max; ipLo++ )
				{
					iso_sp[ipISO][nelem].st[ipLo].nelem() = (int)(nelem+1);
					iso_sp[ipISO][nelem].st[ipLo].IonStg() = (int)(nelem+1-ipISO);

					if( iso_sp[ipISO][nelem].st[ipLo].j() >= 0 )
					{
						iso_sp[ipISO][nelem].st[ipLo].g() = 2.f*iso_sp[ipISO][nelem].st[ipLo].j()+1.f;
					}
					else if( iso_sp[ipISO][nelem].st[ipLo].l() >= 0 )
					{
						iso_sp[ipISO][nelem].st[ipLo].g() = (2.f*iso_sp[ipISO][nelem].st[ipLo].l()+1.f) *
							iso_sp[ipISO][nelem].st[ipLo].S();
					}
					else
					{
						if( ipISO == ipH_LIKE )
							iso_sp[ipISO][nelem].st[ipLo].g() = 2.f*(realnum)POW2( iso_sp[ipISO][nelem].st[ipLo].n() );
						else if( ipISO == ipHE_LIKE )
							iso_sp[ipISO][nelem].st[ipLo].g() = 4.f*(realnum)POW2( iso_sp[ipISO][nelem].st[ipLo].n() );
						else
						{
							/* replace this with correct thing if more sequences are added. */
							TotalInsanity();
						}
					}
					
					char chConfiguration[32];
					long nCharactersWritten = 0;

					ASSERT( iso_sp[ipISO][nelem].st[ipLo].n() < 1000 );

					/* Treat H-like levels as collapsed, for the purposes of
					 * reporting comments in the output of 'save lines labels'.
					 * For the rest, include J only if defined, and for singlets if positive. */
					if( ( ipISO == ipH_LIKE && ipLo > ipH2p ) ||
						iso_sp[ipISO][nelem].st[ipLo].n() > iso_sp[ipISO][nelem].n_HighestResolved_max )
					{
						nCharactersWritten = sprintf( chConfiguration, "n=%3li", 
							iso_sp[ipISO][nelem].st[ipLo].n() );
					}
					else if( ( iso_sp[ipISO][nelem].st[ipLo].j() > 0 &&
						iso_sp[ipISO][nelem].st[ipLo].S() == ipSINGLET ) ||
						( iso_sp[ipISO][nelem].st[ipLo].j() >= 0 &&
						iso_sp[ipISO][nelem].st[ipLo].S() == ipTRIPLET ) )
					{
						nCharactersWritten = sprintf( chConfiguration, "%3li^%li%c_%.1f", 
							iso_sp[ipISO][nelem].st[ipLo].n(), 
							iso_sp[ipISO][nelem].st[ipLo].S(),
							chL[ MIN2( 20, iso_sp[ipISO][nelem].st[ipLo].l() ) ],
							iso_sp[ipISO][nelem].st[ipLo].j() );
					}
					else
					{
						nCharactersWritten = sprintf( chConfiguration, "%3li^%li%c", 
							iso_sp[ipISO][nelem].st[ipLo].n(), 
							iso_sp[ipISO][nelem].st[ipLo].S(),
							chL[ MIN2( 20, iso_sp[ipISO][nelem].st[ipLo].l()) ] );
					}

					if( nCharactersWritten < 0 || nCharactersWritten > 31 )
						TotalInsanity();

					iso_sp[ipISO][nelem].st[ipLo].chConfig() = chConfiguration;
				}
			}
		}
	}

	long ipISO = ipH_LIKE;
	for( long nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			for( ipLo=ipH1s; ipLo < iso_sp[ipISO][nelem].numLevels_max; ipLo++ )
			{
				iso_sp[ipISO][nelem].stJ05[ipLo].nelem() = (int)(nelem+1);
				iso_sp[ipISO][nelem].stJ05[ipLo].IonStg() = (int)(nelem+1-ipISO);
				iso_sp[ipISO][nelem].stJ15[ipLo].nelem() = (int)(nelem+1);
				iso_sp[ipISO][nelem].stJ15[ipLo].IonStg() = (int)(nelem+1-ipISO);

				/* stack of states for extra lyman lines
				   starting with j=1/2 */
				if( iso_sp[ipISO][nelem].stJ05[ipLo].j() >= 0 )
				{
					iso_sp[ipISO][nelem].stJ05[ipLo].g() = 2.f*iso_sp[ipISO][nelem].stJ05[ipLo].j()+1.f;
				}
				else if( iso_sp[ipISO][nelem].stJ05[ipLo].l() >= 0 )
				{
					iso_sp[ipISO][nelem].stJ05[ipLo].g() = (2.f*iso_sp[ipISO][nelem].stJ05[ipLo].l()+1.f) *
						iso_sp[ipISO][nelem].stJ05[ipLo].S();
				}
				else
				{
					iso_sp[ipISO][nelem].stJ05[ipLo].g() = 2.f*(realnum)POW2( iso_sp[ipISO][nelem].stJ05[ipLo].n() );
				}

				/* j = 3/2 */
				if( iso_sp[ipISO][nelem].stJ15[ipLo].j() >= 0 )
				{
					iso_sp[ipISO][nelem].stJ15[ipLo].g() = 2.f*iso_sp[ipISO][nelem].stJ15[ipLo].j()+1.f;
				}
				else if( iso_sp[ipISO][nelem].stJ15[ipLo].l() >= 0 )
				{
					iso_sp[ipISO][nelem].stJ15[ipLo].g() = (2.f*iso_sp[ipISO][nelem].stJ15[ipLo].l()+1.f) *
						iso_sp[ipISO][nelem].stJ15[ipLo].S();
				}
				else
				{
					iso_sp[ipISO][nelem].stJ15[ipLo].g() = 2.f*(realnum)POW2( iso_sp[ipISO][nelem].stJ15[ipLo].n() );
				}

				ASSERT( iso_sp[ipISO][nelem].stJ05[ipLo].n() < 1000 );
				ASSERT( iso_sp[ipISO][nelem].stJ15[ipLo].n() < 1000 );
			}
		}
	}
	return;
}

STATIC void iso_assign_extralyman_levels()
{
	DEBUG_ENTRY( "iso_assign_extralyman_levels()" );

	long ipISO = ipH_LIKE;

	for( long nelem=ipISO; nelem < LIMELM; ++nelem )
	{
		/* only grab core for elements that are turned on */
		if( dense.lgElmtOn[nelem] )
		{
			ExtraLymanLinesJ05[nelem].states() = &iso_sp[ipISO][nelem].stJ05;
			/* We will be indexing by principal quantum number, so the first two elements must not be used. */
			ExtraLymanLinesJ05[nelem][0].Junk();
			ExtraLymanLinesJ05[nelem][1].Junk();
			for( long nHi=2; nHi < iso_ctrl.nLymanHLike[nelem]; nHi++ )
			{
				ipExtraLymanLinesJ05[nelem][nHi] = nHi;
				ExtraLymanLinesJ05[nelem][nHi].Junk();

				long ipHi;
				if( nHi <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_max ) /* resolved levels */
				{
					ipHi = iso_sp[ipH_LIKE][nelem].QN2Index(nHi,1,2);
				}
				else if( nHi <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + iso_sp[ipH_LIKE][nelem].nCollapsed_max ) /* collapsed levels */
				{
					ipHi = iso_sp[ipH_LIKE][nelem].QN2Index(nHi,-1,-1);
				}
				else
				{
					ipHi = iso_sp[ipISO][nelem].numLevels_max + nHi - (iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + iso_sp[ipH_LIKE][nelem].nCollapsed_max) - 1;
				}

				ExtraLymanLinesJ05[nelem][nHi].setHi(ipHi);
				/* lower level is just ground state of the ion */
				ExtraLymanLinesJ05[nelem][nHi].setLo(0);
				ExtraLymanLinesJ05[nelem][nHi].AddLine2Stack();
				ASSERT( ExtraLymanLinesJ05[nelem][nHi].Lo()->g() > 0. );
			}

			ExtraLymanLinesJ15[nelem].states() = &iso_sp[ipISO][nelem].stJ15;
			ExtraLymanLinesJ15[nelem][0].Junk();
			ExtraLymanLinesJ15[nelem][1].Junk();
			for( long nHi=2; nHi < iso_ctrl.nLymanHLike[nelem]; nHi++ )
			{
				ipExtraLymanLinesJ15[nelem][nHi] = nHi;
				ExtraLymanLinesJ15[nelem][nHi].Junk();

				long ipHi;
				if( nHi <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_max ) /* resolved levels */
				{
					ipHi = iso_sp[ipH_LIKE][nelem].QN2Index(nHi,1,2);
				}
				else if( nHi <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + iso_sp[ipH_LIKE][nelem].nCollapsed_max ) /* collapsed levels */
				{
					ipHi = iso_sp[ipH_LIKE][nelem].QN2Index(nHi,-1,-1);
				}
				else
				{
					ipHi = iso_sp[ipISO][nelem].numLevels_max + nHi - (iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + iso_sp[ipH_LIKE][nelem].nCollapsed_max) - 1;
				}

				ExtraLymanLinesJ15[nelem][nHi].setHi(ipHi);
				/* lower level is just ground state of the ion */
				ExtraLymanLinesJ15[nelem][nHi].setLo(0);
				ExtraLymanLinesJ15[nelem][nHi].AddLine2Stack();
				ASSERT( ExtraLymanLinesJ15[nelem][nHi].Lo()->g() > 0. );
			}
		}
	}

	ipISO = ipHE_LIKE;

	for( long nelem=ipISO; nelem < LIMELM; ++nelem )
	{
		/* only grab core for elements that are turned on */
		if( dense.lgElmtOn[nelem] )
		{
			ExtraLymanLinesHeLike[nelem].states() = &iso_sp[ipISO][nelem].st;
			unsigned int nExtraLyman = 0;
			for( long ipHi=2; ipHi < iso_ctrl.nLyman_alloc[ipISO]; ipHi++ )
			{
				ipExtraLymanLinesHeLike[nelem][ipHi] = nExtraLyman;
				ExtraLymanLinesHeLike[nelem][nExtraLyman].Junk();
				long ipHi_offset = iso_sp[ipISO][nelem].numLevels_max + ipHi - 2;
				if( iso_ctrl.lgDielRecom[ipISO] )
					ipHi_offset += 1;
				ExtraLymanLinesHeLike[nelem][nExtraLyman].setHi(ipHi_offset);
				/* lower level is just ground state of the ion */
				ExtraLymanLinesHeLike[nelem][nExtraLyman].setLo(0);
				ExtraLymanLinesHeLike[nelem][nExtraLyman].AddLine2Stack();
				++nExtraLyman;
			}
			ASSERT(ExtraLymanLinesHeLike[nelem].size() == nExtraLyman);
		}
	}
}

#if defined(__ICC) && defined(__i386)
#pragma optimization_level 1
#endif
STATIC void FillExtraLymanLine( const TransitionList::iterator& t, long ipISO, long nelem, long nHi, double j )
{
	double Enerwn, Aul;

	DEBUG_ENTRY( "FillExtraLymanLine()" );

	(*(*t).Hi()).status() = LEVEL_INACTIVE;

	/* atomic number or charge and stage: */
	(*(*t).Hi()).nelem() = (int)(nelem+1);
	(*(*t).Hi()).IonStg() = (int)(nelem+1-ipISO);
	
	(*(*t).Hi()).n() = nHi;

	/* statistical weight is same as statistical weight of corresponding LyA. */
	if (ipISO == ipH_LIKE)
		(*(*t).Hi()).g() = 2.*j+1.;
	else
		(*(*t).Hi()).g() = iso_sp[ipISO][nelem].st[iso_ctrl.nLyaLevel[ipISO]].g();

	/* \todo add correct configuration, or better still link to standard level */
	if( ipISO == ipH_LIKE )
	{
		char chConfiguration[32];
		int j_numerator = j*2;
		sprintf( chConfiguration, "%li^2P%i/2", nHi, j_numerator );
		(*(*t).Lo()).chConfig() = "1^2S1/2";
		(*(*t).Hi()).chConfig() = chConfiguration;
	}
	else
		(*(*t).Hi()).chConfig()  = "ExtraLyman level (probably duplicate)";

	/* energies */
	/********************************************************************************/
	/* Adding Fine Structure Corrections to nP energy levels: 			*/
	/* For all one-electron atoms the following equations     			*/
	/* are being used to find the j-resolved energies.        			*/
	/*										*/
	/*										*/
	/* Total nP_j level energy:                               			*/
	/* Enerwn = E_0_n + E_FS_nj + E_LS_nlj + E_M_nj           			*/
	/*										*/
	/*										*/
	/* Unperturbed Energy:								*/
	/*         -µ Z^2 Ry								*/
	/* E_0_n = ----------								*/
	/*            n^2 								*/
	/*										*/
	/* Dirac-point nucleas energy (Fine Structure correction): 			*/
	/*                                     α Z               			*/
	/* E_FS_nj = m_e c^2 [ 1 + ( ------------------------- )^2 ]^-1/2 - m_e c^2   	*/
	/*                           n - k + √ (k^2 - α^2 Z^2)     			*/
	/*										*/
	/* Lamb-Shift Correction (Resolving orbital-angular momentum l):		*/
	/*			   8 Z^4 α^3               Z^2 Ry        3    c_lj   	*/
	/* E_LS_nlj = ----------- Ry [log ( ----------- ) + --- --------]  		*/
	/*			    3 π n^3   		      K_0(n,l)       8   2l + 1	*/
	/*										*/
	/*	where  c_lj = -l^{-1} if j = l-1/2; c_lj = (l+1)^{-1} if j=l+1/2	*/
	/*	       K_0(n,l) = Bethe logarithm; we have used a fit function instead	*/
	/*										*/
	/* Nuclear Mass Recoil Correction:						*/
	/*		     m_e (α Z)^2	   m_e	 (α Z)^2 			*/
	/* E_M_nj = m_e c^2 ------------ - µ c^2 (-----) --------			*/
	/*		      m_N 2 N^2 	   m_N	   2n^2				*/
	/*										*/
	/* where  N = (( n - k + √ (k^2 - α^2 Z^2) )^2 + α^2 Z^2)^1/2 			*/
	/*										*/
	/* m_e = electron mass 								*/
	/* m_N = nuclear mass 								*/
	/* µ = Reduced mass 								*/
	/* c = speed of light 								*/
	/* α = fine structure constant 							*/
	/* Z = atomic number 								*/
	/* n = principle quantum number 						*/
	/* l = orbital angular quantum number 						*/
	/* j = total quantum number (l+s) 						*/
	/* k = j + 1/2 									*/
	/* Ry = hc R∞ = Rydberg unit of energy 						*/
	/* 										*/
	/* See Gunasekera et. al 2023 for further detail.				*/
	/* 										*/
	/********************************************************************************/
	if( ipISO == ipH_LIKE )
	{
		double reduced_mass = (nelem+1)*PROTON_MASS*ELECTRON_MASS/(((nelem+1)*PROTON_MASS)+ELECTRON_MASS);
		double ELECTRON_REST_ENGY_cm = 4.1214844832e9;
		double mc2 = ELECTRON_REST_ENGY_cm*reduced_mass/ELECTRON_MASS;
		double Eo = -iso_sp[ipISO][nelem].IonPot;
		double Za = FINE_STRUCTURE*(nelem+1.);

		/*Resolving j levels: Dirac equation*/
		double DiracN = nHi - (j + 0.5) + pow( pow2(j + 0.5) - pow2(Za), 0.5);
		double DiracEnergy = ELECTRON_REST_ENGY_cm*(pow(1 + pow2(Za/DiracN), -0.5)  - 1.) - Eo;

		/* Mass Recoil Correction */
		double N = sqrt( pow2( nHi - (j + 0.5) + sqrt( pow2(j + 0.5) - pow2(Za)) ) + pow2(Za) );
		double mass_recoil_correction_term1 = ELECTRON_REST_ENGY_cm*(ELECTRON_MASS/(PROTON_MASS*(nelem+1.)))*pow2(Za)/(2*pow2(N));
		double mass_recoil_correction = mass_recoil_correction_term1 - mc2*pow2(ELECTRON_MASS/(PROTON_MASS*(nelem+1.)))*pow2(Za)/(2*pow2(nHi));

		/* Resolving l levels: Lamb Shift */
		double m=0.04950538, t=0.54290526, b=0.95342613;
		double Bethe_Log = ( m*exp(-nHi*t) + b );
		Bethe_Log = Bethe_Log * exp( pow(1./nHi,3.) * pow((nelem)/(nelem+1.),4.) * log(Bethe_Log) );
		double Lamb_correction_factor = 8.*pow(nelem+1.,4.)*pow(FINE_STRUCTURE,3.)*RYD_INF/(3.*M_PI*pow(nHi,3.));
		double Lamb_correction = Lamb_correction_factor*log(1./Bethe_Log);
		if ( j==0.5 )
		{
			Lamb_correction += Lamb_correction_factor*(-1./3.)*(3./8.);
		}
		else if ( j==1.5 )
		{
			Lamb_correction += Lamb_correction_factor*(-1./6.)*(3./8.);
		}

		if (nelem > 0)
		{
			Enerwn = DiracEnergy + (mass_recoil_correction*0.5) + Lamb_correction;
		}
		else
		{
			Enerwn = DiracEnergy + mass_recoil_correction + Lamb_correction;
		}
	}
	else
	{
		Enerwn = iso_sp[ipISO][nelem].fb[0].xIsoLevNIonRyd * RYD_INF * (  1. - 1./POW2((double)nHi) );
	}

	/* transition energy in various units:*/
	(*t).EnergyWN() = (realnum)(Enerwn);
	(*t).WLangVac() = wn2angVac( Enerwn );

	(*(*t).Hi()).energy().set( Enerwn, "cm^-1" );

	if( ipISO == ipH_LIKE )
	{
		Aul = H_Einstein_A( nHi, 1, 1, 0, nelem+1 );
	}
	else
	{
		if( nelem == ipHELIUM )
		{
			/* A simple fit for the calculation of Helium lyman Aul's.	*/
			Aul = (1.508e10) / pow((double)nHi,2.975);
		}
		else 
		{
			/* Fit to values given in 
			 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
			 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
			/* originally astro.ph. 0201454  */
			Aul = 1.375E10 * pow((double)nelem, 3.9) / pow((double)nHi,3.1);
		}
	}

	(*t).Emis().Aul() = (realnum)Aul;

	(*(*t).Hi()).lifetime() = iso_state_lifetime( ipISO, nelem, nHi, 1 );

	(*t).Emis().dampXvel() = (realnum)( 1.f / (*(*t).Hi()).lifetime() / PI4 / (*t).EnergyWN() );

	if(ipISO == ipH_LIKE)
	{
		long ipHi;
		if( nHi > iso_sp[ipH_LIKE][nelem].n_HighestResolved_local ) /* collapsed levels */
		{
			ipHi = iso_sp[ipH_LIKE][nelem].QN2Index(nHi,-1,-1);
		}
		else /* resolved levels */
		{
			ipHi = iso_sp[ipH_LIKE][nelem].QN2Index(nHi,1,2);
		}

		if( ipHi == iso_ctrl.nLyaLevel[ipISO] )
		{
			long redis = iso_ctrl.ipLyaRedist[ipISO];
			// H LyA has a special redistribution function
			if( ipISO==ipH_LIKE && nelem==ipHYDROGEN )
				redis = ipLY_A;
			(*t).Emis().iRedisFun() = redis;
		}
		else
		{
			/* these are rest of Lyman lines,
			 * complete redistribution, doppler core only, K2 core, default ipCRD */
			(*t).Emis().iRedisFun() = iso_ctrl.ipResoRedist[ipISO];
		}
	}
	else
		(*t).Emis().iRedisFun() = iso_ctrl.ipResoRedist[ipISO];

	(*t).Emis().gf() = (realnum)(GetGF((*t).Emis().Aul(),	(*t).EnergyWN(), (*(*t).Hi()).g()));

	/* derive the abs coef, call to function is Emis().gf(), wl (A), g_low */
	(*t).Emis().opacity() = (realnum)(abscf((*t).Emis().gf(), (*t).EnergyWN(), (*(*t).Lo()).g()));

	/* create array indices that will blow up */
	(*t).ipCont() = INT_MIN;
	(*t).Emis().ipFine() = INT_MIN;

	{
		/* option to print particulars of some line when called
			* a prettier print statement is near where chSpin is defined below
			* search for "pretty print" */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"%li\t%li\t%.2e\t%.2e\n",
				nelem+1,
				nHi,
				(*t).Emis().Aul() , 	
				(*t).Emis().opacity()
				);
		}
	}
	return;
}

/* calculate radiative lifetime of an individual iso state */
double iso_state_lifetime( long ipISO, long nelem, long n, long l )
{
	/* >>refer hydro lifetimes	Horbatsch, M. W., Horbatsch, M. and Hessels, E. A. 2005, JPhysB, 38, 1765 */

	double tau, t0, eps2;
	/* mass of electron */
	double m = ELECTRON_MASS;
	/* nuclear mass */
	double M = (double)dense.AtomicWeight[nelem] * ATOMIC_MASS_UNIT;
	double mu = (m*M)/(M+m);
	long z = 1;
	long Z = nelem + 1 - ipISO;
	
	DEBUG_ENTRY( "iso_state_lifetime()" );

	/* this should not be used for l=0 per the Horbatsch et al. paper */
	ASSERT( l > 0 );

	eps2 = 1. - ( l*l + l + 8./47. - (l+1.)/69./n ) / POW2( (double)n );

	t0 = 3. * H_BAR * powi( (double)n, 5 ) / 
		( 2. * pow4( (double)( z * Z ) ) * powi( FINE_STRUCTURE, 5 ) * mu * pow2( SPEEDLIGHT ) ) *
		pow2( (m + M)/(Z*m + z*M) );

	tau = t0 * ( 1. - eps2 ) / 
		( 1. + 19./88.*( (1./eps2 - 1.) * log( 1. - eps2 ) + 1. - 
		0.5 * eps2 - 0.025 * eps2 * eps2 ) );

	if( ipISO == ipHE_LIKE )
	{
		/* iso_state_lifetime is not spin specific, must exclude helike triplet here. */
		tau /= 3.;
		/* this is also necessary to correct the helike lifetimes */
		tau *= 1.1722 * pow( (double)nelem, 0.1 );
	}

	/* would probably need a new lifetime algorithm for any other iso sequences. */
	ASSERT( ipISO <= ipHE_LIKE );
	ASSERT( tau > 0. );

	return tau;
}

/* calculate cascade probabilities, branching ratios, and associated errors. */
void iso_cascade( long ipISO, long nelem )
{
	long int i, j, ipLo, ipHi;

	DEBUG_ENTRY( "iso_cascade()" );

	/* The sum of all A's coming out of a given n,
	 * Below we assert a monotonic trend. */
	vector<double> SumAPerN(iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max + 1);

	/* Initialize some ground state stuff, easier here than in loops.	*/
	iso_sp[ipISO][nelem].CascadeProb[0][0] = 1.;
	if( iso_ctrl.lgRandErrGen[ipISO] )
	{
		iso_sp[ipISO][nelem].fb[0].SigmaAtot = 0.;
		iso_sp[ipISO][nelem].ex[0][0].SigmaCascadeProb = 0.;
	}

	/***************************************************************************/
	/****** Cascade probabilities, Branching ratios, and associated errors *****/
	/***************************************************************************/
	for( ipHi=1; ipHi<iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
	{
		double SumAs = 0.;

		/** Cascade probabilities are as defined in Robbins 68,
		 * generalized here for cascade probability for any iso sequence.	
		 * >>refer He triplets	Robbins, R.R. 1968, ApJ 151, 497R	
		 * >>refer He triplets	Robbins, R.R. 1968a, ApJ 151, 511R	*/

		/* initialize variables. */
		iso_sp[ipISO][nelem].CascadeProb[ipHi][ipHi] = 1.;
		iso_sp[ipISO][nelem].CascadeProb[ipHi][0] = 0.;
		iso_sp[ipISO][nelem].BranchRatio[ipHi][0] = 0.;

		if( iso_ctrl.lgRandErrGen[ipISO] )
		{
			iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot = 0.;
			iso_sp[ipISO][nelem].ex[ipHi][ipHi].SigmaCascadeProb = 0.;
		}

		long ipLoStart = 0;
		if( opac.lgCaseB && L_(ipHi)==1 && (ipISO==ipH_LIKE || S_(ipHi)==1) )
			ipLoStart = 1;

		for( ipLo=ipLoStart; ipLo<ipHi; ipLo++ )
		{
			SumAs += iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul();
		}

		for( ipLo=ipLoStart; ipLo<ipHi; ipLo++ )
		{
			if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
			{
				iso_sp[ipISO][nelem].CascadeProb[ipHi][ipLo] = 0.;
				iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] = 0.;
				continue;
			}

			iso_sp[ipISO][nelem].CascadeProb[ipHi][ipLo] = 0.;
			iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] = 
				iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() / SumAs;

			ASSERT( iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] <= 1.0000001 );

			SumAPerN[N_(ipHi)] += iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul();

			/* there are some negative energy transitions, where the order
			 * has changed, but these are not optically allowed, these are
			 * same n, different L, forbidden transitions */
			ASSERT( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() > 0. ||
				iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA );

			if( iso_ctrl.lgRandErrGen[ipISO] )
			{
				ASSERT( iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[IPRAD] >= 0. );
				/* Uncertainties in A's are added in quadrature, square root is taken below. */
				iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot += 
					pow2( iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[IPRAD] * 
					(double)iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() );
			}
		}

		if( iso_ctrl.lgRandErrGen[ipISO] )
		{
			/* Uncertainties in A's are added in quadrature above, square root taken here. */
			iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot = sqrt( iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot );
		}

		/* cascade probabilities */
		for( i=0; i<ipHi; i++ )
		{
			for( ipLo=0; ipLo<=i; ipLo++ )
			{
				iso_sp[ipISO][nelem].CascadeProb[ipHi][ipLo] += iso_sp[ipISO][nelem].BranchRatio[ipHi][i] * iso_sp[ipISO][nelem].CascadeProb[i][ipLo];
			}
		}

		if( iso_ctrl.lgRandErrGen[ipISO] )
		{
			for( ipLo=0; ipLo<ipHi; ipLo++ )
			{
				double SigmaCul = 0.;
				for( i=ipLo; i<ipHi; i++ )
				{
					if( iso_sp[ipISO][nelem].trans(ipHi,i).Emis().Aul() > iso_ctrl.SmallA )
					{
						/* Uncertainties in A's and cascade probabilities */
						double SigmaA = iso_sp[ipISO][nelem].ex[ipHi][i].Error[IPRAD] * 
							iso_sp[ipISO][nelem].trans(ipHi,i).Emis().Aul();
						SigmaCul += 
							pow2(SigmaA*iso_sp[ipISO][nelem].CascadeProb[i][ipLo]*iso_sp[ipISO][nelem].st[ipHi].lifetime()) +
							pow2(iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot*iso_sp[ipISO][nelem].BranchRatio[ipHi][i]*
							iso_sp[ipISO][nelem].CascadeProb[i][ipLo]*iso_sp[ipISO][nelem].st[ipHi].lifetime()) +
							pow2(iso_sp[ipISO][nelem].ex[i][ipLo].SigmaCascadeProb*iso_sp[ipISO][nelem].BranchRatio[ipHi][i]);
					}
				}
				SigmaCul = sqrt(SigmaCul);
				iso_sp[ipISO][nelem].ex[ipHi][ipLo].SigmaCascadeProb = SigmaCul;
			}
		}
	}

	/************************************************************************/
	/*** Allowed decay conversion probabilities. See Robbins68b, Table 1. ***/
	/************************************************************************/
	{
		enum {DEBUG_LOC=false};

		if( DEBUG_LOC && (nelem == ipHELIUM) && (ipISO==ipHE_LIKE) )
		{
			/* To output Bm(n,l; ipLo), set ipLo, hi_l, and hi_s accordingly.	*/
			long int hi_l,hi_s;
			double Bm;

			/* these must be set for following output to make sense
			 * as is, a dangerous bit of code - set NaN for safety */
			hi_s = -100000;
			hi_l = -100000;
			ipLo = -100000;
			/* tripS to 2^3P	*/
			//hi_l = 0, hi_s = 3, ipLo = ipHe2p3P0;

			/* tripD to 2^3P	*/
			//hi_l = 2, hi_s = 3, ipLo = ipHe2p3P0;

			/* tripP to 2^3S	*/
			//hi_l = 1, hi_s = 3, ipLo = ipHe2s3S;	

			ASSERT( hi_l != iso_sp[ipISO][nelem].st[ipLo].l() );

			fprintf(ioQQQ,"Bm(n,%ld,%ld;%ld)\n",hi_l,hi_s,ipLo);
			fprintf(ioQQQ,"m\t2\t\t3\t\t4\t\t5\t\t6\n");

			for( ipHi=ipHe2p3P2; ipHi<iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max; ipHi++ )
			{
				/* Pick out excitations from metastable 2tripS to ntripP.	*/
				if( (iso_sp[ipISO][nelem].st[ipHi].l() == 1) && (iso_sp[ipISO][nelem].st[ipHi].S() == 3) )
				{
					fprintf(ioQQQ,"\n%ld\t",iso_sp[ipISO][nelem].st[ipHi].n());
					j = 0;
					Bm = 0;
					for( i = ipLo; i<=ipHi; i++)
					{
						if( (iso_sp[ipISO][nelem].st[i].l() == hi_l) && (iso_sp[ipISO][nelem].st[i].S() == hi_s)  )
						{
							if( (ipLo == ipHe2p3P0) && (i > ipHe2p3P2) )
							{
								Bm += iso_sp[ipISO][nelem].CascadeProb[ipHi][i] * ( iso_sp[ipISO][nelem].BranchRatio[i][ipHe2p3P0] + 
									iso_sp[ipISO][nelem].BranchRatio[i][ipHe2p3P1] + iso_sp[ipISO][nelem].BranchRatio[i][ipHe2p3P2] );
							}
							else
								Bm += iso_sp[ipISO][nelem].CascadeProb[ipHi][i] * iso_sp[ipISO][nelem].BranchRatio[i][ipLo];

							if( (i == ipHe2p3P0) || (i == ipHe2p3P1) || (i == ipHe2p3P2) )
							{
								j++;
								if(j == 3)
								{
									fprintf(ioQQQ,"%2.4e\t",Bm);
									Bm = 0;
								}
							}
							else
							{
								fprintf(ioQQQ,"%2.4e\t",Bm);
								Bm = 0;
							}
						}
					}
				}
			}
			fprintf(ioQQQ,"\n\n");
		}
	}

	/******************************************************/
	/***  Lifetimes should increase monotonically with  ***/
	/***  increasing n...Make sure the A's decrease.    ***/
	/******************************************************/
	for( i=2; i < iso_sp[ipISO][nelem].n_HighestResolved_max; ++i)
	{
		ASSERT( (SumAPerN[i] > SumAPerN[i+1]) || opac.lgCaseB );
	}

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC /* && (ipISO == ipH_LIKE) && (nelem == ipHYDROGEN) */)
		{
			for( i = 2; i<= (iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max); ++i)
			{
				fprintf(ioQQQ,"n %ld\t lifetime %.4e\n", i, 1./SumAPerN[i]);
			}
		}
	}

	return;
}

/** \todo	2	say where these come from	*/	
/* For double-ionization discussions, see Lindsay, Rejoub, & Stebbings 2002	*/
/* Also read Itza-Ortiz, Godunov, Wang, and McGuire 2001.	*/
STATIC void iso_satellite( void )
{
	DEBUG_ENTRY( "iso_satellite()" );

	for( long ipISO = ipHE_LIKE; ipISO < NISO; ipISO++ )
	{
		for( long nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] && iso_ctrl.lgDielRecom[ipISO] ) 
			{
				for( long i=0; i<iso_sp[ipISO][nelem].numLevels_max; i++ )
				{
					TransitionList::iterator tr = SatelliteLines[ipISO][nelem].begin()+ipSatelliteLines[ipISO][nelem][i];
					(*tr).Zero();
			
					/* Make approximation that all levels have energy of H-like 2s level */
					/* Lines to 1s2s have roughly energy of parent Ly-alpha.  So lines to 1snL will have an energy
					 * smaller by the difference between nL and 2s energies.  Therefore, the following has
					 * energy of parent Ly-alpha MINUS the difference between daughter level and daughter n=2 level. */
					(*tr).WLangVac() = (realnum)(RYDLAM/
						((iso_sp[ipISO-1][nelem].fb[0].xIsoLevNIonRyd - iso_sp[ipISO-1][nelem].fb[1].xIsoLevNIonRyd) -
						 (iso_sp[ipISO][nelem].fb[1].xIsoLevNIonRyd- iso_sp[ipISO][nelem].fb[i].xIsoLevNIonRyd)) );

					(*tr).EnergyWN() = 1.e8f / (*tr).WLangVac();

					(*tr).Emis().iRedisFun() = ipCRDW;
					/* this is not the usual nelem, is it atomic not C scale. */
					(*(*tr).Hi()).nelem() = nelem + 1;
					(*(*tr).Hi()).IonStg() = nelem + 1 - ipISO;
					fixit("what should the stat weight of the upper level be? For now say 2.");
					(*(*tr).Hi()).g() = 2.f;

					/* \todo add correct configuration, or better still link to standard level */
					(*(*tr).Hi()).chConfig() = "Satellite level (probably duplicate)";

					// The lower level must already be initialized.  
					ASSERT( (*(*tr).Lo()).g() == iso_sp[ipISO][nelem].st[i].g() );
					//(*(*tr).Lo()).g() = iso_sp[ipISO][nelem].st[i].g();
					(*tr).Emis().PopOpc() = 
						(*(*tr).Lo()).Pop();

					(*tr).Emis().pump() = 0.;

				}
			}
		}
	}

	return;
}

void iso_satellite_update( long nelem )
{
	double ConBoltz, LTE_pop=SMALLFLOAT+FLT_EPSILON, factor1, ConvLTEPOP;
	
	DEBUG_ENTRY( "iso_satellite_update()" );

	for( long ipISO = ipHE_LIKE; ipISO < MIN2(NISO,nelem+1); ipISO++ )
	{
		if( dense.lgElmtOn[nelem] && iso_ctrl.lgDielRecom[ipISO] ) 
		{
			/* This Boltzmann factor is exp( +ioniz energy / Te ).  For simplicity, we make
			 * the fair approximation that all of the autoionizing levels have an energy
			 * equal to the parents n=2. */
			ConBoltz = dsexp(iso_sp[ipISO-1][nelem].fb[1].xIsoLevNIonRyd/phycon.te_ryd);
				
			for( long i=0; i<iso_sp[ipISO][nelem].numLevels_max; i++ )
			{
				double dr_rate = iso_sp[ipISO][nelem].fb[i].DielecRecomb * iso_ctrl.lgDielRecom[ipISO];
				
				TransitionList::iterator tr = SatelliteLines[ipISO][nelem].begin()+ipSatelliteLines[ipISO][nelem][i];
				(*tr).Emis().xObsIntensity() = 
				(*tr).Emis().xIntensity() = 
					dr_rate * dense.eden * dense.xIonDense[nelem][nelem+1-ipISO] *
					ERG1CM * (*tr).EnergyWN();
				
				/* We set line intensity above using a rate, but here we need a transition probability.  
				 * We can obtain this by dividing dr_rate by the population of the autoionizing level.
				 * We assume this level is in statistical equilibrium. */
				factor1 = HION_LTE_POP*dense.AtomicWeight[nelem]/
					(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);
				
				/* term in () is stat weight of electron * ion */
				ConvLTEPOP = powpq(factor1,3,2)/(2.*iso_ctrl.stat_ion[ipISO])/phycon.te32;
				
				if( ConBoltz >= SMALLDOUBLE )
				{
					/* The energy used to calculate ConBoltz above
					 * should be negative since this is above the continuum, but 
					 * to be safe we calculate ConBoltz with a positive energy above
					 * and multiply by it here instead of dividing.  */
					LTE_pop = (*(*tr).Hi()).g() * ConBoltz * ConvLTEPOP;
				}
				
				LTE_pop = max( LTE_pop, 1e-30f );
				
				/* Now the transition probability is simply dr_rate/LTE_pop. */
				(*tr).Emis().Aul() = (realnum)(dr_rate/LTE_pop);
				(*tr).Emis().Aul() = 
					max( iso_ctrl.SmallA, (*tr).Emis().Aul() );
				
				(*tr).Emis().gf() = (realnum)GetGF(
					(*tr).Emis().Aul(),
					(*tr).EnergyWN(),
					(*(*tr).Hi()).g());
				
				(*tr).Emis().gf() = 
					max( 1e-20f, (*tr).Emis().gf() );
				
				(*(*tr).Hi()).Pop() = LTE_pop * dense.xIonDense[nelem][nelem+1-ipISO] * dense.eden;
				// In the approximation used here, DepartCoef is unity by definition.	
				(*(*tr).Hi()).DepartCoef() = 1.;
				
				(*tr).Emis().PopOpc() = 
					(*(*tr).Lo()).Pop() - 
					(*(*tr).Hi()).Pop() * 
					(*(*tr).Lo()).g()/(*(*tr).Hi()).g();
				
				(*tr).Emis().opacity() = 
					(realnum)(abscf((*tr).Emis().gf(),
										 (*tr).EnergyWN(),
										 (*(*tr).Lo()).g()));
				
				/* a typical transition probability is of order 1e10 s-1 */
				double lifetime = 1e-10;
				
				(*tr).Emis().dampXvel() = (realnum)( 
					(1.f/lifetime)/PI4/(*tr).EnergyWN());
			}
		}
	}
	
	return;
}

long iso_get_total_num_levels( long ipISO, long nmaxResolved, long numCollapsed )
{
	DEBUG_ENTRY( "iso_get_total_num_levels()" );

	long tot_num_levels;

	/* return the number of levels up to and including nmaxResolved PLUS 
	 * the number (numCollapsed) of collapsed n-levels */		

	if( ipISO == ipH_LIKE )
	{
		tot_num_levels = (long)( nmaxResolved * 0.5 *( nmaxResolved + 1 ) ) + numCollapsed;
	}
	else if( ipISO == ipHE_LIKE )
	{
		tot_num_levels = nmaxResolved*nmaxResolved + nmaxResolved + 1 + numCollapsed;
	}
	else
		TotalInsanity();

	return tot_num_levels;
}

void iso_update_num_levels( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_update_num_levels()" );

	/* This is the minimum resolved nmax. */
	ASSERT( iso_sp[ipISO][nelem].n_HighestResolved_max >= 3 );

	iso_sp[ipISO][nelem].numLevels_max = 
		iso_get_total_num_levels( ipISO, iso_sp[ipISO][nelem].n_HighestResolved_max, iso_sp[ipISO][nelem].nCollapsed_max );

	if( iso_sp[ipISO][nelem].numLevels_max > iso_sp[ipISO][nelem].numLevels_alloc )
	{
		fprintf( ioQQQ, "The number of levels for ipISO %li, nelem %li, has been increased since the initial coreload.\n",
			ipISO, nelem );
		fprintf( ioQQQ, "This cannot be done.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* set local copies to the max values */
	iso_sp[ipISO][nelem].numLevels_local = iso_sp[ipISO][nelem].numLevels_max;
	iso_sp[ipISO][nelem].nCollapsed_local = iso_sp[ipISO][nelem].nCollapsed_max;
	iso_sp[ipISO][nelem].n_HighestResolved_local = iso_sp[ipISO][nelem].n_HighestResolved_max;

	/* find the largest number of levels in any element in all iso sequences
	 * we will allocate one matrix for ionization solver, and just use a piece of that memory
	 * for smaller models. */
	max_num_levels = MAX2( max_num_levels, iso_sp[ipISO][nelem].numLevels_max);

	return;
}

#if	0
STATIC void Prt_AGN_table( void )
{
	/* the designation of the levels, chLevel[n][string] */
	multi_arr<char,2> chLevel(max_num_levels,10);

	/* create spectroscopic designation of labels */
	for( long ipLo=0; ipLo < iso_sp[ipISO][ipISO].numLevels_max-iso_sp[ipISO][ipISO].nCollapsed_max; ++ipLo )
	{
		long nelem = ipISO;
		sprintf( &chLevel.front(ipLo) , "%li %li%c", N_(ipLo), S_(ipLo), chL[MIN2(20,L_(ipLo))] );
	}

	/* option to print cs data for AGN */
	/* create spectroscopic designation of labels */
	{
		/* option to print particulars of some line when called */
		enum {AGN=false};
		if( AGN )
		{
#			define NTEMP 6
			double te[NTEMP]={6000.,8000.,10000.,15000.,20000.,25000. };
			double telog[NTEMP] ,
				cs ,
				ratecoef;
			long nelem = ipHELIUM;
			fprintf( ioQQQ,"trans");
			for( long i=0; i < NTEMP; ++i )
			{
				telog[i] = log10( te[i] );
				fprintf( ioQQQ,"\t%.3e",te[i]);
			}
			for( long i=0; i < NTEMP; ++i )
			{
				fprintf( ioQQQ,"\t%.3e",te[i]);
			}
			fprintf(ioQQQ,"\n");

			for( long ipHi=ipHe2s3S; ipHi< iso_sp[ipHE_LIKE][ipHELIUM].numLevels_max; ++ipHi )
			{
				for( long ipLo=ipHe1s1S; ipLo < ipHi; ++ipLo )
				{

					/* deltaN = 0 transitions may be wrong because 
					 * COLL_CONST below is only correct for electron colliders */
					if( N_(ipHi) == N_(ipLo) )
						continue;

					/* print the designations of the lower and upper levels */
					fprintf( ioQQQ,"%s - %s",
						 &chLevel.front(ipLo) , &chLevel.front(ipHi) );

					/* print the interpolated collision strengths */
					for( long i=0; i < NTEMP; ++i )
					{
						phycon.alogte = telog[i];
						/* print cs */
						cs = HeCSInterp( nelem , ipHi , ipLo, ipELECTRON );
						fprintf(ioQQQ,"\t%.2e", cs );
					}

					/* print the rate coefficients */
					for( long i=0; i < NTEMP; ++i )
					{
						phycon.alogte = telog[i];
						phycon.te = exp10(telog[i] );
						tfidle(false);
						cs = HeCSInterp( nelem , ipHi , ipLo, ipELECTRON );
						/* collisional deexcitation rate */
						ratecoef = cs/sqrt(phycon.te)*COLL_CONST/iso_sp[ipHE_LIKE][nelem].st[ipLo].g() *
							sexp( iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo).EnergyK() / phycon.te );
						fprintf(ioQQQ,"\t%.2e", ratecoef );
					}
					fprintf(ioQQQ,"\n");
				}
			}
			cdEXIT(EXIT_FAILURE);
		}
	}

	return;
}
#endif
