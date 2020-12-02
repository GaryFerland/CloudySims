/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "iso.h"

#include "two_photon.h"
#include "freebound.h"

t_isoCTRL iso_ctrl;

t_iso_sp iso_sp[NISO][LIMELM];

long int max_num_levels;

void t_isoCTRL::zero()
{
	DEBUG_ENTRY( "t_isoCTRL::zero()" );

	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		/* option to disable continuum lowering */
		lgContinuumLoweringEnabled[ipISO] = true;

		/* flag set by compile he-like command, says to regenerate table of recombination coef */
		lgCompileRecomb[ipISO] = false;
		lgNoRecombInterp[ipISO] = false;

		/* how the gbar cs will be treated - set with atom he-like gbar command */
		/** \todo	2	change this to CS_new */
		lgCS_Vriens[ipISO] = true;
		lgCS_Vrinceanu[ipISO] = false;
		lgCS_PS64[ipISO] = true;
		lgCS_PSClassic[ipISO] = false;
		lgCS_VOS12[ipISO] = false; // lgCS_Vrinceanu[ipISO] == true overrides
		lgCS_VOS12QM[ipISO] = false;
		/*Seaton  M. J. 1962, Proc. Phys. Soc. 79, 1105 treatment for l<=3 is default */
		lgCS_Seaton[ipISO] = true;
		lgCS_B72[ipISO] = false;
		lgCS_PSdeg[ipISO] = true;

		fixit("make this the default for ipH_LIKE if not too slow.");
		lgCS_Vrinceanu[ipH_LIKE] = false;
		lgCS_PS64[ipH_LIKE] = true;

		lgCS_therm_ave[ipISO] = false;
		lgCS_None[ipISO] = false;
		/* when set try actually set to 1 or 2, depending on which fit is to be used,
		 * 1 is the broken power law fit */
		/* >>chng 02 dec 21, change to broken power law fit */
		nCS_new[ipISO] = 1;
		/* This flag says whether the density is high enough that helium is sufficiently l-mixed. */
		lgCritDensLMix[ipISO] = true;
		/* flag saying whether to include fine-structure mixing in spontaneous decays	
		 * set with SPECIES HE-LIKE FSM command */
		lgFSM[ipISO] = 0;
		/* This is the flag saying whether to generate errors.  false means don't.	*/
		lgRandErrGen[ipISO] = false;
		/* this is the flag saying whether we should include excess recombination in the
		 * helike sequence.  Should only be off if testing effect of top off approximations. */
		lgTopoff[ipISO] = true;
		/* Dielectronic recombination for helike ions is on by default.	*/
		lgDielRecom[ipISO] = true;

		/* number of Lyman lines to include in opacities, this can be vastly larger
		 * than the number of actual levels in the model atom */
		nLyman[ipISO] = 100;
		nLyman_max[ipISO] = 100;
		nLyman_malloc[ipISO] = 100;

		/* controls whether l-mixing and collisional ionization included */
		lgColl_l_mixing[ipISO] = true;
		lgColl_excite[ipISO] = true;
		lgColl_ionize[ipISO] = true;
		lgLTE_levels[ipISO] = false;
		lgPrintNumberOfLevels = false;

		for (long nelem=0; nelem<LIMELM; ++nelem)
		{
			RRC_TeUsed[ipISO][nelem]=0.;
		}
	}

	/* Dielectronic recombination forming hydrogen-like ions does not exist. */
	lgDielRecom[ipH_LIKE] = false;

	/* smallest transition probability allowed */
	SmallA = 1e-30f;

	/* reset with SET IND2 command, turns on/off induced two photon */
	lgInd2nu_On = false;

	/* hydrogen redistribution functions */
	ipLyaRedist[ipH_LIKE] = ipPRD;
	ipResoRedist[ipH_LIKE] = ipCRD;
	ipSubRedist[ipH_LIKE] = ipCRDW;

	/* this is the upper level for each Lya, which uses the special ipLY_A */
	nLyaLevel[ipH_LIKE] = ipH2p;
	nLyaLevel[ipHE_LIKE] = ipHe2p1P;

	/* he-like redistribution functions */
	ipLyaRedist[ipHE_LIKE] = ipPRD;
	ipResoRedist[ipHE_LIKE] = ipCRD;
	ipSubRedist[ipHE_LIKE] = ipCRDW;

	lgPessimisticErrors = false;

	/* do not average collision strengths - evaluate at kT 
	 * set true with command SET COLLISION STRENGHTS AVERAGE */
	lgCollStrenThermAver = false;	
}

void t_iso_sp::Reset()
{
	// this is flag indicating which type of model atom to use 
	strcpy( chTypeAtomUsed , "none" );
	CaseBCheck = 0.;
	/* a first guess at the recombination coefficients */
	RadRec_caseB = 1e-13;
	lgLevelsLowered = false;
	lgLevelsEverLowered = false;
	lgMustReeval = false;
	lgPopsRescaled = false;
	/* error generation done yet? false means not done.	*/
	lgErrGenDone = false;
	for( vector<two_photon>::iterator it = TwoNu.begin(); it != TwoNu.end(); ++it )
		(*it).Reset();
	for( vector<freeBound>::iterator it = fb.begin(); it != fb.end(); ++it )
		(*it).Reset();
}
