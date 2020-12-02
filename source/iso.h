/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ISO_H_
#define ISO_H_

/**\file iso.h - information for isoelectronic sequences */
#include "module.h"
#include "transition.h"
#include "container_classes.h"

class two_photon;
class freeBound;

extern long int max_num_levels;

/** This macro is used to zero any radiative process with photon energy below
 * the plasma frequency.  The energy must be in Rydbergs!	*/
#define KILL_BELOW_PLASMA(E_)		( (rfield.lgPlasNu && ((E_)<rfield.plsfrq) ) ? 0.:1. )

/** these macros are just an easy way to return the quantum numbers of a given level. */
#define N_(A_)	(iso_sp[ipISO][nelem].st[A_].n())
#define L_(A_)	(iso_sp[ipISO][nelem].st[A_].l())
#define S_(A_)	(iso_sp[ipISO][nelem].st[A_].S())
#define J_(A_)	(iso_sp[ipISO][nelem].st[A_].j())

/** some levels for hydrogenic species */

const int ipH1s = 0;
const int ipH2s = 1;
const int ipH2p = 2;
const int ipH3s = 3;
const int ipH3p = 4;
const int ipH3d = 5;
const int ipH4s = 6;
const int ipH4p = 7;
const int ipH4d = 8;
const int ipH4f = 9;

/** some levels for he-like species */

/* level 1 */
const int ipHe1s1S = 0;

/* level 2 */
const int ipHe2s3S = 1;
const int ipHe2s1S = 2;
const int ipHe2p3P0 = 3;
const int ipHe2p3P1 = 4;
const int ipHe2p3P2 = 5;
const int ipHe2p1P = 6;

/* level 3 */
const int ipHe3s3S = 7;
const int ipHe3s1S = 8;
const int ipHe3p3P = 9;
const int ipHe3d3D = 10;
const int ipHe3d1D = 11;
const int ipHe3p1P = 12;

/** these are array indices for isoelectronic sequences,
 * same as element but used for array addressing to make
 * context totally clear */
const int ipH_LIKE = 0;
const int ipHE_LIKE = 1;
const int ipLI_LIKE = 2;
const int ipBE_LIKE = 3;
const int ipB_LIKE = 4;
const int ipC_LIKE = 5;
const int ipN_LIKE = 6;
const int ipO_LIKE = 7;
const int ipF_LIKE = 8;
const int ipNE_LIKE = 9;
const int ipNA_LIKE = 10;
const int ipMG_LIKE = 11;
const int ipAL_LIKE = 12;
const int ipSI_LIKE = 13;
const int ipP_LIKE = 14;
const int ipS_LIKE = 15;
const int ipCL_LIKE = 16;
const int ipAR_LIKE = 17;

enum {
	ipSINGLET = 1, ipDOUBLET = 2, ipTRIPLET = 3, 
	ipMULTIPLET_END, ipMULTIPLET_BEGIN=ipSINGLET
};

#define IPRAD	0
#define IPCOLLIS	1
/*define IPENERGY	2*/

/* following two macros used to define recombination coef arrays */
/* Max n desired in RRCoef file.	*/
/** this is the number of levels used with the
 * atom xx-like levels large command */
/* Hydrogen and helium atoms will have precompiled recombination coefficients up to these maximum n.	*/
#define RREC_MAXN	40

/** Ions of the sequences will go up to this n, h-like He will get same as iso roots.	*/
#define LIKE_RREC_MAXN( A_ )		( A_ == ipHELIUM ? 40 : 20 )

#define N_ISO_TE_RECOMB		41

/** This is the n to go up to when calculating total recombination.	Any change 
 * here will not be reflected in total recomb until "compile xxlike" is run	*/
#define SumUpToThisN	1000
/** the magic number for the table of recombination coefficients, YYMMDD */
#define RECOMBMAGIC		(130216)

/**iso_cascade - calculate cascade probabilities, branching ratios, and associated errors
\param ipISO
\param nelem
*/
void iso_cascade( long ipISO, long nelem );

/**iso_charge_transfer_update - update rate coefficients for CT of H and He with everything else
 */
void iso_charge_transfer_update( long nelem );

/**iso_collapsed_Aul_update - update decays from collapsed levels 
\param ipISO
\param nelem
*/
void iso_collapsed_Aul_update( long ipISO, long nelem );

/**iso_collapsed_lifetimes_update - update lifetimes of collapsed levels 
\param ipISO
\param nelem
*/
void iso_collapsed_lifetimes_update( long ipISO, long nelem );

/** iso_collide - calculate collision data for ipISO, nelem  
\param ipISO
\param nelem
*/
void iso_collide( long ipISO, long nelem );

/** iso_collisional_ionization - calculate collisional ionization rate for ipISO, nelem  
\param ipISO
\param nelem
*/
void iso_collisional_ionization( long ipISO, long nelem );

/** iso_continuum_lower - limit max prin. quan. no. due to continuum lowering processes 
\param ipISO
\param nelem
*/
void iso_continuum_lower( long ipISO , long nelem );

/**iso_cool compute net heating/cooling due to hydrogenc atom species 
\param ipISO the isoelectronic sequence, 0 for H 
\param nelem is element, so 0 for H itself
*/
void iso_cool( long ipISO , long nelem );

/** iso_setRedisFun assign the line redistribution function type
 * \param ipISO isoelectronic sequence
 * \param nelem element index
 * \param ipLo  index to lower state
 * \param ipHi  index to upper state
 */
void iso_setRedisFun (long ipISO, long nelem, long ipLo, long ipHi);

/** iso_setOpacity compute line opacity
 * \param ipISO isoelectronic sequence
 * \param nelem element index
 * \param ipLo  index to lower state
 * \param ipHi  index to upper state
 */
void iso_setOpacity (long ipISO, long nelem, long ipLo, long ipHi);

/**iso_create create storage space data for iso sequences, 1 one time per coreload 
*/
void iso_create( void );

/**iso_cross_section get cross section for a particular level of an iso sequence ion
\param ERyd 
\param EthRyd
\param n
\param l
\param S
\param Z 
\param ipISO
*/
double iso_cross_section( double ERyd , double EthRyd, long n, long l, long S, long globalZ, long globalISO );

/**iso_departure_coefficients - calculate departure coefficients
\param ipISO
\param nelem
*/
void iso_departure_coefficients( long ipISO, long nelem );

/**iso_dielec_recomb_rate - get state-specific dielectronic recombination rate 
\param ipISO
\param nelem
\param ipLo
*/
double iso_dielec_recomb_rate( long ipISO, long nelem, long ipLo );

/**iso_error_generation generate gaussian errors 
\param ipISO
\param nelem
*/
void iso_error_generation( long ipISO, long nelem );

/**iso_get_total_num_levels - get total number of levels with the given number of resolved and collapsed
\param ipISO
\param nmaxResolved
\param numCollapsed
*/
long iso_get_total_num_levels( long ipISO, long nmaxResolved, long numCollapsed );

/**IonHydro this controls hydrogen atomic and molecular crosstalk
*/
void IonHydro( );

/**iso_ionize_recombine evaluate state specific creation and destruction processes 
\param ipISO
\param nelem
*/
void iso_ionize_recombine( long ipISO , long nelem );

/**iso_level solve for iso-sequence ionization balance 
 * \param ipISO		iso-sequence index
 * \param nelem		element index
 * \param renorm	renormalization parameter
 * \param lgPrtMatrix	boolean to print the matrix to be solved
*/
void iso_level( const long ipISO, const long nelem, double& renorm,
		bool lgPrtMatrix );

/**iso_photo do photoionization rates for element nelem on the ipISO isoelectronic sequence 
\param ipISO
\param nelem
*/
void iso_photo( long ipISO , long nelem );

/**iso_prt_pops routine to print level pops or departure coefficients for iso sequences 
\param ipISO
\param nelem
\param lgPrtDeparCoef
*/
void iso_prt_pops( long ipISO, long nelem, bool lgPrtDeparCoef );

/**iso_put_error put an error bar on a piece of data, to be used with Gaussian random noise gen 
\param ipISO
\param nelem
\param ipHi
\param ipLo
\param whichData
\param errorOpt
\param errorPess
*/
void iso_put_error( long ipISO, long nelem, long ipHi, long ipLo, long whichData, realnum errorOpt, realnum errorPess);

/**iso_radiative_recomb - get rad recomb rate coefficients for iso sequences.
\param ipISO
\param nelem
*/
void iso_radiative_recomb( long ipISO, long nelem );

/**iso_radiative_recomb_effective - get effective recomb rate coefficients into each level (including indirect)
\param ipISO
\param nelem
*/
void iso_radiative_recomb_effective( long ipISO, long nelem );

/**iso_recomb_check - called by SanityCheck to confirm that recombination coef are ok,
 * return value is relative error between new calculation of recom, and interp value 
 \param ipISO
 \param nelem the chemical element, 1 for He
 \param level the level, 0 for ground
 \param temperature the temperature to be used
*/
double iso_recomb_check( long ipISO, long nelem, long level, double temperature );

/** iso_recomb_auxiliary_free - free up some auxiliary space associated with iso recombination tables.
*/
void iso_recomb_auxiliary_free( void );

/** iso_recomb_malloc - malloc space needed for iso recombination tables.
*/
void iso_recomb_malloc( void );

/** iso_recomb_setup - read in or compile iso recombination tables.
\param ipISO
*/
void iso_recomb_setup( long ipISO );

/** iso_RRCoef_Te - interpolate iso recomb coeff as function of temperature
\param ipISO
\param nelem
\param temp - the electron temperature
\param n
*/
double iso_RRCoef_Te( long ipISO, long nelem, double temp, long n );

/**iso_satellite_update - update iso satellite line information 
*/
void iso_satellite_update( long nelem );

/* calculate radiative lifetime of an individual iso state 
\param ipISO
\param nelem
\param n
\param l
*/
double iso_state_lifetime( long ipISO, long nelem, long n, long l );

/**iso_solve - main routine to call iso_level and determine iso level balances
\param ipISO
*/
void iso_solve( long ipISO, long nelem, double &maxerr );

/**iso_suprathermal - calculate secondary excitation by suprathermal electrons for iso sequences 
\param ipISO
\param nelem
*/
void iso_suprathermal( long ipISO, long nelem );

/**iso_update_num_levels - update level informations for iso sequences 
\param ipISO
\param nelem
*/
void iso_update_num_levels( long ipISO, long nelem );

/**iso_update_rates routine to set up iso rates, level balance is done elsewhere 
*/
void iso_update_rates( void );

void iso_collapsed_update( void );

void iso_set_ion_rates( long ipISO, long nelem);

class t_isoCTRL : public module
{
public:
	void zero();
	void comment(t_warnings&) {}

	const char *chName() const
	{
		return "iso_ctrl";
	}

	bool lgPrintNumberOfLevels;

	const char *chISO[NISO];

	/** number of Lyman lines to include only as opacity sources, in each iso seq,
	 * all now set to 100 in zero.c */
	long int nLyman[NISO],
		/** max number of levels to consider - probably greater than above */
		nLyman_max[NISO],
		/** number of levels actually malloc'd - probably greater than above */
		nLyman_malloc[NISO];

	/** option to turn off l-mixing collisions */
	bool lgColl_l_mixing[NISO];

	/** option to turn off collisional excitation */
	bool lgColl_excite[NISO];

	/** option to turn off collisional ionization */
	bool lgColl_ionize[NISO];

	bool lgLTE_levels[NISO];

	/** do thermal average of collision strengths if true, false by default,
	 * set true with SET COLLISION STRENGTHS AVERAGE command */
	bool lgCollStrenThermAver;

	/** flag saying whether induced two photon is included
	 * in the level pops for H- and He-like */
	bool lgInd2nu_On;

	/* option to disable continuum lowering due to stark broadening, particle packing, etc. */
	bool lgContinuumLoweringEnabled[NISO];

	/** statistical weight of the ground state of the parent ions for each
	 * species, used for Milne relation and recombination */
	realnum stat_ion[NISO];

	/** tells whether dielectronic recombination is turned on	*/
	bool lgDielRecom[NISO];

	/** this is the rate for the Aul given to bogus transitions,
	 * set to 1e-30 in zero */
	/** >>chng 04 may 17, esd 1e-20, changed to 1e-30 to allow
	 * rydberg levels to be treated with their small As */
	realnum SmallA;

	/** types of redistribution functions for Lya, other resonances, and subordinate lines */
	int ipLyaRedist[NISO] , ipResoRedist[NISO] , ipSubRedist[NISO];

	/** this is the upper level for Lya */
	int nLyaLevel[NISO];

	/** flag set by compile he-like command, says to regenerate table of recombination coef */
	bool lgCompileRecomb[NISO];

	/** flag set by atom he-like no recomb interp command,
	 * says to generate recombination coefficients
	 * on the fly */
	bool lgNoRecombInterp[NISO];

	/** parameters for changing gbar - set with set hegbar command */
	bool lgCS_Vriens[NISO] ,
		lgCS_None[NISO] ,
		lgCS_Seaton[NISO],
		lgCS_B72[NISO],
		lgCS_PSdeg[NISO],
		lgCS_Vrinceanu[NISO],
		lgCS_PS64[NISO],
		lgCS_PSClassic[NISO],
		lgCS_VOS12[NISO],
		lgCS_VOS12QM[NISO],
		lgCS_therm_ave[NISO];
	int nCS_new[NISO];//vals are 0, 1, and 2

	/** used to print warning if density too low for first collapsed level to be l-mixed	*/
	bool lgCritDensLMix[NISO];

 	/** flag saying whether to include fine-structure mixing in spontaneous decays
	 * set with SPECIES HE-LIKE FSM command */
	bool lgFSM[NISO];

	/** This flag is set to true if the rates should be treated with a randomly generated error,
	 * on the range specifically set for each rate, before being entered into the rate matrix.	*/
	bool lgRandErrGen[NISO];


	bool lgPessimisticErrors;

	bool lgTopoff[NISO];

	/** This is the used to set a unique seed in parallel gaussian runs */
	int modelRank[NISO];

	double RRC_TeUsed[NISO][LIMELM];
};

extern t_isoCTRL iso_ctrl;

class extra_tr
{
public:
	/** stark broadening in Puetter formalism */
	double pestrk;
	double pestrk_up;

	/* NB NB NB ---  Error and ErrorFactor need one more slot than all the rest of these! */
	/* and the last dimension can just be hardwired to 3 */

	/** This is the array in which uncertainties are stored if lgRandErrGen is set. */
	/* first dimension is upper level,
	 * second is lower level,
	 * third is for radiative, collisional, or energy errors.
	 * MACROS are used for the last dimension: IPRAD, IPCOLLIS, and IPENERGY. */
	realnum Error[3];

	/** This is the array in which gaussian errors are generated, using the values in
	 * the Error array above as the standard deviations */
	realnum ErrorFactor[3];

	/** total brancing ratio and standard deviation in it */
	double SigmaCascadeProb;
};

class t_iso_sp
{
public:
	TransitionProxy trans( const long ipHi, const long ipLo ) 
	{
		return (*tr)[ ipTrans[ipHi][ipLo] ];
	}
	multi_arr<long,2> ipTrans;
	multi_arr<extra_tr,2> ex;
	multi_arr<double,2> CascadeProb;
	multi_arr<double,2> BranchRatio;
	vector<freeBound> fb;
	qList st;
	TransitionList* tr;

	/** Find index given quantum numbers
	 * Since separate j levels within a triplet term are only resolved in the case of 2tripP,
	 * allocating memory for a j dimension is unwarranted.  Instead 
	 * iso.QuantumNumbers2Index[ipISO][nelem][2][1][1] will point to 2^3P2, with 2^3P0 and 2^3P1
	 * easily accessed by subtracting 2 and 1 respectively from the returned index.	*/
	multi_arr<long,3> QuantumNumbers2Index;
	// This is the same as QuantumNumbers2Index, except indexed as if collapsed levels were resolved.
	multi_arr<long,3> IndexIfAllResolved;

	/** the ratio of ion to atom for all iso species
	 * xIonSimple is simple estimate, should agree at low density */
	double xIonSimple;

	/** option to print departure coefficients */
	bool lgPrtDepartCoef;

	/** option to print level populations */
	bool lgPrtLevelPops;

	/** true if the number of levels is currently lowered */
	bool lgLevelsLowered;

	/** This variable is set to true if the continuum was lowered at any point in the calculation.
	 * Necessary because some models will lowered continuum at intermediate points but not last zone. */
	bool lgLevelsEverLowered;

	/* flag that says we must reevaluate everything about this ion */
	bool lgMustReeval;

	/* set true if "element ionization" forces rescaling of pops */
	bool lgPopsRescaled;

	/** the number of collapsed levels, these lie on top of resolved levels */
	long int nCollapsed_max;
	long int nCollapsed_local;

	/** total number of collapsed and resolve levels, 
	 * numLevels_max is derived from total resolved and collapsed levels 
	 * it is the maximum number of levels ever to be used in this core load. */
	long int numLevels_max;

	/** total number of levels with continuum pressure lowering included 
	 * this varies from zone to zone, and from model to model, but cannot
	 * exceed numLevels_max  */
	long int numLevels_local;

	/** number of levels malloc'd in the core load, can't go over that later 
	 * in later sims can lower number of levels but not raise them  */
	long int numLevels_malloc;

	/** principal quantum number n of the highest resolved level */
	long int n_HighestResolved_max;
	/** the local (pressure lowered) version of the above */
	long int n_HighestResolved_local;

	/** difference between actual case b photons in rtdiffuse, and correct case b */
	realnum CaseBCheck;

	/** case b recombination rate coefficient */
	double RadRec_caseB;

	/** the total effective radiative recombination rate coefficient (cm3 s-1), 
	 * radiative rate with correction for absorption and ionization */
	double RadRec_effec;

	/** ratio of collisional recombination rate to recom from all processes */
	double RecomCollisFrac;

	/** true is all lte populations positive for Hydrogenic atoms */
	bool lgPopLTE_OK;

	/** net free bound cooling for this element */
	double FreeBnd_net_Cool_Rate;

	/** net cooling due to collisional ionization */
	double coll_ion;

	/** net cooling due to collisional excit of higher lines */
	double cRest_cool;

	/** net cooling due to total collisional excit of lines */
	double xLineTotCool;

	/** deriv of net cooling due to total collisional excit of lines */
	double dLTot;

	/** net cooling due to rad rec */
	double RadRecCool;

	/** net cooling due to collisional excit of balmer lines */
	double cBal_cool;

	/** net cooling due to collisional excit of higher lyman lines */
	double cLyrest_cool;

	/** net cooling due to collisional excit of Lya */
	double cLya_cool;

	/** the actual induced recom cooling rate, erg cm-3 s-1 */
	double RecomInducCool_Rate;

	/** flag to set which type of solution was used for level pops, "zero" or "popul" */
	char chTypeAtomUsed[10];

	/** this is flag saying that random gaussians have already been set...they should only
	 * be done once per model, and this must be reset to false at the beginning of each model.	*/
	bool lgErrGenDone;

	/** the effective collisional rate from 2S, for h-like and he-like sequences */
	double qTot2S;

	multi_arr<realnum,3> CachedAs;
	void Reset();
	vector<two_photon> TwoNu;

	vector<double> HighestLevelOpacStack;

	/** print Matrix input to solver */
	bool lgPrtMatrix;
};

extern t_iso_sp iso_sp[NISO][LIMELM];

/** iso_renorm - renormalize H-like so that it agrees with the ionization balance */
void iso_renorm( long nelem, long ipISO, double& renorm );

/** update multiplet opacities */
void iso_multiplet_opacities( void );

/* iso_comment_tran_levels - prepare comment string for entry to the line stack.
 * The comment has the form 'H-like, 1 3, 1^2S - 2^2P', where '1 3' are the energy
 * level indices, 1 being ground.
 *
 * \param ipISO		iso-sequence index
 * \param nelem		element index
 * \param ipLo, ipHi	lower and upper level indices
 * \return		comment string
 */
string iso_comment_tran_levels( long ipISO, long nelem, long ipLo, long ipHi );

#endif /* ISO_H_ */
