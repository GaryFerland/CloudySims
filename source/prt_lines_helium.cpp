/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_helium put He-like iso sequence into line intensity stack */
/*TempInterp interpolates on a grid of values to produce predicted value at current Te.*/
#include "cddefines.h"
#include "dense.h"
#include "prt.h"
#include "helike.h"
#include "iso.h"
#include "atmdat.h"
#include "lines.h"
#include "phycon.h"
#include "taulines.h"
#include "thirdparty.h"
#include "trace.h"
#include "freebound.h"
#include "two_photon.h"
#include "lines_service.h"

#define NUMTEMPS	21
#define NUMDENS		14		
typedef struct 
{
	/* index for upper and lower levels of line */
	long int ipHi;
	long int ipLo;

	char label[5];
} stdLines;

STATIC void GetStandardHeLines(void);
STATIC double TempInterp2( double* TempArray , double* ValueArray, long NumElements, double Te );
STATIC void DoSatelliteLines( long nelem );

static bool lgFirstRun = true;
static double CaBDensities[NUMDENS];
static double CaBTemps[NUMTEMPS];
static long NumLines;
static double ****CaBIntensity;
static stdLines **CaBLines;



STATIC void insert_trans( vector<TransitionProxy> &trList, TransitionProxy tr )
{
	if( tr.ipCont() < 1 )
		return;
	trList.push_back( tr );
	return;
}



STATIC void multiplet_sum( vector<TransitionProxy> &trList, realnum &av_wl,
			double &sum_inten, double &sum_obs_inten, double &sum_cool, double &sum_heat )
{
	av_wl = -1.;
	sum_inten = sum_obs_inten = sum_cool = sum_heat = 0.;

	if( trList.size() == 0)
		return;

	vector<realnum> wl;
	for( vector<TransitionProxy>::iterator tr = trList.begin(); tr != trList.end(); tr++ )
	{
		sum_obs_inten += (*tr).Emis().xObsIntensity();
		sum_inten += (*tr).Emis().xIntensity();
		sum_cool += (*tr).Coll().cool();
		sum_heat += (*tr).Coll().heat();
		wl.push_back( (*tr).WLAng() );
	}

	if( wl.size() == 0 )
		return;

	std::sort( wl.begin(), wl.end(), std::less<realnum>() );
	if( wl.size() >= 1)
		av_wl = wl[1];
	else	av_wl = wl[0];

	return;
}


STATIC inline void PutLineSum ( const TransitionProxy &tr, const realnum av_wl,
				const double sumxInt, const double sumxObsInt,
				const double sumcool, const double sumheat,
				const char *chComment )
{
	if( av_wl < 0. )
		return;

	TransitionProxy tr2 = tr;

	tr2.WLAng() = av_wl;
	tr2.Emis().xIntensity() = sumxInt;
	tr2.Emis().xObsIntensity() = sumxObsInt;
	tr2.Coll().cool() = sumcool;
	tr2.Coll().heat() = sumheat;

	PutLine( tr2, chComment );

	return;
}

STATIC inline realnum get_multiplet_wl( vector<TransitionProxy> &multiHe, long ipHi, long ipLo )
{
	realnum wl = -1.;

	for( vector<TransitionProxy>::iterator tr = multiHe.begin(); tr != multiHe.end(); tr++ )
	{
		if( (*tr).ipHi() == ipHi && (*tr).ipLo() == ipLo )
		{
			wl = (*tr).WLAng();
			break;
		}
	}

	return	wl;
}

STATIC inline void randomize_inten( t_iso_sp* sp, long ipLo, long ipHi ) 
{ 
	sp->trans(ipHi,ipLo).Emis().xIntensity() *= sp->ex[ipHi][ipLo].ErrorFactor[IPRAD]; 
	sp->trans(ipHi,ipLo).Emis().xObsIntensity() *= sp->ex[ipHi][ipLo].ErrorFactor[IPRAD]; 
	return; 
} 


void lines_helium(void)
{
	long ipISO = ipHE_LIKE;
	string chLabel="    ";

	double log10_eden = log10(dense.eden);

	DEBUG_ENTRY( "lines_helium()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "   prt_lines_helium called\n" );

	// this can be changed with the atom levels command but must be at
	// least 3.
	ASSERT( iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max >= 3 );

	long i = StuffComment( "He-like iso-sequence" );
	linadd( 0., (realnum)i , "####", 'i',
		" start He-like iso sequence");

	/* read in Case A and B lines from data file	*/
	if( lgFirstRun )
	{
		GetStandardHeLines();
		lgFirstRun = false;
	}

	/* this is the main printout, where line intensities are entered into the stack */
	for( long nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		vector<TransitionProxy> multiplet;
		vector<TransitionProxy> HeTrList;
		double sumxInt = 0., sumxObsInt = 0., sumcool = 0., sumheat = 0.;
		realnum	av_wl = 0.;


		if( dense.lgElmtOn[nelem] )
		{
			t_iso_sp* sp = &iso_sp[ipHE_LIKE][nelem];

			ASSERT( sp->n_HighestResolved_max >= 3 );

			// add two-photon details here
			if( LineSave.ipass == 0 )
			{
				/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" 
				 * the result, chLable, is only used when ipass == 0, can be undefined otherwise */
				chLabel = chIonLbl(nelem+1, nelem+1-ipISO);
			}
			for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
			{
				fixit("This was multiplied by Pesc when treated as a line, now what?  Only used for printout?");
				fixit("below should be 'i' instead of 'r' ?");

				string tpc_comment = "";
				if( LineSave.ipass == 0 )
				{
					tpc_comment = " two photon continuum, " +
						iso_comment_tran_levels( ipISO, nelem, (*tnu).ipLo, (*tnu).ipHi );
				}

				linadd(	tnu->AulTotal * tnu->E2nu * EN1RYD * (*tnu->Pop), 
					2. * wn2ang( (*sp).trans( (*tnu).ipHi, (*tnu).ipLo ).EnergyWN() ),
					chLabel.c_str(), 'r', tpc_comment.c_str() );
			}

			/* here we will create an entry for the three lines 
			 * coming from 2 3P to 1 1S - the entry called TOTL will
			 * appear before the lines of the multiplet */

			for( long i=ipHe2p3P0; i <= ipHe2p3P2; i++ )
			{
				set_xIntensity( sp->trans(i, ipHe1s1S) );
				insert_trans( multiplet, sp->trans(i, ipHe1s1S) );
				if( nelem == ipHELIUM )
					HeTrList.push_back( sp->trans(i, ipHe1s1S) );
			}

			multiplet_sum( multiplet, av_wl, sumxInt, sumxObsInt, sumcool, sumheat );
			multiplet.resize( 0 );

			linadd(sumxObsInt, sp->trans(ipHe2p3P1,ipHe1s1S).WLAng(), "TOTL", 'i',
				" total emission in He-like intercombination lines from 2p3P to ground ");


			/* set number of levels we want to print, first is default,
			 * only print real levels, second is set with "print line
			 * iso collapsed" command */
			long int nLoop  = sp->numLevels_max - sp->nCollapsed_max;
			if( prt.lgPrnIsoCollapsed )
				nLoop  = sp->numLevels_max;

			/* now do real permitted lines */
			/* NB NB - low and high must be in this order so that all balmer, paschen,
			 * etc series line up correctly in final printout */
			/* >>chng 01 jun 13, bring 23P lines back together */
			for( long ipLo=0; ipLo < ipHe2p3P0; ipLo++ )
			{
				vector<long> EnterTheseLast;
				for( long ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* >>chng 01 may 30, do not add fake he-like lines (majority) to line stack */
					/* >>chng 01 dec 11, use variable for smallest A */
					if( sp->trans(ipHi,ipLo).ipCont() < 1 )
						continue;

					/* recombine fine-structure lines since the energies are
					 * not resolved anyway.	*/
					if( iso_ctrl.lgFSM[ipISO] && ( abs(sp->st[ipHi].l() -
						sp->st[ipLo].l())==1 )
						&& (sp->st[ipLo].l()>1) 
						&& (sp->st[ipHi].l()>1) 
						&& ( sp->st[ipHi].n() ==
						sp->st[ipLo].n() ) )
					{
						/* skip if both singlets. */
						if( (sp->st[ipHi].S()==1) 
							&& (sp->st[ipLo].S()==1) )
						{
							continue;
						}
						else if( (sp->st[ipHi].S()==3) 
							&& (sp->st[ipLo].S()==3) )
						{
							string comment_trans = "";

							/* singlet to singlet*/
							insert_trans( multiplet, sp->trans(ipHi  ,ipLo+1) );
							insert_trans( multiplet, sp->trans(ipHi+1,ipLo+1) );
							if( nelem == ipHELIUM )
							{
								HeTrList.push_back( sp->trans(ipHi  , ipLo+1) );
								HeTrList.push_back( sp->trans(ipHi+1, ipLo+1) );
							}
							if( multiplet.size() == 0 ) continue;

							multiplet_sum( multiplet, av_wl, sumxInt, sumxObsInt, sumcool, sumheat );
							multiplet.resize( 0 );

							if( LineSave.ipass == 0 )
							{
								comment_trans = "singlet to singlet: ";
								comment_trans += iso_comment_tran_levels( ipISO, nelem, ipLo+1, ipHi );
								comment_trans += "; ";
								comment_trans += iso_comment_tran_levels( ipISO, nelem, ipLo+1, ipHi+1 );
							}
							PutLineSum( sp->trans(ipHi+1,ipLo+1),
								sp->trans(ipHi+1,ipLo+1).WLAng(), sumxInt, sumxObsInt, sumcool, sumheat,
								comment_trans.c_str() );

							/* triplet to triplet */
							insert_trans( multiplet, sp->trans(ipHi  ,ipLo) );
							insert_trans( multiplet, sp->trans(ipHi+1,ipLo) );
							if( nelem == ipHELIUM )
							{
								HeTrList.push_back( sp->trans(ipHi  , ipLo) );
								HeTrList.push_back( sp->trans(ipHi+1, ipLo) );
							}
							if( multiplet.size() == 0 ) continue;

							multiplet_sum( multiplet, av_wl, sumxInt, sumxObsInt, sumcool, sumheat );
							multiplet.resize( 0 );

							if( LineSave.ipass == 0 )
							{
								comment_trans = "triplet to triplet: ";
								comment_trans += iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
								comment_trans += "; ";
								comment_trans += iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi+1 );
							}
							PutLineSum( sp->trans(ipHi,ipLo),
								sp->trans(ipHi,ipLo).WLAng(), sumxInt, sumxObsInt, sumcool, sumheat,
								comment_trans.c_str() );
						}
					}

					else if( ipLo==ipHe2s3S && ipHi == ipHe2p3P0 )
					{
						/* here we will create an entry for the three lines 
						 * coming from 2 3P to 2 3S - the entry called TOTL will
						 * appear before the lines of the multiplet 
						 * for He I this is 10830 */

						for( long i=ipHe2p3P0; i <= ipHe2p3P2; i++ )
						{
							set_xIntensity( sp->trans(i, ipLo) );

							/* >>chng 13-jun-06
							 * correct for isotropic continuum before applying error randomization
							 * to avoid shutting off emission lines (if the correction is applied _after_
							 * the error is computed, it is likely to be higher than the updated intensity)
							 */
							if( iso_ctrl.lgRandErrGen[ipISO] )
							{
								randomize_inten( sp, ipLo, i );
							}

							insert_trans( multiplet, sp->trans(i, ipLo) );
							if( nelem == ipHELIUM )
								HeTrList.push_back( sp->trans(i, ipLo) );
						}

						if( multiplet.size() == 0 ) continue;
						multiplet_sum( multiplet, av_wl, sumxInt, sumxObsInt, sumcool, sumheat );
						multiplet.resize( 0 );


						linadd(sumxObsInt, av_wl, "TOTL", 'i',
							"total emission in He-like lines, use average of three line wavelengths " );


						if( nelem == ipHELIUM )
						{
							string comment_trans = "";
							if( LineSave.ipass == 0 )
							{
								comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
							}
							PutLineSum( sp->trans(ipHi,ipLo), av_wl, sumxInt, sumxObsInt, sumcool, sumheat,
									comment_trans.c_str() );
						}
						else
						{
							string comment_trans = "";
							if( LineSave.ipass == 0 )
							{
								comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
							}
							PutLine( sp->trans(ipHi,ipLo), comment_trans.c_str() );

							/* also add this with the regular label, so it is correctly picked up by assert case-b command */
							linadd(sumxObsInt, av_wl, chLabel.c_str(), 'i',
								"total emission in He-like lines, use average of three line wavelengths " );
						}
					}
					else if( ipLo==ipHe2s3S && (ipHi == ipHe2p3P1 || ipHi==ipHe2p3P2) )
					{
						/* >>chng 13-aug-01
						 * If He, do nothing.
						 * The transitions have been computed already in the loop above.
						 * N.B. If this block is removed, the transitions will be rentered
						 *	by the block below.
						 */
						if( nelem > ipHELIUM )
						{
							string comment_trans = "";
							if( LineSave.ipass == 0 )
							{
								comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
							}
							PutLine( sp->trans(ipHi,ipLo), comment_trans.c_str() );
						}
					}
					else
					{
						set_xIntensity( sp->trans(ipHi,ipLo) );

						if( iso_ctrl.lgRandErrGen[ipISO] )
						{
							randomize_inten( sp, ipLo, ipHi );
						}

						if( abs( L_(ipHi) - L_(ipLo) ) != 1 )
						{
							EnterTheseLast.push_back( ipHi );
							continue;
						}

						/* 
						fprintf(ioQQQ,"1 loop %li %li %.1f\n", ipLo,ipHi, 
							sp->trans(ipHi,ipLo).WLAng() ); */
						string comment_trans = "";
						if( LineSave.ipass == 0 )
						{
							comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
						}
						PutLine( sp->trans(ipHi,ipLo), comment_trans.c_str() );

						{
							/* option to print particulars of some line when called
							 * a prettier print statement is near where chSpin is defined below*/
							enum {DEBUG_LOC=false};
							if( DEBUG_LOC )
							{
								if( nelem==1 && ipLo==0 && ipHi==1 )
								{
									fprintf(ioQQQ,"he1 626 %.2e %.2e \n", 
										sp->trans(ipHi,ipLo).Emis().TauIn(),
										sp->trans(ipHi,ipLo).Emis().TauTot()
										);
								}
							}
						}
					}
				}

				// enter these lines last because they are generally weaker quadrupole transitions.
				for( vector<long>::iterator it = EnterTheseLast.begin(); it != EnterTheseLast.end(); it++ )
				{
					string comment_trans = "";
					if( LineSave.ipass == 0 )
					{
						comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, *it );
					}
					PutLine( sp->trans(*it,ipLo), comment_trans.c_str() );
				}
			}

			/* this sum will bring together the three lines going to J levels within 23P */
			for( long ipHi=ipHe2p3P2+1; ipHi < nLoop; ipHi++ )
			{
				for( long ipLo=ipHe2p3P0; ipLo <= ipHe2p3P2; ++ipLo )
				{
					if( sp->trans(ipHi,ipLo).ipCont() < 1 ) 
						continue;

					set_xIntensity( sp->trans(ipHi, ipLo) );

					if( iso_ctrl.lgRandErrGen[ipISO] )
					{
						randomize_inten( sp, ipLo, ipHi );
					}

					insert_trans( multiplet, sp->trans(ipHi, ipLo) );
					if( nelem == ipHELIUM )
						HeTrList.push_back( sp->trans(ipHi, ipLo) );
				}

				if( sp->trans(ipHi,ipHe2p3P2).ipCont() < 1 )
				{
					multiplet.resize( 0 );
					continue;
				}

				if( multiplet.size() == 0 ) continue;
				multiplet_sum( multiplet, av_wl, sumxInt, sumxObsInt, sumcool, sumheat );
				multiplet.resize( 0 );

				string comment_trans = "";
				if( LineSave.ipass == 0 )
				{
					comment_trans = iso_comment_tran_levels( ipISO, nelem, ipHe2p3P2, ipHi );
				}
				PutLineSum( sp->trans(ipHi,ipHe2p3P2), av_wl, sumxInt, sumxObsInt, sumcool, sumheat,
						comment_trans.c_str() );
			}
			for( long ipLo=ipHe2p3P2+1; ipLo < nLoop-1; ipLo++ )
			{
				vector<long> EnterTheseLast;
				for( long ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* skip non-radiative lines */
					if( sp->trans(ipHi,ipLo).ipCont() < 1 ) 
						continue;

					set_xIntensity( sp->trans(ipHi,ipLo) );

					if( iso_ctrl.lgRandErrGen[ipISO] )
					{
						randomize_inten( sp, ipLo, ipHi );
					}

					if( N_(ipHi) > sp->n_HighestResolved_max || abs( L_(ipHi) - L_(ipLo) ) != 1 )
					{
						EnterTheseLast.push_back( ipHi );
						continue;
					}

					string comment_trans = "";
					if( LineSave.ipass == 0 )
					{
						comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
					}
					PutLine(sp->trans(ipHi,ipLo), comment_trans.c_str() );
				}

				// enter these lines last because they are generally weaker quadrupole transitions.
				for( vector<long>::iterator it = EnterTheseLast.begin(); it != EnterTheseLast.end(); it++ )
				{
					string comment_trans = "";
					if( LineSave.ipass == 0 )
					{
						comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, *it );
					}
					PutLine( sp->trans(*it,ipLo), comment_trans.c_str() );
				}
			}

			/* Now put the satellite lines in */
			if( iso_ctrl.lgDielRecom[ipISO] )
				DoSatelliteLines(nelem);

			if( nelem == ipHELIUM )
			{
				for( long i=0; i< NumLines; i++ )
				{
					double intens_at_Te[NUMDENS];
					for( long ipDens = 0; ipDens < NUMDENS; ++ipDens )
						intens_at_Te[ipDens] = TempInterp2( CaBTemps, CaBIntensity[nelem][i][ipDens], NUMTEMPS, phycon.te );
					double intens = linint( CaBDensities, intens_at_Te, NUMDENS, log10_eden );
					intens = exp10(  intens ) * dense.xIonDense[nelem][nelem+1-ipISO]*dense.eden;
					ASSERT( intens >= 0. );

					realnum wvlg = get_multiplet_wl( HeTrList, CaBLines[nelem][i].ipHi, CaBLines[nelem][i].ipLo );
					if( wvlg <= 0.0) wvlg = atmdat.CaseBWlHeI[i];
					linadd( intens, wvlg, CaBLines[nelem][i].label, 'i', "Case B intensity " );
				}
			}
		}
	}


	if( iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max >= 4 &&
		( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max + iso_sp[ipH_LIKE][ipHYDROGEN].nCollapsed_max ) >=8 )
	{
		t_iso_sp* sp = &iso_sp[ipHE_LIKE][ipHELIUM];
		const long ipHe4s3S = iso_sp[ipHE_LIKE][ipHELIUM].QuantumNumbers2Index[4][0][3];
		const long ipHe4p3P = iso_sp[ipHE_LIKE][ipHELIUM].QuantumNumbers2Index[4][1][3];

		/* For info only, add the total photon flux in 3889 and 7065,
		* and 3188, 4713, and 5876. */
		double photons_3889_plus_7065 =
			/* to 2p3P2 */
			phots( sp->trans(ipHe3s3S,ipHe2p3P2) ) +
			phots( sp->trans(ipHe3d3D,ipHe2p3P2) ) +
			phots( sp->trans(ipHe4s3S,ipHe2p3P2) ) +
			/* to 2p3P1 */
			phots( sp->trans(ipHe3s3S,ipHe2p3P1) ) +
			phots( sp->trans(ipHe3d3D,ipHe2p3P1) ) +
			phots( sp->trans(ipHe4s3S,ipHe2p3P1) ) +
			/* to 2p3P0 */
			phots( sp->trans(ipHe3s3S,ipHe2p3P0) ) +
			phots( sp->trans(ipHe3d3D,ipHe2p3P0) ) +
			phots( sp->trans(ipHe4s3S,ipHe2p3P0) ) +
			/* to 2s3S */
			phots( sp->trans(ipHe3p3P,ipHe2s3S) ) +
			phots( sp->trans(ipHe4p3P,ipHe2s3S) ) ;

		long upperIndexofH8 = iso_sp[ipH_LIKE][ipHYDROGEN].QuantumNumbers2Index[8][1][2];

		/* Add in photon flux of H8 3889 */
		photons_3889_plus_7065 += 
			phots( iso_sp[ipH_LIKE][ipHYDROGEN].trans(upperIndexofH8,1) );

		/* now multiply by ergs of normalization line, so that relative flux of
		* this line will be ratio of photon fluxes. */
		if( LineSave.WavLNorm > 0 )
			photons_3889_plus_7065 *= (ERG1CM*1.e8)/LineSave.WavLNorm;
		linadd( photons_3889_plus_7065, 3889., "Pho+", 'i',
			"photon sum given in Porter et al. 2007 (astro-ph/0611579)");
	}

	/* ====================================================
	 * end helium
	 * ====================================================*/

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_helium returns\n" );
	}
	return;
}

// Searches through for next non-comment line, returns NULL on failure
inline char* read_data_line( char *chLine, int size, FILE *ioDATA )
{
	char *b;
	do
	{
		b = read_whole_line( chLine , size , ioDATA );
		if ( b == NULL)
			break;
	}
	while (chLine[0] == '#');
	return b;
}

STATIC void GetStandardHeLines(void)
{
	FILE *ioDATA;
	bool lgEOL;
	long i, i1, i2;

#	define chLine_LENGTH 1000
	char chLine[chLine_LENGTH];

	DEBUG_ENTRY( "GetStandardHeLines()" );

	if( trace.lgTrace )
		fprintf( ioQQQ," lines_helium opening he1_case_b.dat\n");

	ioDATA = open_data( "he1_case_b.dat", "r" );

	/* check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " lines_helium could not read first line of he1_case_b.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	/* first number is magic number, second is number of lines in file	*/
	i1 = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	i2 = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	NumLines = i2;

	/* the following is to check the numbers that appear at the start of he1_case_b.dat */
	if( i1 !=CASEBMAGIC )
	{
		fprintf( ioQQQ, 
			" lines_helium: the version of he1_case_ab.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			" lines_helium: I expected to find the number %i and got %li instead.\n" ,
			CASEBMAGIC, i1 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* get the array of densities */
	if ( read_data_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " lines_helium could not find line of densities in he1_case_ab.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}		
	
	lgEOL = false;
	i = 1;
	for( long j=0; !lgEOL && j < NUMDENS; ++j)
	{
		CaBDensities[j] = FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	}

	/* get the array of temperatures */
	if ( read_data_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " lines_helium could not find line of temperatures in he1_case_ab.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	lgEOL = false;
	i = 1;
	for (long j=0; !lgEOL && j < NUMTEMPS; ++j)
	{
		CaBTemps[j] = FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	}

	/* create space for array of values, if not already done */
	CaBIntensity = (double ****)MALLOC(sizeof(double ***)*(unsigned)LIMELM );
	/* create space for array of values, if not already done */
	CaBLines = (stdLines **)MALLOC(sizeof(stdLines *)*(unsigned)LIMELM );

	for( long nelem=ipHELIUM; nelem<LIMELM; ++nelem )
	{
		/** \todo	2	- this structure is currently only used for helium itself...
		 * stuff numbers in for other elements, or drop the [nelem] dimension off
		 * of CaBLines	*/
		if( nelem != ipHELIUM )
			continue;

		/* only grab core for elements that are turned on */
		if( nelem == ipHELIUM || dense.lgElmtOn[nelem])
		{
			/* create space for array of values, if not already done */
			CaBIntensity[nelem] = (double ***)MALLOC(sizeof(double **)*(unsigned)(i2) );
			CaBLines[nelem] = (stdLines *)MALLOC(sizeof(stdLines )*(unsigned)(i2) );

			/* avoid allocating 0 bytes, some OS return NULL pointer, PvH 
			CaBIntensity[nelem][0] = NULL;*/
			for( long j = 0; j < i2; ++j )
			{
				CaBIntensity[nelem][j] = (double **)MALLOC(sizeof(double*)*(unsigned)NUMDENS );

				CaBLines[nelem][j].ipHi = -1;
				CaBLines[nelem][j].ipLo = -1;
				strcpy( CaBLines[nelem][j].label , "    " );

				for( long k = 0; k < NUMDENS; ++k )
				{
					CaBIntensity[nelem][j][k] = (double *)MALLOC(sizeof(double)*(unsigned)NUMTEMPS );
					for( long l = 0; l < NUMTEMPS; ++l )
					{
						CaBIntensity[nelem][j][k][l] = 0.;
					}
				}
			}
		}
	}

	/* now read in the case B line data */
	long nelem = ipHELIUM;
	int line = 0;
	while( read_data_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* only look at lines without '#' in first col */
		if( (chLine[0] == ' ') || (chLine[0]=='\n') )
			break;
		
		/* get lower and upper level index */
		long j = 1;
		// the first number is the wavelength, which is only used if the
		// model atom is too small to include this transition
		realnum wl = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
		long ipLo = (long)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
		long ipHi = (long)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
		CaBLines[nelem][line].ipLo = ipLo;
		CaBLines[nelem][line].ipHi = ipHi;
		
		ASSERT( CaBLines[nelem][line].ipLo < CaBLines[nelem][line].ipHi );
		
		strcpy( CaBLines[nelem][line].label, "Ca B" );
		
		t_iso_sp* sp = &iso_sp[ipHE_LIKE][nelem];
		if( ipHi < sp->numLevels_max - sp->nCollapsed_max )
			atmdat.CaseBWlHeI.push_back( sp->trans(ipHi,ipLo).WLAng() );
		else
			atmdat.CaseBWlHeI.push_back( wl );
		
		for( long ipDens = 0; ipDens < NUMDENS; ++ipDens )
		{
			if( read_whole_line( chLine, (int)sizeof(chLine) , ioDATA ) == NULL )
			{
				fprintf( ioQQQ, " lines_helium could not scan case B lines, current indices: %li %li\n",
							CaBLines[nelem][line].ipHi,
							CaBLines[nelem][line].ipLo );
				cdEXIT(EXIT_FAILURE);
			}
			
			char *chTemp = chLine;
			j = 1;
			long den = (long)FFmtRead(chTemp,&j,sizeof(chTemp),&lgEOL);
			if( den != ipDens + 1 )
				TotalInsanity();
			for( long ipTe = 0; ipTe < NUMTEMPS; ++ipTe )
			{
				double b;
				if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
				{
					fprintf( ioQQQ, " lines_helium could not scan case B lines, current indices: %li %li\n",
								CaBLines[nelem][line].ipHi,
								CaBLines[nelem][line].ipLo );
					cdEXIT(EXIT_FAILURE);
				}
				++chTemp;
				sscanf( chTemp, "%le" , &b );
				CaBIntensity[nelem][line][ipDens][ipTe] = b;
			}
		}
		line++;
	}

	ASSERT( line == NumLines );
	ASSERT( atmdat.CaseBWlHeI.size() == (unsigned)line );

	/* close the data file */
	fclose( ioDATA );
	return;
}

/** \todo	there is a virtually identical routine in helike_recom.cpp -> combine */
STATIC double TempInterp2( double* TempArray , double* ValueArray, long NumElements, double Te )
{
	long int ipTe=-1;
	double rate = 0.;
	long i0;

	DEBUG_ENTRY( "TempInterp2()" );

	/* te totally unknown */
	if( Te <= TempArray[0] )
	{
		return ValueArray[0] + log10( TempArray[0]/Te );
	}
	else if( Te >= TempArray[NumElements-1] )
	{
		return ValueArray[NumElements-1];
	}

	/* now search for temperature */
	ipTe = hunt_bisect( TempArray, NumElements, Te );			

	ASSERT( (ipTe >=0) && (ipTe < NumElements-1)  );

	/* Do a four-point interpolation */
	const int ORDER = 3; /* order of the fitting polynomial */
	i0 = max(min(ipTe-ORDER/2,NumElements-ORDER-1),0);
	rate = lagrange( &TempArray[i0], &ValueArray[i0], ORDER+1, Te );

	return rate;
}

/** \todo	2	say where these come from	*/	
/* For double-ionization discussions, see Lindsay, Rejoub, & Stebbings 2002	*/
/* Also read Itza-Ortiz, Godunov, Wang, and McGuire 2001.	*/
STATIC void DoSatelliteLines( long nelem )
{
	long ipISO = ipHE_LIKE;
	
	DEBUG_ENTRY( "DoSatelliteLines()" );

	ASSERT( iso_ctrl.lgDielRecom[ipISO] && dense.lgElmtOn[nelem] );

	for( long i=0; i < iso_sp[ipISO][nelem].numLevels_max; i++ )
	{
		double dr_rate = iso_sp[ipISO][nelem].fb[i].DielecRecomb;
		TransitionProxy tr = SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][i]];

		tr.Emis().xObsIntensity() = 
		tr.Emis().xIntensity() = 
			dr_rate * dense.eden * dense.xIonDense[nelem][nelem+1-ipISO] * ERG1CM * tr.EnergyWN();
		tr.Emis().pump() = 0.;

		PutLine( tr, "iso satellite line" );
	}

	return;
}
