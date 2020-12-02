/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_hydro put H-like iso sequence into line intensity stack */
#include "cddefines.h"
#include "atmdat.h"
#include "dense.h"
#include "prt.h"
#include "hydrogenic.h"
#include "iso.h"
#include "rfield.h"
#include "geometry.h"
#include "lines.h"
#include "phycon.h"
#include "radius.h"
#include "secondaries.h"
#include "trace.h"
#include "two_photon.h"
#include "lines_service.h"
#include "elementnames.h"

static bool isSkipTransSet = false;
static vector< vector< pair<long, long> > > skipTrans( LIMELM );

/* This loop is identical to the one in lines_hydro() that brings together
 * the Balmer series, and resets all but one of the fine structure lines
 * to a zero intensity. */
static void collectSkipTrans( void )
{
	long ipISO = ipH_LIKE;

	for( long nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.IonHigh[nelem] == nelem + 1 )
		{
			vector< pair<long, long> > thisElm_skipTrans;

			/* bring nL - n'L' emission together as n-n' emission. */
			for( long ipHi=1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
			{
				long index_of_nHi_P;

				/* is ipHi is collapsed level, index_of_nHi_P is ipHi */
				if( N_(ipHi) > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max )
					index_of_nHi_P = ipHi;
				else
					index_of_nHi_P = iso_sp[ipH_LIKE][nelem].QuantumNumbers2Index[ N_(ipHi) ][1][2];

				/* only need to consider resolved lower level here */
				for( long ipLo=0; ipLo < ipHi; ipLo++ )
				{	
					long index_of_nLo_S = iso_sp[ipH_LIKE][nelem].QuantumNumbers2Index[ N_(ipLo) ][0][2];

					/* jump out if ipLo is collapsed 
					 * NB this must be up to n_HighestResolved_local and not n_HighestResolved_max */
					if( N_(ipLo) > iso_sp[ipH_LIKE][nelem].n_HighestResolved_local || N_(ipLo) == N_(ipHi) )
						break;

					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					/* add everything into nP - n'S, skip if current indices are those levels. */
					if( ipHi == index_of_nHi_P && ipLo == index_of_nLo_S )
						continue;
					else
					{
						thisElm_skipTrans.push_back( pair<long, long>( ipLo, ipHi ) );
						if( 0 )
						{
							fprintf( ioQQQ,
								"%2s\t index_of_nHi_P= %2ld (%5s)\t index_of_nLo_S= %2ld (%5s)\t"
								" ipHi= %2ld (%5s)\t ipLo= %2ld (%5s)\t %10g %10g\t Dwl= %g\n",
								elementnames.chElementSym[ nelem ],
								index_of_nHi_P,
								iso_sp[ipISO][nelem].st[index_of_nHi_P].chConfig().c_str(),
								index_of_nLo_S,
								iso_sp[ipISO][nelem].st[index_of_nLo_S].chConfig().c_str(),
								ipHi,
								iso_sp[ipISO][nelem].st[ipHi].chConfig().c_str(),
								ipLo,
								iso_sp[ipISO][nelem].st[ipLo].chConfig().c_str(),
								iso_sp[ipISO][nelem].trans(index_of_nHi_P,index_of_nLo_S).WLAng(),
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLAng(),
								iso_sp[ipISO][nelem].trans(index_of_nHi_P,index_of_nLo_S).WLAng() -
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLAng() );
						}
					}
				}
			}
			skipTrans[ nelem ] = thisElm_skipTrans;
		}
	}

	isSkipTransSet = true;
}

static long this_ipLo, this_ipHi;
static bool isSkipped( const pair< long, long > level_indices )
{
	return level_indices.first == this_ipLo &&
		level_indices.second == this_ipHi;
}

static bool isTransSkipped( long nelem, long ipLo, long ipHi )
{
	this_ipLo = ipLo;
	this_ipHi = ipHi;

	vector< pair<long, long> >::iterator sti;
	sti = find_if( skipTrans[ nelem ].begin(), skipTrans[ nelem ].end(), isSkipped );

	bool skippedTrans = true;
	if( sti == skipTrans[ nelem ].end() )
		skippedTrans = false;
	return skippedTrans;
}

string iso_comment_tran_levels( long ipISO, long nelem, long ipLo, long ipHi )
{
	string isoSeq = "H-like, ";
	if( ipISO == ipHE_LIKE )
	        isoSeq = "He-like, ";

	// Terrible hack to get transitions into n=2 to not be
	// reported in the comments of 'save lines labels' as
	// into only 2^2S or 2^2P.
	if( ipISO == ipH_LIKE && (ipLo == 1 || ipLo == 2 ) )
	{
		return isoSeq + db_comment_tran_levels( ipLo+1, ipHi+1 ) + ", n=  2 - "
			+ (*iso_sp[ipISO][nelem].trans(ipHi,ipLo).Hi()).chConfig(); 
	}

	return isoSeq + db_comment_tran_levels( ipLo+1, ipHi+1 ) + ", "
		+ GenerateTransitionConfiguration( iso_sp[ipISO][nelem].trans(ipHi,ipLo) );
}

void lines_hydro(void)
{
	long ipISO = ipH_LIKE;
	long int i, nelem, ipHi, ipLo;
	string chLabel="    ";

	double hbetab, 
		em , 
		caseb;

	DEBUG_ENTRY( "lines_hydro()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "   lines_hydro called\n" );

	// this can be changed with the atom levels command but must be at least 3
	ASSERT( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max >= 3 );
	ASSERT( iso_sp[ipH_LIKE][ipHELIUM].n_HighestResolved_max >= 3 );

	i = StuffComment( "H-like iso-sequence" );
	linadd( 0., (realnum)i , "####", 'i',
		" start H -like iso sequence ");

	/*fprintf(ioQQQ," debugg\t%.2e\t%.2e\t%.2e\n", 
		radius.drad,
		iso_sp[ipH_LIKE][ipHYDROGEN].xLineTotCool , 
		iso_sp[ipH_LIKE][ipHYDROGEN].cLya_cool);*/

	/* >>chng 95 jun 25 changed from info to cooling to pick this up in primal.in   */
	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].cLya_cool),1215.67,"Cool",'i',
		"collisionally excited La cooling ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].cLya_cool),1215.67,"Heat",'i',
		"  collisionally de-excited La heating ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].cLyrest_cool),960,"Crst",'i',
		"  cooling due to n>2 Lyman lines ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].cLyrest_cool),960,"Hrst",'i',
		"  heating due to n>2 Lyman lines ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].cBal_cool),4861.33,"Crst",'i',
		"  cooling due to n>3 Balmer lines ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].cBal_cool),4861.33,"Hrst",'i',
		"  heating due to n>3 Balmer lines ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].cRest_cool),0,"Crst",'i',
		"  cooling due to higher Paschen lines ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].cRest_cool),0,"Hrst",'i',
		"  heating due to higher Paschen lines ");

	/* remember largest fractional ionization of H due to secondaries */
	secondaries.SecHIonMax = MAX2( secondaries.SecHIonMax , secondaries.sec2total );

	/* remember fraction of H ionizations due to ct */
	atmdat.HIonFracMax = MAX2( atmdat.HIonFracMax, atmdat.HIonFrac);

	/* remember largest fraction of thermal collisional ionization of H ground state */
	hydro.HCollIonMax = 
		(realnum)MAX2( hydro.HCollIonMax , hydro.H_ion_frac_collis );

	linadd(secondaries.x12tot*iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop()*1.634e-11,1215.67,"LA X" ,'i',
		"Lya contribution from suprathermal secondaries from ground ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].coll_ion),0,"CION",'c',
		"collision ionization cooling of hydrogen ");

	linadd(MAX2(-iso_sp[ipH_LIKE][ipHYDROGEN].coll_ion,0.),0,"3bHt",'h',
		"  this is the heating due to 3-body recombination ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHELIUM].coll_ion),0,"He2C",'c',
		"collision ionization cooling of He+ ");

	linadd(MAX2(-iso_sp[ipH_LIKE][ipHELIUM].coll_ion,0.),0,"He2H",'h',
		"  this is the heating due to 3-body recombination onto He+");

	fixit("why is there a zero here?");
	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop()*0.*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ipH2p][ipH1s].pestrk*1.634e-11,1215.67,"Strk",'i',
	  "  Stark broadening contribution to line ");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH3s].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ipH3s][ipH2p].pestrk*3.025e-12,
	  6562.81,"Strk",'i',
	  "  Stark broadening contribution to line ");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH4s].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ipH4s][ipH2p].pestrk*4.084e-12,
	  4861.33,"Strk",'i',
	  "Stark broadening contribution to line ");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH4p].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ipH4p][ipH3s].pestrk*1.059e-12,
	  18751,"Strk",'i',
		   " Stark broadening contribution to line ");

	/* pestrk[5,4] is A[4,5]*pest[4,5] 
	 * Stark broadening contribution to line */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max >= 5 )
	{
		long ip5p = iso_sp[ipH_LIKE][ipHYDROGEN].QuantumNumbers2Index[5][1][2];
		linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ip5p].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ip5p][ipH4s].pestrk*4.900e-13,40512,"Strk",'i',
			"Stark broadening part of line");
	}
	/* this can fail if RT_line_all never updates the ots rates, a logic error,
	 * but only assert this during actual calculation (ipass>0), */
	ASSERT( LineSave.ipass  <1 ||
		iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().ots()>= 0.);

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).EnergyErg(), 1215.67,"Dest",'i',
		"  portion of line lost due to absorp by background opacity ");

	/* portion of line lost due to absorb by background opacity */
	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH3p,ipH2s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH3p,ipH2s).EnergyErg(), 6562.81,"Dest",'i',
		"Ha destroyed by background opacity");

	/* portion of line lost due to absorp by background opacity */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max >= 5 )
	{
		long ip5p = iso_sp[ipH_LIKE][ipHYDROGEN].QuantumNumbers2Index[5][1][2];
		linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip5p,ipH4s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip5p,ipH4s).EnergyErg(),40516, "Dest",'i',
			"portion of line lost due to absorb by background opacity");
	}

	/* portion of line lost due to absorb by background opacity */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_max > ipH4p )
		linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH2s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH2s).EnergyErg(), 4861.33,"Dest",'i',
			"portion of line lost due to absorb by background opacity");

	/* portion of line lost due to absorb by background opacity */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_max > ipH4p )
		linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH3s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH3s).EnergyErg() ,18751, "Dest",'i',
			"portion of line lost due to absorb by background opacity");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Aul()*
		hydro.dstfe2lya*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).EnergyErg() , 1215.67 , "Fe 2" , 'i',
		"Ly-alpha destroyed by overlap with FeII " );

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].RadRec_caseB*dense.xIonDense[ipHYDROGEN][1]*dense.eden * 1.64e-11,1215.67,"Ca B",'i',
		" simple high-density case b intensity of Ly-alpha, no two photon ");

	/* these entries only work correctly if the APERTURE command is not in effect */
	if( geometry.iEmissPower == 2 )
	{
		/* H-beta computed from Q(H) and specified covering factor */
		if( nzone == 1 )
		{
			/* evaluate the case b emissivity by interpolating on the hummer & storey tables */
			caseb = rfield.qhtot*
				atmdat_HS_caseB( 4 , 2 , 1 , phycon.te , dense.eden, 'b' ) / iso_sp[ipH_LIKE][ipHYDROGEN].RadRec_caseB;
			/* the atmdat_HS_caseB returned -1 if the physical conditions were outside range of validity.  
			 * In this case use simple approximation with no temperature or density dependence */
			if( caseb < 0 )
			{
				caseb = rfield.qhtot*4.75e-13;
			}
			LineSave.lines[LineSave.nsum].SumLineZero();
		}
		else
		{
			caseb = 0.;
		}
		/* H-beta computed from Q(H) and specified covering factor */
		linadd( caseb/radius.dVeffAper*geometry.covgeo , 4861.33 , "Q(H)" , 'i' ,
			"Case B H-beta computed from Q(H) and specified covering factor");

		if( nzone == 1 )
		{
			// the cast to double prevents an FPE with Solaris Studio 12.4 (limit_compton_hi_t.in)
			caseb = rfield.qhtot*double(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).EnergyErg());
			LineSave.lines[LineSave.nsum].SumLineZero();
		}
		else
		{
			caseb = 0.;
		}
		/* >>chng 02 nov 05, better approximation for Lya for temperature of first zone */
		linadd( caseb/radius.dVeffAper*geometry.covgeo , 1215.67 , "Q(H)" , 'i',
			"Ly-alpha from Q(H), high-dens lim, specified covering factor" );
	}

	/* this is the main printout, where line intensities are entered into the stack */
	for( nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			ASSERT( iso_sp[ipH_LIKE][nelem].n_HighestResolved_max >= 3 );

			for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
			{
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					set_xIntensity( iso_sp[ipISO][nelem].trans(ipHi,ipLo) );
				}
			}
		}
	}

	if( ! isSkipTransSet )
	{
		collectSkipTrans();
	}

	/* create emissivity or intensity for hydrogenic species,
	 * first combine/bring balmer series together */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.IonHigh[nelem] == nelem + 1 )
		{
			/* bring nL - n'L' emission together as n-n' emission. */
			for( ipHi=1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
			{
				long index_of_nHi_P;

				/* is ipHi is collapsed level, index_of_nHi_P is ipHi */
				if( N_(ipHi) > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max )
					index_of_nHi_P = ipHi;
				else
					index_of_nHi_P = iso_sp[ipH_LIKE][nelem].QuantumNumbers2Index[ N_(ipHi) ][1][2];

				/* only need to consider resolved lower level here */
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{	
					long index_of_nLo_S = iso_sp[ipH_LIKE][nelem].QuantumNumbers2Index[ N_(ipLo) ][0][2];

					/* jump out if ipLo is collapsed 
					 * NB this must be up to n_HighestResolved_local and not n_HighestResolved_max */
					if( N_(ipLo) > iso_sp[ipH_LIKE][nelem].n_HighestResolved_local || N_(ipLo) == N_(ipHi) )
						break;

					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					/* add everything into nP - n'S, skip if current indices are those levels. */
					if( ipHi == index_of_nHi_P && ipLo == index_of_nLo_S )
						continue;
					else
					{
						/* add resolved line to nP - n'S */
						iso_sp[ipH_LIKE][nelem].trans(index_of_nHi_P,index_of_nLo_S).Emis().xIntensity() +=
							iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().xIntensity();
						iso_sp[ipH_LIKE][nelem].trans(index_of_nHi_P,index_of_nLo_S).Emis().xObsIntensity() +=
							iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().xObsIntensity();
						/* zero out the resolved line */
						iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().xIntensity() = 0;
						iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().xObsIntensity() = 0;
						//ASSERT( iso_sp[ipH_LIKE][nelem].trans(index_of_nHi_P,index_of_nLo_S).Emis().xIntensity() >= 0. );
						//ASSERT( iso_sp[ipH_LIKE][nelem].trans(index_of_nHi_P,index_of_nLo_S).Emis().xObsIntensity() >= 0. );
					}
				}
			}
		}
	}

	/* H beta recombination, assuming old case B */
	hbetab = (double)((exp10(-20.89 - 0.10612*POW2(phycon.alogte - 4.4)))/
	  phycon.te);
	/* need to pass this assert if CaBo is to have valid array indices for ipCont */
	/* 06 aug 28, from numLevels_max to _local. */
	/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
	ASSERT( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_max > 4 );
	hbetab *= dense.xIonDense[ipHYDROGEN][1]*dense.eden;

	lindst(hbetab, -4861.33 ,"CaBo",
		1 ,'i',false,
		" this is old case b based on Ferland (1980) PASP ");

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* need to pass this assert if CaBo is to have valid array indices for ipCont */
		/* 06 aug 28, from numLevels_max to _local. */
		/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
		ASSERT( iso_sp[ipH_LIKE][ipHELIUM].numLevels_max > 4 );
		/* 1640 1640 1640 */
		em = 2.03e-20/(phycon.te70*phycon.te10*phycon.te03);
		em *= dense.xIonDense[ipHELIUM][2]*dense.eden;

		lindst(em,-1640,"CaBo",
			1,'i',false,
			" old prediction of He II 1640, Case B at low densities");

		/* hydrogenic helium */
		/* old prediction of He II 4686, case B */
		em = 2.52e-20/(pow(phycon.te,1.05881));
		em *= dense.xIonDense[ipHELIUM][2]*dense.eden;

		lindst(em,-4685.64,"CaBo",	1,'i',false,
			   " old prediction of He II 4686, Case B at low densities");
	}

	/* predict case b intensities of hydrogen lines */
	if( LineSave.ipass <= 0 )
	{
		for(nelem=0; nelem<HS_NZ; ++nelem )
		{
			atmdat.lgHCaseBOK[0][nelem] = true;
			atmdat.lgHCaseBOK[1][nelem] = true;
		}
	}
	/* this is the main printout, where line intensities are entered into the stack */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* HS_NZ is limit to charge of elements in HS predictions, now 8 == oxygen */
			/* but don't do the minor elements, Li, Be, B - these were not read in and so should not be
			 * printed - remove equivalent if statement in atmdat_readin.cpp to read them in */
			if( nelem < HS_NZ && (nelem<2 || nelem>4) )
			{
				int iCase;
				for( iCase=0; iCase<2; ++iCase )
				{
					char chAB[2]={'A','B'};
					char chLab[5]="Ca  ";

					/* adding iCase means start from n=1 for case A, n=2 for Case B,
					 * note that principal quantum number is on physics scale, not C */
					/* 06 aug 28, both of these from numLevels_max to _local. */
					/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
					/* HS data file extends to ipHi == atmdat.ncut[iCase][nelem]==25, report all possible values for comparisons */
					for( ipLo=1+iCase; ipLo<MIN2(atmdat.ncut[iCase][nelem]-1,iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + iso_sp[ipH_LIKE][nelem].nCollapsed_max); ++ipLo )
					{
						for( ipHi=ipLo+1; ipHi< MIN2(atmdat.ncut[iCase][nelem],iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + iso_sp[ipH_LIKE][nelem].nCollapsed_max+1); ++ipHi )
						{
							realnum wl;
							double case_b_Intensity;
							long int ipCHi , ipCLo;
							/* Put case b predictions into line stack
							 * NB NB NB each Hummer & Storey case b line must be 
							 * explicitly clobbered by hand in routine final if 
							 * atmdat.lgHCaseBOK[iCase][nelem] flag is set false
							 * since this indicates that we exceeded bounds of table,
							 * DO NOT want to print lines in that case */

							/* first do case b emissivity of balmer lines */

							/* get HS predictions */
							case_b_Intensity = atmdat_HS_caseB( ipHi,ipLo , nelem+1, phycon.te , dense.eden, chAB[iCase] );
							if( case_b_Intensity<=0. )
							{
								atmdat.lgHCaseBOK[iCase][nelem] = false;
								case_b_Intensity = 0.;
							}

							case_b_Intensity *= dense.xIonDense[nelem][nelem+1-ipISO]*dense.eden;

							if( iCase==0 && ipLo==1 )
							{
								/* get physical scal prin quant numbers onto cloudy c scale */
								ipCHi = ipHi;
								ipCLo = 0;
							}
							else
							{
								/* get physical scal prin quant numbers onto cloudy c scale */
								ipCHi = ipHi;
								ipCLo = ipLo;
							}

							/* make label either Ca A or Ca B */
							chLab[3] = chAB[iCase];

							/* new treatment is different from old for indices greater than 2. */
							if( ipCHi > 2 )
							{
								if( ipCLo >  2 )
								{
									/* if both indices above two, just treat as nP to n'S transition. */
									ipCHi = iso_sp[ipH_LIKE][nelem].QuantumNumbers2Index[ipCHi][1][2];
									ipCLo = iso_sp[ipH_LIKE][nelem].QuantumNumbers2Index[ipCLo][0][2];
								}
								else if(  ipCLo ==  2 )
								{
									/* treat as nS to 2P transition. */
									ipCHi = iso_sp[ipH_LIKE][nelem].QuantumNumbers2Index[ipCHi][0][2];
								}
								else if( ipCLo ==  1 || ipCLo == 0 )
								{
									/* treat as nP to n'S transition. */
									ipCHi = iso_sp[ipH_LIKE][nelem].QuantumNumbers2Index[ipCHi][1][2];
								}
							}

							/* this is wavelength of interpolated case b from HS tables */
							wl = iso_sp[ipH_LIKE][nelem].trans(ipCHi,ipCLo).WLAng();
							atmdat.WaveLengthCaseB[nelem][ipHi][ipLo] = wl;

							lindst(case_b_Intensity,wl,chLab,iso_sp[ipH_LIKE][nelem].trans(ipCHi,ipCLo).ipCont(),'i',false,
								" case a or case b from Hummer & Storey tables" );
						}
					}
				}
			}

			// add two-photon details here
			if( LineSave.ipass == 0 )
			{
				/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" 
				 * the result, chLable, is only used when ipass == 0, can be undefined otherwise */
				chLabel = chIonLbl(nelem+1, nelem+1-ipISO);
			}
			for( vector<two_photon>::iterator tnu = iso_sp[ipH_LIKE][nelem].TwoNu.begin(); tnu != iso_sp[ipH_LIKE][nelem].TwoNu.end(); ++tnu )
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
					2. * wn2ang( iso_sp[ipH_LIKE][nelem].trans( (*tnu).ipHi, (*tnu).ipLo ).EnergyWN() ),
					chLabel.c_str(), 'r', tpc_comment.c_str() );
			}

			/* NB NB - low and high must be in this order so that all balmer, paschen,
			 * etc series line up correctly in final printout */
			for( ipLo=ipH1s; ipLo < iso_sp[ipH_LIKE][nelem].numLevels_max-1; ipLo++ )
			{
				/* don't bother with decays to 2p since we set them to zero above */
				if( ipLo==ipH2p )
					continue;

				/* set number of levels we want to print, first is default,
				 * only print real levels, second is set with "print line
				 * iso collapsed" command */
				long int nLoop  = iso_sp[ipH_LIKE][nelem].numLevels_max - iso_sp[ipH_LIKE][nelem].nCollapsed_max;
				if( prt.lgPrnIsoCollapsed )
					nLoop  = iso_sp[ipH_LIKE][nelem].numLevels_max;

				for( ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					// skip non-radiative transitions
					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).ipCont() < 1 )
						continue;

					// skip 2s-1s, so that 2p-1s comes first and cdLine finds LyA instead of the M2 transition.	
					if( ipHi==1 && ipLo==0 )
						continue;

					if( isTransSkipped( nelem, ipLo, ipHi ) )
						continue;

					string comment_trans = "";
					if( LineSave.ipass == 0 )
					{
						comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
					}
					PutLine(iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo), comment_trans.c_str());
				}
			}
		}
	}
	
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_hydro returns\n" );
	}
	return;
}

