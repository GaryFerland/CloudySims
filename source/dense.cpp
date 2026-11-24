/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"

#include "dense.h"

#include "abund.h"
#include "colden.h"
#include "conv.h"
#include "dynamics.h"
#include "elementnames.h"
#include "deuterium.h"
#include "hmi.h"
#include "phycon.h"
#include "radius.h"
#include "struc.h"
#include "thermal.h"
#include "trace.h"
#include "iso.h"
#include "h2.h"
#include "mole.h"
#include "save.h"

void t_dense::updateXMolecules()
{
	total_molecule_elems(m_xMolecules);
}

void t_dense::zero()
{
	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		dense.SetGasPhaseDensity( nelem, 0. );
		for( long ion=0; ion < LIMELM+1; ion++ )
		{
			dense.xIonDense[nelem][ion] = 0.;
		}
	}
	for (long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem)
		m_xMolecules[nelem] = 0.f;
	// needed for TempChange to work but arrays needed for EdenChange to
	// work are not yet defined
	dense.eden = 1.;

	dense.xMassTotal = 0.;	
	dense.xNucleiTotal = 1.;
	/* WJH */
	dense.xMassDensity0 = -1.0f; 

	/* this is a scale factor that changes the n(H0)*1.7e-4 that is added to the
	* electron density to account for collisions with atomic H.  it is an order of
	* magnitude guess, so this command provides ability to see whether it affects results */
	dense.HCorrFac = 1.f;
}

void ScaleAllDensities(realnum factor)
{
	double edensave = dense.eden;

	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		if (dense.lgElmtOn[nelem])
		{
			ScaleIonDensities( nelem, factor );
			dense.SetGasPhaseDensity( nelem, dense.gas_phase[nelem] * factor);
		}
	}

	EdenChange( dense.eden * factor );
	
	if( trace.lgTrace && trace.lgNeBug )
	{
		fprintf( ioQQQ, " EDEN change PressureChange from to %10.3e %10.3e %10.3e\n", 
					edensave, dense.eden, edensave/dense.eden );
	}

	hmi.H2_total *= factor;
	h2.ortho_density *= factor;
	h2.para_density *= factor;
	for( long mol=0; mol < mole_global.num_total; mol++ )
	{
		mole.species[mol].den *= factor;
	}
	dense.updateXMolecules();
	deut.updateXMolecules();

	ASSERT(lgElemsConserved());
}

void ScaleIonDensities( const long nelem, const realnum factor )
{
	double renorm;
	for( long ion=0; ion < (nelem + 2); ion++ )
	{
		dense.xIonDense[nelem][ion] *= factor;
		if (nelem-ion >= 0 && nelem-ion < NISO)
			iso_renorm(nelem, nelem-ion, renorm);
	}
	
	if( nelem==ipHYDROGEN && deut.lgElmtOn )
		ScaleDensitiesDeuterium( factor );

	return;
}

void t_dense::SetGasPhaseDensity( const long nelem, const realnum density )
{
	gas_phase[nelem] = density;

	if( nelem==ipHYDROGEN && deut.lgElmtOn )
	{
		// fractionation is taken care of inside called routine
		SetGasPhaseDeuterium( density );
	}

	return;
}

bool lgElemsConserved (void)
{
	bool lgOK=true;

	/* this checks all molecules and ions add up to abundance */
	for( ChemNuclideList::iterator nuc = nuclide_list.begin(); nuc != nuclide_list.end(); ++nuc)
	{
		long nelem = (*nuc)->el()->Z-1;
		
		double sum_monatomic = 0.;
		realnum sum_molecules = 0.f;
		realnum sum_gas_phase = 0.f;

		if( (*nuc)->label() == "D" )
		{
			sum_monatomic = deut.xIonDense[0] + deut.xIonDense[1];
			sum_molecules = deut.xMolecules();
			sum_gas_phase = deut.gas_phase;

		}
		/* check that species add up if element is turned on */
		else if( dense.lgElmtOn[nelem] )
		{
			/* this sum of over the molecules did not include the atom and first
			 * ion that is calculated in the molecular solver */
			for( long ion=0; ion<nelem+2; ++ion )
			{
				sum_monatomic += dense.xIonDense[nelem][ion] * (*nuc)->frac;
			}
			sum_molecules = dense.xMolecules(nelem) * (*nuc)->frac;
			sum_gas_phase = dense.gas_phase[nelem] * (*nuc)->frac;
		}
		else
			continue;


		if( sum_monatomic + sum_molecules <= SMALLFLOAT && 
			sum_gas_phase > SMALLFLOAT)
		{
			fprintf(ioQQQ,"PROBLEM non-conservation of nuclei %s\tions %g moles %g error %g of %g\n",
				(*nuc)->label().c_str(),
				sum_monatomic,
				sum_molecules,
				sum_monatomic + sum_molecules-sum_gas_phase,
				sum_gas_phase );
			lgOK=false;
		}

		/* check that sum agrees with total gas phase abundance */
		/** \todo	2	this assert is not passed if error made much smaller.  This
		 * error should be related to a check on convergence of the molecular
		 * networks and their sum rules, with a criteria used here and there */
		if( fabs( sum_monatomic + sum_molecules - sum_gas_phase ) > 
			conv.GasPhaseAbundErrorAllowed * sum_gas_phase )
		{
			fprintf(ioQQQ,"PROBLEM non-conservation of nuclei %s\t nzone %li atoms %.12e moles %.12e sum %.12e tot gas %.12e rel err %.3e\n",
				(*nuc)->label().c_str(),
				nzone,
				sum_monatomic,
				sum_molecules,
				sum_monatomic + sum_molecules,
				sum_gas_phase,
				(sum_monatomic + sum_molecules - sum_gas_phase)/sum_gas_phase );
			fprintf(ioQQQ," Consider increasing value of conv.GasPhaseAbundErrorAllowed with SET DENSITY TOLERANCE as a workaround. "
					"Its current value is %.2e\n" , conv.GasPhaseAbundErrorAllowed );
			fprintf(ioQQQ," Temperature %.2e H density %.2e zone %ld \n", phycon.te, dense.gas_phase[ipHYDROGEN], nzone);
			lgOK=false;
		}
		
		if( 0 )
		{
			// print reactions involving the neutral 
			mole_print_species_reactions( findspecies( elementnames.chElementSym[nelem] ) );
		}
	}

	return lgOK;
}

void lgStatesConserved ( long nelem, long ionStage, qList states, long numStates, realnum err_tol, long loop_ion )
{
	if( 0 && conv.lgConvIoniz() && 
		dense.lgElmtOn[nelem] && 
		dense.xIonDense[nelem][ionStage]/dense.eden > conv.EdenErrorAllowed/10. &&
		dense.xIonDense[nelem][ionStage] > SMALLFLOAT )
	{
		double abund = 0.;
		for (long ipLev=0; ipLev<numStates; ipLev++)
		{
			abund += states[ipLev].Pop();
		}

		double rel_err = fabs(abund-dense.xIonDense[nelem][ionStage])/
			SDIV(dense.xIonDense[nelem][ionStage]);

		//fprintf(ioQQQ,"Iso consistency error %ld %ld %g\n",nelem,ionStage,rel_err);

		if( rel_err > err_tol ) 
		{
			/* we'll set the flag for non-convergence, but don't bother complaining in printout for first sweeps */
			if( 0 && nzone > 0 && loop_ion > 0 )
			{	
				fprintf(ioQQQ,"PROBLEM Inconsistent states/stage pops nzone %3ld loop_ion %2ld nelem %2ld ion %2ld states = %e stage = %e error = %e\n",
					nzone,
					loop_ion,
					nelem,
					ionStage,
					abund,
					dense.xIonDense[nelem][ionStage],
					rel_err );
			}
			ostringstream chConvIoniz;
			chConvIoniz << "States!=stage pops nelem " << nelem << " ion " << ionStage;
			conv.setConvIonizFail( chConvIoniz.str().c_str(), abund,
					       dense.xIonDense[nelem][ionStage] );
		}
	}
}

void SumDensities( void )
{
	realnum den_mole,
		DenseAtomsIons;

	/* find total number of atoms and ions
	 * at this point does not include molecules */
	DenseAtomsIons = 0.;
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem]  )
		{
			for( long ion=0; ion<=nelem+1; ++ion )
				DenseAtomsIons += dense.xIonDense[nelem][ion];
		}
	}

	/* find number of molecules of the heavy elements in the gas phase.
	 * Atoms and ions were counted above.  Do not include ices, solids
	 * on grains */
	den_mole = total_molecules_gasphase();
	
	/* total number of atoms, ions, and molecules in gas phase per unit vol */
	dense.xNucleiTotal = den_mole + DenseAtomsIons;
	if( dense.xNucleiTotal > BIGFLOAT )
	{
		fprintf(ioQQQ, "PROBLEM DISASTER SumDensities has found "
			"dense.xNucleiTotal with an insane density.\n");
		fprintf(ioQQQ, "The density was %.2e\n",
			dense.xNucleiTotal);
		TotalInsanity();
	}
	ASSERT( dense.xNucleiTotal > 0. );

	/* particle density that enters into the pressure includes electrons
	 * cm-3 */
	dense.pden = (realnum)(dense.eden + dense.xNucleiTotal);
        
	/* dense.wmole is molecular weight, AtomicMassUnit per particle */
	dense.wmole = 0.;
	for( long i=0; i < LIMELM; i++ )
	{
		dense.wmole += dense.gas_phase[i]*(realnum)dense.AtomicWeight[i];
	}
        /* dense.wmole is now molecular weight, AtomicMassUnit per particle */
        dense.wmole /= dense.pden;
        ASSERT( dense.wmole > 0. && dense.pden > 0. );
        
        /* xMassDensity is density in grams per cc */
        /** \todo       2       - should this include mass in grains? */
        /** \todo       2       - should this include mass in grain mantle ice deposits? */
        dense.xMassDensity = (realnum)(dense.wmole*ATOMIC_MASS_UNIT*dense.pden);
        
        /*>>chng 04 may 25, introduce following test on xMassDensity0 
         * WJH 21 may 04, this is the mass density that corresponds to the hden 
         * specified in the init file. It is used to calculate the mass flux in the 
         * dynamics models. It may not necessarily be the same as struc.DenMass[0],
         * since the pressure solver may have jumped to a different density at the
         * illuminated face from that specified.*/   
        if( dense.xMassDensity0 < 0.0 )
                 dense.xMassDensity0 = dense.xMassDensity;

	return;
}

namespace {
	/** >>chng 26 nov 21, added helper function for fixing an element table bug
	 * 1. Scales all molecules by the factor of their most restrictive constituent
	 * 2. Refreshes cached molecular totals via dense.updateXMolecules()
	 * 3. Ionic species take the remainder, being scaled for conservation purposes
	 *
	 * Anything calling RebalanceAbundances mustn't first modify xIonDense 
	 **/
	void RebalanceAbundances( const double abund_ratio[] ) 
	{
		DEBUG_ENTRY( "RebalanceAbundances()" );

		ASSERT( !deut.lgElmtOn );  // molecules containing deuterium not handled here

		bool lgModified[LIMELM];

		for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
			lgModified[nelem] = ( abund_ratio[nelem] != 1.0 );

		for( long i=0; i < mole_global.num_total; ++i )
		{
			if( mole.species[i].location != NULL )
				continue;

			molecule::nNucsMap& nuclides = mole_global.list[i]->nNuclide;

			double factor = 1.0;  // running minimum

			for( const auto& entry : nuclides )
				factor = MIN2(factor, abund_ratio[entry.first->el()->Z-1]);

			if( factor != 1.0 )  // abund_ratio should be precisely == 1.0 for an unchanged element 
			{
				mole.species[i].den *= factor;
				for( const auto& atom : nuclides )
					lgModified[atom.first->el()->Z-1] = true;
			}
		}
		
		mole.set_isotope_abundances();
		dense.updateXMolecules();
		deut.updateXMolecules();

		for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem)
		{
			if( !lgModified[nelem] || !dense.lgElmtOn[nelem] )
				continue;
			
			double ion_new = dense.gas_phase[nelem] - dense.xMolecules(nelem);
			double ion_old = 0.0;

			for( long ion=0; ion < (nelem + 2); ++ion )
			{
				ion_old += dense.xIonDense[nelem][ion];  // callers mustn't modify xIonDense
			}

			ASSERT( ion_new >= 0.0 );
			ASSERT( ion_old > 0.0 );  // avoid divide-by-zero; ion_new being 0 is physically possible

			ScaleIonDensities(nelem, (realnum)(ion_new / ion_old));  // minor dtype loss of precision here
		}
	}
}

bool AbundChange( )
{
	DEBUG_ENTRY( "AbundChange()" );

	/* The routine here will probably lead to material being dumped into the ionic species,
	 * but this should slosh back at the next molecular network solve. The current procedure is 
	 * at least a) conservative; b) consistent with the pure-ion case. */	

	/* this will be set true if density or abundances change in this zone */
	bool lgChange = false;

	/* >>chng 97 jun 03, added variable abundances for element table command
	 * option to get abundances off a table with element read command */
	/* >>chng 26 nov 22, updated scaling according to the proposed improved procedure
	 * to avoid lgElemsConserved issues with the element table command & improve accuracy */
	if( abund.lgAbTaON )
	{
		lgChange = true;

		double abund_ratio[LIMELM];
		// initialise ratios to 1.0 then skip hydrogen in the loop so it stays 1
		std::fill_n(abund_ratio, LIMELM, 1.0);

		for( long nelem=ipHELIUM; nelem < LIMELM; ++nelem)
		{
			// skip turned-off elements to avoid NaN/inf denominator
			if( abund.AbunTab[nelem].nvals() > 0 && dense.lgElmtOn[nelem])
			{
				double abun = abund.AbunTab[nelem].tabval(radius.Radius,radius.depth)*
				  dense.gas_phase[ipHYDROGEN];

				// store abundance ratio before setting gas phase density
				abund_ratio[nelem] = abun / dense.gas_phase[nelem];

				// rounding ends up in the ions (see above)
				dense.SetGasPhaseDensity( nelem, (realnum)abun );
			}
		}

		RebalanceAbundances( abund_ratio );
	}
	/* this is set false if fluctuations abundances command entered,
	 * when density variations are in place this is true,
	 * and dense.chDenseLaw is "SINE" */
	double FacAbun = 1.;
	if( !dense.lgDenFlucOn )
	{
		static double FacAbunSav;
		double OldAbun = 0.0;
		/* abundance variations are in place */
		lgChange = true;
		if( nzone > 1 )
		{
			OldAbun = FacAbunSav;
		}

		/* rapid abundances fluctuation */
		if( dense.lgDenFlucRadius )
		{
			/* cycle over radius */
			FacAbunSav = dense.cfirst*cos(radius.depth*dense.flong+
				dense.flcPhase) + dense.csecnd;
		}
		else
		{
			/* cycle over column density */
			FacAbunSav = dense.cfirst*cos( colden.colden[ipCOL_HTOT]*dense.flong+
				dense.flcPhase) + dense.csecnd;
		}

		if( nzone > 1 )
		{
			FacAbun = FacAbunSav/OldAbun;
		}
	}

	if( FacAbun != 1. )  // >>chng 26 nov 22, update this code path to use the new helper
	{
		ASSERT(!dynamics.lgAdvection); // Local abundance change doesn't make sense with advective transport

		double abund_ratio[LIMELM];  
		std::fill_n(abund_ratio, LIMELM, 1.0);  // initialise ratios to 1.0, hydrogen + helium stay 1.0

		/* Li on up is affect by abundance variations */
		for( long nelem=ipLITHIUM; nelem < LIMELM; nelem++ )
		{
			if (dense.lgElmtOn[nelem])
			{
				abund_ratio[nelem] = FacAbun;
				dense.SetGasPhaseDensity(nelem, dense.gas_phase[nelem] * (realnum) FacAbun);
			}
		}
		RebalanceAbundances( abund_ratio );
	}

	if (lgChange)
	{
		/* must call TempChange since ionization has changed, there are some
		 * terms that affect collision rates (H0 term in electron collision) */
		TempChange(phycon.te , false);
	}

	return lgChange;
}

namespace {
	const int SCALE_NEW = 1;
}

realnum scalingDensity(void)
{
	if (SCALE_NEW)
		return dense.xMassDensity/(realnum)ATOMIC_MASS_UNIT;
	else
		return dense.gas_phase[ipHYDROGEN];
}
realnum scalingZoneDensity(long i)
{
	if (SCALE_NEW)
		return struc.DenMass[i]/(realnum)ATOMIC_MASS_UNIT;
	else
		return struc.hden[i];
}
