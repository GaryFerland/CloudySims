/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef SPECIES_H_
#define SPECIES_H_

/*The species structure is used to hold information about a particular atom,ion or molecule
mentioned in the species.ini file.The name of the atom/ion/molecule is used to obtain the density
of molecules in the case of the Leiden Database and along with atomic number and ion stage the 
density of atoms/ions in the case of the CHIANTI database */
class species
{
public:
	species()
	{
		chLabel = NULL;
		database = "";
		index = INT_MAX;
		numLevels_max = -INT_MAX;
		setLevels = -1;
		numLevels_local = -INT_MAX;
		numLevels_masterlist = -1;
		fmolweight = realnum(0.);
		lgMolecular = false;
		fracType = realnum(0.);
		fracIsotopologue = realnum(0.);
		CoolTotal = 0.;
		lgActive = false;
		maxWN = 0.;
		lgLTE = false;
		lgPrtMatrix = false;
	}
	~species()
	{
		delete[] chLabel;
	}
	/*Name of the atom/ion/ molecule*/
	char *chLabel;
	/* Database origin */
	string database;
	// index in chemistry
	long index;
	/*Actual Number of energy levels in the data file*/
	long numLevels_max;
	/*Number of levels defined in "species" command */
	long setLevels;
	/*Number of energy levels used locally*/
	long numLevels_local;
    /** Masterlist defined minimum number of levels, needed for
     *  species like NI which require significant number of
     *  levels to get pumped contributions to important lines */
	long numLevels_masterlist;
	/*Molecular weight*/
	realnum fmolweight;
	/* is molecular? */
	bool lgMolecular;
	/* fraction in this "type" (e.g. para, ortho) */
	realnum fracType;
	/* chemical fractionation */
	realnum fracIsotopologue;
	/** total cooling due to this species */
	double CoolTotal;
	/** are populations currently defined? */
	bool lgActive;
	/* maximum wavenumber in chianti */
	double maxWN;
	/** should level pops be done in LTE? */
	bool lgLTE;
	/** print Matrix input to solver */
	bool lgPrtMatrix;
};

struct t_pseudo_cont
{
	realnum wlLo;
	realnum wlHi;
	long nBins;
};

extern t_pseudo_cont pseudoContDef;

void parsespect(char* chLabel, long& nelem, long& IonStg);
string makeChemical(long nelem, long ion);
void makeChemical(char* chLabelChemical, long nelem, long ion);
void spectral_to_chemical( char *chLabelChemical, char* chLabel );
void chemical_to_spectral( const string chLabelChem, string &chLabelSpec );

class genericState;

/** getSpecies -- acquire the species matching the input string
 *
 * \param speciesLabel	input species string
 * \param species	output reference to requested species
 */
void getSpecies( const string &speciesLabel, genericState &species );

#endif /* SPECIES_ */
