/**
 * @file atmdat_chianti.cpp
 * @brief Implements routines for reading and processing atomic and molecular data from the STOUT and CHIANTI databases for use in Cloudy.
 *
 * This file contains functions and utilities to parse, validate, and map atomic/molecular energy levels, transition probabilities, and collisional data
 * from external data files (STOUT and CHIANTI formats) into Cloudy's internal data structures. It handles various file formats, data consistency checks,
 * and supports special cases for certain species (e.g., iron). The code is responsible for:
 *   - Reading and sorting energy levels.
 *   - Mapping file indices to internal indices, including handling irregular or non-sequential indices.
 *   - Reading and assigning transition probabilities (A-values, gf-values, line strengths).
 *   - Reading and storing collisional data for multiple colliders.
 *   - Error checking and reporting for malformed or incomplete data files.
 *   - Providing debug output for data verification.
 *
 * Key functions:
 *   - getCode(): Converts a transition type string (e.g., "E1", "M2") to an integer code.
 *   - processIndices(): Maps file-based energy level indices to internal indices, handling both regular and irregular cases.
 *   - atmdat_STOUT_readin(): Reads and processes STOUT data files for a given species.
 *   - atmdat_CHIANTI_readin(): Reads and processes CHIANTI data files for a given species.
 *
 * The file also defines supporting structures and constants, such as LevelInfo and ENERGY_MIN_WN.
 *
 * @author Gary J. Ferland and others
 * @copyright Copyright (C)1978-2025 by Gary J. Ferland and others. For conditions of distribution and use see license.txt.
 */
/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "taulines.h"
#include "trace.h"
#include "thirdparty.h"
#include "atmdat.h"
#include "lines_service.h"
#include "parser.h"
#include "parse_species.h"
#include "rfield.h"
#include <cstdlib>  // for strtol
#include <cerrno>   // for errno

typedef vector< pair<double,long> > DoubleLongPairVector;

struct LevelInfo
{
	double nrg;
	long index;
	double stwt;
	string config;
	bool lgTheo;
	LevelInfo(double e, long i, double w, string c, bool t) : nrg(e), index(i), stwt(w), config(c), lgTheo(t) {}
	bool operator< ( const LevelInfo& l ) const
	{
		if( nrg < l.nrg )
			return true;
		else if( nrg == l.nrg )
			return ( index < l.index );
		else
			return false;
	}
};

// minimum energy difference (wavenumbers) we will accept
const double ENERGY_MIN_WN = 1e-10;

/* convert transition type into integer code */
/**
 * @brief Converts a transition type string to its corresponding code.
 *
 * This function interprets a two-character string representing a transition type
 * (e.g., "E1", "M2", "NS") and returns an integer code corresponding to the type:
 *   - "NS": Returns 0, treating the transition as E1 (default, as per NIST).
 *   - "E1", "E2", "E3": Returns 0, 1, or 2 respectively.
 *   - "M1", "M2", "M3": Returns 3, 4, or 5 respectively.
 * If the input string does not match the expected format, returns -1.
 *
 * @param transType A two-character string representing the transition type.
 * @return int The corresponding code for the transition type, or -1 if invalid.
 */
inline int getCode(const string& transType)
{
	DEBUG_ENTRY( "getCode()" );

	int val = -1;
	if( transType.length() != 2 )
		return -1;
	if( transType == "NS" )
		return 0; // when not set, treat the transition as E1 (this is what NIST does)
	if( transType[0] == 'E' )
		val = 0;
	else if( transType[0] == 'M' )
		val = 3;
	else
		return -1;
	if( transType[1] == '1' )
		(void)0;
	else if( transType[1] == '2' )
		val += 1;
	else if( transType[1] == '3' )
		val += 2;
	else
		return -1;
	return val;
}

/**
 * @brief Processes and maps input energy level indices to internal indices, handling both regular and irregular cases.
 *
 * This function adjusts the lower and upper energy level indices (`ipLo`, `ipHi`) based on whether the indices are regular
 * (sequential and zero-based) or irregular (requiring mapping via `indexold2new`). If the indices are irregular and not found
 * in the mapping, both output indices are set to `LONG_MAX`. The function also ensures that the lower index is less than or
 * equal to the upper index, swapping them if necessary.
 *
 * @param[in]  ipLoInFile      The lower energy level index as read from the file (1-based).
 * @param[in]  ipHiInFile      The upper energy level index as read from the file (1-based).
 * @param[in]  lgIsRegular     Flag indicating if the indices are regular (true) or require mapping (false).
 * @param[in]  indexold2new    Mapping from old (file) indices to new (internal) indices.
 * @param[out] ipLo            The mapped lower energy level index (0-based).
 * @param[out] ipHi            The mapped upper energy level index (0-based).
 */
inline void processIndices(long ipLoInFile, long ipHiInFile, bool lgIsRegular, const map<long,long>& indexold2new,
			   long& ipLo, long& ipHi)
{
	DEBUG_ENTRY( "processIndices()" );

	if( lgIsRegular )
	{
		ipLo = ipLoInFile-1;
		ipHi = ipHiInFile-1;
	}
	else
	{
		// ignore transitions with an index that is not present in the level file
		// this is not an error as the eof marker in the level file may have been moved upward
		auto plo = indexold2new.find(ipLoInFile);
		if( plo == indexold2new.end() )
		{
			ipLo = ipHi = LONG_MAX;
			return;
		}
		auto phi = indexold2new.find(ipHiInFile);
		if( phi == indexold2new.end() )
		{
			ipLo = ipHi = LONG_MAX;
			return;
		}

		/* Account for reordered energy levels */
		ipLo = plo->second;
		ipHi = phi->second;
	}

	ASSERT( ipLo != ipHi );

	// swap indices if energy levels were not correctly sorted
	if( ipHi < ipLo )
	{
		long swap = ipHi;
		ipHi = ipLo;
		ipLo = swap;
	}
}

STATIC void ApplyTPDatum(
	DataParser& d,
	const string& dataType,
	long ipLoInFile,
	long ipHiInFile,
	double tpData,
	const string& transTypeRaw,
	bool lgOverrideExisting,
	long intNS,
	long nMolLevs,
	bool lgIsRegular,
	const map<long,long>& indexFile2Sorted,
	multi_arr<bool,3,C_TYPE>& lgLineStrengthTT,
	multi_arr<bool,2,C_TYPE>& lgResetThisTP,
	const vector<LevelInfo>& dBaseStatesOrg,
	int numTransitionTypes )
{
	DEBUG_ENTRY( "ApplyTPDatum()" );

	long ipLo, ipHi;
	processIndices( ipLoInFile, ipHiInFile, lgIsRegular, indexFile2Sorted, ipLo, ipHi );

	// unmapped indices from irregular level files are ignored
	if( ipLo == LONG_MAX || ipHi == LONG_MAX )
		return;

	if( ipHi >= nMolLevs )
		return;

	if( tpData < 0. )
		d.errorAbort( "transition probability data must be >= 0" );

	if( dataType == "A" && tpData < atmdat.aulThreshold )
		return;

	string transType = transTypeRaw;
	if( transType.empty() )
		transType = "NS";

	if( transType == "UT" )
		transType = "M1+E2";

	if( dataType == "S" )
	{
		if( transType == "NS" )
			d.errorAbort( "line strength specified, but no transition type was found" );
		if( transType.length() != 2 )
			d.errorAbort( "invalid transition type" );
	}

	TransitionList::iterator tr =
		dBaseTrans[intNS].begin() + ipdBaseTrans[intNS][ipHi][ipLo];

	if( lgOverrideExisting && tr->hasEmis() && !lgResetThisTP[ipLo][ipHi] )
	{
		tr->Emis().Aul() = 0.;
		tr->Emis().gf() = 0.;

		for( int k=0; k < numTransitionTypes; ++k )
			lgLineStrengthTT[ipLo][ipHi][k] = false;

		lgResetThisTP[ipLo][ipHi] = true;
	}

	if( !tr->hasEmis() )
	{
		tr->AddLine2Stack();
		tr->Emis().Aul() = 0.;
		tr->Emis().gf() = 0.;
	}

	size_t len = transType.length();
	size_t p = 0;
	while( p < len )
	{
		string tt = transType.substr( p, 2 );
		int ttCode = getCode( tt );
		char c = ( p+2 >= len ) ? '+' : transType[p+2];

		if( ttCode < 0 || c != '+' )
			d.errorAbort( "invalid transition type" );

		if( lgLineStrengthTT[ipLo][ipHi][ttCode] )
			d.errorAbort( "this transition already has an Aul value set" );

		lgLineStrengthTT[ipLo][ipHi][ttCode] = true;
		p += 3;
	}

	if( tr->EnergyWN() > ENERGY_MIN_WN )
	{
		if( dataType == "A" )
		{
			tr->Emis().Aul() += tpData;
			tr->Emis().gf() =
				(realnum)GetGF( tr->Emis().Aul(), tr->EnergyWN(), tr->Hi()->g() );
		}
		else if( dataType == "G" )
		{
			tr->Emis().gf() += tpData;
			tr->Emis().Aul() =
				(realnum)eina( tr->Emis().gf(), tr->EnergyWN(), tr->Hi()->g() );
		}
		else
		{
			tr->Emis().Aul() +=
				S2Aul( tpData, tr->WLangVac(), tr->Hi()->g(), transType );
			tr->Emis().gf() =
				(realnum)GetGF( tr->Emis().Aul(), tr->EnergyWN(), tr->Hi()->g() );
		}
	}

	tr->setComment(
		db_comment_tran_levels(
			dBaseStatesOrg[ipLo].config,
			dBaseStatesOrg[ipHi].config ) );
}
/**
 * @brief Extract numeric block index from a STOUT block filename.
 *
 * This function searches the filename for the substring "_blk_"
 * and extracts the integer that immediately follows it.
 *
 * Example:
 *     "al_8_blk_20.tp"  →  returns 20
 *
 * The extracted number is used to:
 *     - Sort block files numerically
 *     - Determine baseline block (lowest N)
 *     - Control override priority (higher N overrides lower N)
 *
 * If the filename does not contain "_blk_" or no digits follow it,
 * the program aborts with an error.
 */
static long ExtractBlkIndexOrDie( const string& fn, const string& shortName, const string& extension )
{
	DEBUG_ENTRY( "ExtractBlkIndexOrDie()" );

	// v2 block files must use exactly:
	//   <shortName>_blk_<N><extension>
	// e.g.
	//   al_8_blk_10.tp
	//   al_8_blk_10.coll
	const string prefix = shortName + "_blk_";

	if( fn.compare( 0, prefix.size(), prefix ) != 0 )
	{
		fprintf( ioQQQ,
			" PROBLEM: Invalid STOUT v2 block filename '%s'.\n"
			" Expected filename to start with '%s'.\n",
			fn.c_str(), prefix.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	if( fn.size() <= prefix.size() + extension.size() ||
	    fn.compare( fn.size() - extension.size(), extension.size(), extension ) != 0 )
	{
		fprintf( ioQQQ,
			" PROBLEM: Invalid STOUT v2 block filename '%s'.\n"
			" Expected filename format is '%s<N>%s'.\n",
			fn.c_str(), prefix.c_str(), extension.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	const string blockString =
		fn.substr( prefix.size(), fn.size() - prefix.size() - extension.size() );

	if( blockString.empty() ||
	    blockString.find_first_not_of( "0123456789" ) != string::npos )
	{
		fprintf( ioQQQ,
			" PROBLEM: Invalid STOUT v2 block filename '%s'.\n"
			" The block index must be a non-negative integer in '%s<N>%s'.\n",
			fn.c_str(), prefix.c_str(), extension.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	char* endptr = nullptr;
	errno = 0;
	long blockIndex = strtol( blockString.c_str(), &endptr, 10 );

	if( errno != 0 || endptr == blockString.c_str() || *endptr != '\0' )
	{
		fprintf( ioQQQ,
			" PROBLEM: Failed to read block index from STOUT v2 filename '%s'.\n",
			fn.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	return blockIndex;
}
/**
 * @brief Reads and processes STOUT atomic/molecular data files for a given species.
 *
 * This function reads the STOUT data files (.nrg, .tp, .coll) for a specified species,
 * parses the energy levels, transition probabilities, and collisional data, and
 * initializes the corresponding data structures used by the code. It handles
 * various file formats, checks for data consistency, and supports special cases
 * (e.g., Fe species with more levels). The function also performs error checking
 * and outputs debug information if enabled.
 *
 * @param intNS
 *    The index of the species in the dBaseSpecies array.
 * @param chPrefix
 *    The prefix string used to construct the filenames for the STOUT data files.
 *
 * @details
 * The function performs the following steps:
 *   - Reads and validates the energy levels file (.nrg), sorts levels, and initializes state arrays.
 *   - Reads the transition probability file (.tp), processes radiative data, and populates transition arrays.
 *   - Reads the collision data file (.coll), processes collisional strengths/rates, and fills collisional data arrays.
 *   - Handles special cases for Fe species and user-specified level limits.
 *   - Performs extensive error checking and outputs debug information if DEBUGSTATE is enabled.
 *
 * @throws
 *   Exits the program with an error message if any file is malformed, missing required data,
 *   or contains invalid entries.
 */
void atmdat_STOUT_readin( long intNS, const string& chPrefix )
{
	DEBUG_ENTRY( "atmdat_STOUT_readin()" );

	const int MAX_NUM_LEVELS = 999;
	const bool DEBUGSTATE = false;

	// versions one and two magic numbers for the stout files
	const long iyr_v1 = 17, imo_v1 = 9, idy_v1 = 5;
	const long iyr_v2 = 25, imo_v2 = 10, idy_v2 = 15;

	dBaseSpecies[intNS].lgMolecular = false;
	dBaseSpecies[intNS].lgLTE = false;

	setProperties(dBaseSpecies[intNS]);

	string chUnCaps = chPrefix;
	uncaps(chUnCaps);

	if( dBaseSpecies[intNS].dataset.length() > 0 )
		chUnCaps += "_" + dBaseSpecies[intNS].dataset;
	// Construct the filename for the STOUT energy levels file (.nrg) based on the species prefix.
	string chNRGFilename = chUnCaps + ".nrg";
	string chTPFilename  = chUnCaps + ".tp";
	
	// Resolve STOUT .coll filenames from the required directory structure.
	// chUnCaps is something like: "stout/al/al_8/al_8"
	string basedir;
	string shortName;

	string::size_type pos = chUnCaps.find_last_of('/');
	if( pos == string::npos )
	{
		fprintf( ioQQQ,
				" PROBLEM: STOUT data files must follow the required directory structure.\n"
				" The supplied STOUT prefix '%s' does not contain a directory component.\n",
				chUnCaps.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	basedir  = chUnCaps.substr( 0, pos+1 );      // e.g. "stout/al/al_8/"
	shortName = chUnCaps.substr( pos+1 );        // e.g. "al_8"

	vector<string> collFiles;
	string collPattern = basedir + shortName + ".*\\.coll";
	getFileList( collFiles, collPattern );
	if( collFiles.empty() )
	{
		fprintf( ioQQQ,
				" PROBLEM: No STOUT .coll file found matching pattern '%s'\n",
				collPattern.c_str() );
		cdEXIT( EXIT_FAILURE );
	}
	// Split coll files by magic number.
	// V1 magic files are in single-file format and do NOT require _blk_.
	// V2 magic files must be block files and will be ordered by blk index.
	vector<string> collFilesV1;
	vector<string> collFilesV2;

	for( const auto& f : collFiles )
	{
		DataParser dp;
		dp.open( basedir + f, ES_STARS_ONLY );

		long yr, mo, dy;
		dp.getline();
		dp.getToken( yr );
		dp.getToken( mo );
		dp.getToken( dy );

		if( yr == iyr_v1 && mo == imo_v1 && dy == idy_v1 )
		{
			collFilesV1.push_back( f );
		}
		else if( yr == iyr_v2 && mo == imo_v2 && dy == idy_v2 )
		{
			collFilesV2.push_back( f );
		}
		else
		{
			fprintf(ioQQQ, " PROBLEM: Invalid magic number in STOUT .coll file %s\n", (basedir+f).c_str() );
			cdEXIT(EXIT_FAILURE);
		}
	}

	DataParser d;


	/******************************************************
	 ***************** Energy Levels File *****************
	 ******************************************************/
	d.open( chNRGFilename, ES_STARS_ONLY );

	long index = 0;
	double nrg = 0.0;
	double istwt, stwt = 0.0;

	/* first line is a version number - now confirm that it is valid */
	d.getline();
	long myr, mmo, mdy;
	d.getToken(myr);
	d.getToken(mmo);
	d.getToken(mdy);

	bool lgNRGv1 = ( myr == iyr_v1 && mmo == imo_v1 && mdy == idy_v1 );
	bool lgNRGv2 = ( myr == iyr_v2 && mmo == imo_v2 && mdy == idy_v2 );

	if( !lgNRGv1 && !lgNRGv2 )
	{
		d.errorAbort("invalid magic number in STOUT .nrg file");
	}

	/* Create array for holding energies and statistical weights so that
	 * we can put the energies in the correct order before moving them to
	 * dBaseStates */
	vector<LevelInfo> dBaseStatesOrg;

	// check for an end of file sentinel
	bool lgSentinelReached = false;
	while( d.getline() )
	{
		if( d.lgEODMarker() )
		{
			lgSentinelReached = true;
			break;
		}

		d.getToken( index );
		bool lgTheo = d.getLvlEnergy( nrg );
		d.getToken( stwt );
		if( stwt <= 0. || modf(stwt, &istwt) != 0. )
			d.errorAbort( "invalid statistical weight" );
		string config;
		d.getQuoteOptional( config );

		d.checkEOL();

		dBaseStatesOrg.emplace_back(nrg, index, stwt, config, lgTheo);
	}

	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n", chNRGFilename.c_str() );
		fprintf( ioQQQ, " Check that stars (***) appear after the last line of data"
			 " and start in the first column of that line.\n");
		cdEXIT(EXIT_FAILURE);
	}

	long nMolLevs = dBaseStatesOrg.size();
	long HighestIndexInFile = nMolLevs;

	dBaseSpecies[intNS].numLevels_max = HighestIndexInFile;

	if( tolower(dBaseSpecies[intNS].chLabel[0]) == 'f' && tolower(dBaseSpecies[intNS].chLabel[1]) == 'e' )
	{
		// Fe is special case with more levels
		nMolLevs = MIN3(nMolLevs, atmdat.nStoutMaxLevelsFe, MAX_NUM_LEVELS );
	}
	else
	{
		nMolLevs = MIN3(nMolLevs, atmdat.nStoutMaxLevels, MAX_NUM_LEVELS );
	}

	//Consider the masterlist specified number of levels as the min. =1 if not specified
	long numMasterlist = MIN2( dBaseSpecies[intNS].numLevels_masterlist, HighestIndexInFile );
	nMolLevs = MAX2(nMolLevs,numMasterlist);

	if (dBaseSpecies[intNS].setLevels != -1)
	{
		if (dBaseSpecies[intNS].setLevels > HighestIndexInFile)
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel );
			if( dBaseSpecies[intNS].setLevels != LONG_MAX )
				fprintf( ioQQQ, "Using STOUT spectrum %s (species: %s) with %li requested,"
					 " only %li energy levels available.\n",
					 dBaseSpecies[intNS].chLabel, chLabelChemical, dBaseSpecies[intNS].setLevels,
					 HighestIndexInFile );
			nMolLevs = HighestIndexInFile;	  
		}
		else
		{
			nMolLevs = dBaseSpecies[intNS].setLevels;
		}
	}

	ASSERT( nMolLevs <= HighestIndexInFile );

	if( atmdat.lgStoutPrint )
	{
		char chLabelChemical[CHARS_SPECIES] = "";
		spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
		fprintf( ioQQQ, "Using STOUT spectrum %s (species: %s) with %li levels of %li available.\n",
				 dBaseSpecies[intNS].chLabel, chLabelChemical, nMolLevs, HighestIndexInFile );
	}

	dBaseSpecies[intNS].numLevels_max = nMolLevs;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	/*Resize the States array*/
	dBaseStates[intNS].init(dBaseSpecies[intNS].chLabel,nMolLevs);
	/*Zero out the maximum wavenumber for each species */
	dBaseSpecies[intNS].maxWN = 0.;

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nStout Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Energy Level File: %s\n",chNRGFilename.c_str());
		fprintf(ioQQQ,"Number of Energy Levels in File: %li\n",HighestIndexInFile);
		fprintf(ioQQQ,"Number of Energy Levels Cloudy is Currently Using: %li\n",nMolLevs);
		fprintf(ioQQQ,"Species|File Index|Energy(WN)|StatWT\n");
		for( size_t i=0; i < dBaseStatesOrg.size(); ++i )
		{
			fprintf( ioQQQ, "<%s>\t%li\t%.3f\t%.1f\n", dBaseSpecies[intNS].chLabel,
				 dBaseStatesOrg[i].index, dBaseStatesOrg[i].nrg, dBaseStatesOrg[i].stwt );
		}
	}

	// STOUT v1 files may be unsorted, so cloudy sorts them here.
	// STOUT v2 requires the energy level file to already be sorted.
	if( lgNRGv1 )
	{
		sort(dBaseStatesOrg.begin(), dBaseStatesOrg.end());
	}
	else
	{
		for( size_t i = 1; i < dBaseStatesOrg.size(); ++i )
		{
			if( dBaseStatesOrg[i] < dBaseStatesOrg[i-1] )
			{
				fprintf( ioQQQ,
					" PROBLEM: STOUT v2 .nrg file is not sorted by increasing energy.\n"
					" Previous level: index %ld, energy %.10e\n"
					" Current  level: index %ld, energy %.10e\n",
					dBaseStatesOrg[i-1].index,
					dBaseStatesOrg[i-1].nrg,
					dBaseStatesOrg[i].index,
					dBaseStatesOrg[i].nrg );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	// highest level index in file-space that is still retained after truncation
	//Why here: This value is needed by both the v2 .tp and v2 .coll readers.
	long maxRetainedIndexInFile = -1;
	for( long i = 0; i < nMolLevs; ++i )
		maxRetainedIndexInFile = MAX2( maxRetainedIndexInFile, dBaseStatesOrg[i].index );

	// regular files have the level indices running as 1, 2, 3, ...
	// (without gaps) and the energies are also already sorted
	bool lgIsRegular = true;
	for( size_t i=0; i < dBaseStatesOrg.size(); i++ )
	{
		if( dBaseStatesOrg[i].index != long(i+1) )
		{
			lgIsRegular = false;
			break;
		}
	}

	map<long, long> indexold2new;
	if( !lgIsRegular ) {
		for( size_t i=0; i < dBaseStatesOrg.size(); i++ )
		{
			long old = dBaseStatesOrg[i].index;
			auto p = indexold2new.find(old);
			if( p == indexold2new.end() )
				indexold2new[old] = i;
			else
			{
				fprintf(ioQQQ, "Duplicate index %ld found in file %s\n", old, chNRGFilename.c_str());
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\n\n***Energy levels have been sorted in order of increasing energy.***\n");
		fprintf(ioQQQ,"Species|File Index|Sorted Index|Energy(WN)|StatWT\n");
	}

	/* Store sorted energies in their permanent home */
	for( long i=0; i < nMolLevs; i++ )
	{
		double nrg = dBaseStatesOrg[i].nrg;
		long oldindex = dBaseStatesOrg[i].index;
		double stwt = dBaseStatesOrg[i].stwt;
		bool lgTheory = dBaseStatesOrg[i].lgTheo;

		if( DEBUGSTATE )
			fprintf( ioQQQ, "<%s>\t%li\t%li\t%.3f\t%.1f\n", dBaseSpecies[intNS].chLabel,
				 oldindex, i+1, nrg, stwt );

		dBaseStates[intNS][i].energy().set(nrg,"cm^-1");
		dBaseStates[intNS][i].g() = stwt;
		dBaseStates[intNS][i].ipOrg() = oldindex;
		dBaseStates[intNS][i].theory() = lgTheory;
	}

	/* allocate the Transition array*/
	ipdBaseTrans[intNS].reserve(nMolLevs);
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
		ipdBaseTrans[intNS].reserve(ipHi,ipHi);
	ipdBaseTrans[intNS].alloc();
	dBaseTrans[intNS].resize(ipdBaseTrans[intNS].size());
	dBaseTrans[intNS].states() = &dBaseStates[intNS];
	dBaseTrans[intNS].chLabel() = dBaseSpecies[intNS].chLabel;
	dBaseSpecies[intNS].database = "Stout";

	//This is creating transitions that we don't have collision data for
	//Maybe use gbar or keep all of the Fe 2 even if it was assumed (1e-5)
	int ndBase = 0;
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
	{
		for( long ipLo = 0; ipLo < ipHi; ipLo++)
		{
			ipdBaseTrans[intNS][ipHi][ipLo] = ndBase;
			dBaseTrans[intNS][ndBase].Junk();
			dBaseTrans[intNS][ndBase].setLo(ipLo);
			dBaseTrans[intNS][ndBase].setHi(ipHi);
			++ndBase;
		}
	}

	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for( auto tr=dBaseTrans[intNS].begin(); tr != dBaseTrans[intNS].end(); ++tr )
	{
		int ipHi = tr->ipHi();
		int ipLo = tr->ipLo();
		ASSERT( ipHi > ipLo &&
			dBaseStates[intNS][ipHi].energy().WN() >= dBaseStates[intNS][ipLo].energy().WN() );
		double fenergyWN = dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN();
		tr->EnergyWN() = fenergyWN;
		if( rfield.isEnergyBound( Energy( fenergyWN, "cm^-1" ) ) )
		{
			tr->WLangVac() = wn2angVac( fenergyWN );
			dBaseSpecies[intNS].maxWN = MAX2(dBaseSpecies[intNS].maxWN,tr->EnergyWN());
		}
		else
			tr->WLangVac() = 1e30;
	}

/******************************************************
************* Transition Probability File ************
******************************************************/

/*
*V1 vs V2 STOUT .tp handling:
*
*   - V1 STOUT (magic 17 09 05):
*       * File is version-1-format: "<shortName>.tp" (no _blk_)
*       * Each data row starts with dataType:  A|G|S
*       * End-of-data MUST be the "***" sentinel (exactly as before)
*
*   - V2 STOUT (magic 25 10 15):
*       * File(s) may be split as "<shortName>_blk_<N>.tp"
*       * After magic line there is a global dataType line: A|G|S
*       * Each data row is: ipLo ipHi value [transType]
*       * End-of-data can be "***" OR an empty line
*       * Multiple blocks are applied in blk order; later blocks OVERRIDE earlier ones
*/

// ---- discover tp candidates (version 1 and/or version 2) ----
	vector<string> tpFiles;
	string tpPattern = basedir + shortName + ".*\\.tp";
	getFileList( tpFiles, tpPattern );

	vector<string> tpFilesv2; // contains _blk_
	vector<string> tpFilesv1; // version 1 (no _blk_)
	for( const auto& fn : tpFiles )
	{
		if( fn.find("_blk_") != string::npos )
			tpFilesv2.push_back(fn);
		else
			tpFilesv1.push_back(fn);
	}

	struct TPFileInfo
	{
		string filename;
		long blockIndex;
		bool lgBlockFile;
	};

	vector<TPFileInfo> tpFilesToProcess;
	map<long,string> tpBlockFiles;
	long baselineTpBlk = LONG_MAX;

	if( !tpFilesv2.empty() )
	{
		// Version 2 format uses block files.
		for( const string& tpFile : tpFilesv2 )
		{
			long blockIndex = ExtractBlkIndexOrDie( tpFile, shortName, ".tp" );

			if( tpBlockFiles.find( blockIndex ) != tpBlockFiles.end() )
			{
				fprintf( ioQQQ,
					" PROBLEM: Duplicate STOUT v2 .tp block index %ld found:\n"
					"  %s\n"
					"  %s\n",
					blockIndex,
					tpBlockFiles[blockIndex].c_str(),
					tpFile.c_str() );
				cdEXIT( EXIT_FAILURE );
			}

			tpBlockFiles[blockIndex] = tpFile;
		}

		baselineTpBlk = tpBlockFiles.begin()->first;

		for( const auto& tpBlock : tpBlockFiles )
		{
			TPFileInfo tpInfo;
			tpInfo.filename = tpBlock.second;
			tpInfo.blockIndex = tpBlock.first;
			tpInfo.lgBlockFile = true;
			tpFilesToProcess.push_back( tpInfo );
		}
	}
	else
	{
		// Version 1 format uses a single non-block file.
		if( tpFilesv1.empty() )
		{
			fprintf( ioQQQ,
				" PROBLEM: No STOUT .tp file found matching pattern '%s'\n",
				tpPattern.c_str() );
			cdEXIT( EXIT_FAILURE );
		}

		string tpFile;

		if( tpFilesv1.size() == 1 )
		{
			tpFile = tpFilesv1.front();
		}
		else
		{
			const string preferred = shortName + ".tp";
			auto it = find( tpFilesv1.begin(), tpFilesv1.end(), preferred );

			if( it == tpFilesv1.end() )
			{
				fprintf( ioQQQ,
					" PROBLEM: Multiple version 1 STOUT .tp files found for '%s' but none matches '%s'\n",
					shortName.c_str(), preferred.c_str() );
				cdEXIT( EXIT_FAILURE );
			}

			tpFile = *it;
		}

		TPFileInfo tpInfo;
		tpInfo.filename = tpFile;
		tpInfo.blockIndex = LONG_MAX;
		tpInfo.lgBlockFile = false;
		tpFilesToProcess.push_back( tpInfo );
	}

	// ---- line-strength bookkeeping (shared across all tp blocks/files) ----
	static const int intNumCols = 6;
	multi_arr<bool,3,C_TYPE> lgLineStrengthTT(nMolLevs, nMolLevs, intNumCols);
	lgLineStrengthTT = false;

	// ---- track which transitions have been "reset" for overrides (per file/block pass) ----
	multi_arr<bool,2,C_TYPE> lgResetThisTP(nMolLevs, nMolLevs);
	lgResetThisTP = false;


	// ---- now actually read each tp file ----
	for( const TPFileInfo& tpInfo : tpFilesToProcess )
	{
		// Reset per-block override bookkeeping so each block can reset each transition once.
		lgResetThisTP = false;

		const string thisTP = basedir + tpInfo.filename;

		// In v2 block mode, blocks after the baseline block override earlier data.
		const bool lgOverrideExisting =
			( tpInfo.lgBlockFile && tpInfo.blockIndex > baselineTpBlk );
	

		d.open( thisTP, ES_STARS_ONLY );

		// read and validate magic
		d.getline();
		d.getToken(myr);
		d.getToken(mmo);
		d.getToken(mdy);

		bool lgTPv1 = (myr == iyr_v1 && mmo == imo_v1 && mdy == idy_v1);
		bool lgTPv2 = (myr == iyr_v2 && mmo == imo_v2 && mdy == idy_v2);

		if( !lgTPv1 && !lgTPv2 )
			d.errorAbort("invalid magic number in STOUT .tp file");

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"\nStout Species: %s\n", dBaseSpecies[intNS].chLabel);
			fprintf(ioQQQ,"Radiative Data File: %s\n", thisTP.c_str());
		}

		// ==========================================================
		// FORMAT Version 1  (magic number 17 09 05) --
		// ==========================================================
		if( lgTPv1 )
		{
			if( DEBUGSTATE )
				fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Data Type (A,G,S)|Data\n");

			lgSentinelReached = false;
			while( d.getline() )
			{
				if( d.lgEODMarker() )
				{
					lgSentinelReached = true;
					break;
				}

				string dataType;
				d.getToken( dataType );

				if( dataType != "A" && dataType != "G" && dataType != "S" )
					d.errorAbort( "invalid data type" );

				long ipLoInFile, ipHiInFile;
				d.getToken( ipLoInFile );
				d.getToken( ipHiInFile );

				double tpData;
				d.getToken( tpData );

				string transType;
				if( !d.getTokenOptional( transType ) )
					transType = "NS";

				d.checkEOL(); // legacy strictness stays

				ApplyTPDatum(
					d,
					dataType,
					ipLoInFile,
					ipHiInFile,
					tpData,
					transType,
					false,
					intNS,
					nMolLevs,
					lgIsRegular,
					indexold2new,
					lgLineStrengthTT,
					lgResetThisTP,
					dBaseStatesOrg,
					intNumCols );
			}

			if( !lgSentinelReached )
			{
				fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n", thisTP.c_str() );
				fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data"
					" and start in the first column of that line.");
				cdEXIT(EXIT_FAILURE);
			}
		}
		// ==========================================================
		// Format Version 2 (magic number 25 10 15): block-based TP format.
		// ==========================================================
		else
		{
			// line after magic contains global data type: A|G|S
			d.getline();
			string globalType;
			d.getToken( globalType );
			if( globalType != "A" && globalType != "G" && globalType != "S" )
				d.errorAbort( "invalid global data type in STOUT .tp file" );
			d.checkEOL();
			
			lgSentinelReached = false;
			bool lgEarlyStop = false;
			long prevHiInFile = -1;
			long prevLoInFile = -1;
		
			while( d.getline() )
			{
				if( d.lgEODMarker() )
				{
					lgSentinelReached = true;
					break;
				}

				// allow blank line as end-of-data in v2 format
				long ipLoInFile;
				if( !d.getTokenOptional(ipLoInFile) )
				{
					lgSentinelReached = true;
					break;
				}
				long ipHiInFile;
				d.getToken( ipHiInFile );
				// verify v2 ordering by (upper, lower) file indices
				if( ipHiInFile < prevHiInFile ||
			    	( ipHiInFile == prevHiInFile && ipLoInFile < prevLoInFile ) )
				{
					d.errorAbort( "STOUT v2 .tp file is not sorted by (upper, lower) level index" );
				}
				prevHiInFile = ipHiInFile;
				prevLoInFile = ipLoInFile;

				// safe early termination: all later rows are beyond retained levels
				if( ipHiInFile > maxRetainedIndexInFile )
				{
					lgEarlyStop = true;
					break;
				}

				double tpData;
				d.getToken( tpData );

				string transType;
				if( !d.getTokenOptional( transType ) )
					transType = "NS";

				d.checkEOL();

				ApplyTPDatum(
					d,
					globalType,
					ipLoInFile,
					ipHiInFile,
					tpData,
					transType,
					lgOverrideExisting,
					intNS,
					nMolLevs,
					lgIsRegular,
					indexold2new,
					lgLineStrengthTT,
					lgResetThisTP,
					dBaseStatesOrg,
					intNumCols );
			}

			if( !lgSentinelReached && !lgEarlyStop )
			{
				fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n", thisTP.c_str() );
				fprintf( ioQQQ, " STOUT v2 .tp expects stars (*****) or an empty line at end-of-data.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	} // end loop over tpFilesToProcess




	/******************************************************
	 ************* Collision Data File ********************
	******************************************************/
	int numpoints = 0;
	vector<double> temps;
	long ipCollider = -1;

	/****** Could add ability to count number of temperature changes, electron CS, and proton CS ****/
	/* Allocate space for collision strengths ONCE per species (before reading any blocks/files) */
	StoutCollData[intNS].alloc(nMolLevs, nMolLevs, ipNCOLLIDER);
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			for( long k=0; k<ipNCOLLIDER; k++ )
			{
				/* initialize all spline variables */
				StoutCollData[intNS].junk(ipHi, ipLo, k);
			}
		}
	}

	/*
	* Decide which collision set to process:
	*   - V2-magic (25 10 15): supports block files named *_blk_<N>.coll
	*   - V1-magic (17 09 05): single-file format like al_1.coll (NO _blk_)
	*
	* IMPORTANT:
	*   We must NOT require _blk_ for version-1-format files.
	*   Numeric blk ordering and override behavior only apply to block files.
	*/
	vector<string> collToProcess;
	bool lgUsingV2CollBlocks = false;

	if( !collFilesV2.empty() )
	{
		/* We have at least one *_blk_<N>.coll candidate -> treat as v2 block mode */
		collToProcess = collFilesV2; /* only blk files have v2 magic, so this is safe */
		lgUsingV2CollBlocks = true;
	}
	else
	{
		/* No block files -> version-1-format mode */
		collToProcess = collFilesV1;
		lgUsingV2CollBlocks = false;
	}

	/*
	* In v1 mode we should not “merge” multiple files silently.
	* If multiple v1 candidates exist, it likely indicates a naming problem.
	* We can relax this if you truly want to read multiple V1-format files.
	*/
	if( !lgUsingV2CollBlocks && collToProcess.size() > 1 )
	{
	/* Prefer the exact "<shortName>.coll" if present, otherwise abort */
		string preferred = shortName + ".coll";
		auto it = find( collToProcess.begin(), collToProcess.end(), preferred );
		if( it != collToProcess.end() )
		{
			collToProcess.clear();
			collToProcess.push_back( preferred );
		}
		else
		{
			fprintf( ioQQQ,
				" PROBLEM: Multiple version-1-format STOUT .coll files match pattern '%s' for species %s\n",
				(basedir + shortName + ".*\\.coll").c_str(),
				dBaseSpecies[intNS].chLabel );
			fprintf( ioQQQ,
				" None is a unique default (expected exactly one version-1-format file like '%s').\n",
				preferred.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/*
	* Build processing order:
	*   - V2 block mode: sort by numeric blk index (ascending), baseline = lowest blk
	*   - V1 mode: just process the single file (or lexicographic order if you later allow >1)
	*/
	vector<size_t> collOrder(collToProcess.size());
	for( size_t i=0; i<collOrder.size(); ++i )
		collOrder[i] = i;

	vector<long> collBlk; /* only used in block mode */
	long baselineBlk = LONG_MIN;

	if( lgUsingV2CollBlocks )
	{
		collBlk.assign( collToProcess.size(), LONG_MAX );
		for( size_t i=0; i<collToProcess.size(); ++i )
		{
			collBlk[i] = ExtractBlkIndexOrDie( collToProcess[i], shortName, ".coll" ); /* requires _blk_ */
		}

		sort( collOrder.begin(), collOrder.end(),
			[&](size_t a, size_t b)
			{
				return collBlk[a] < collBlk[b];
			});

		baselineBlk = collBlk[ collOrder.front() ];
	}
	else
	{
	/* version-1-format mode: lexicographic is fine (usually only 1 file anyway) */
		sort( collOrder.begin(), collOrder.end(),
			[&](size_t a, size_t b)
			{
				return collToProcess[a] < collToProcess[b];
			});
	}

	/* Loop over selected .coll files */
	for( size_t ii=0; ii < collOrder.size(); ++ii )
	{
		size_t iFile = collOrder[ii];

		string chCOLLFilename = basedir + collToProcess[iFile];

	/*
	* Override semantics:
	*   - Only in the block mode (Version 2 magic number): higher blk overrides lower blk
	*   - In version 1 mode: never override (duplicates are errors as before)
	*/
		bool lgOverrideExisting = false;
		if( lgUsingV2CollBlocks )
			lgOverrideExisting = ( collBlk[iFile] > baselineBlk );

		d.open( chCOLLFilename, ES_STARS_ONLY );

		/* first line is a version number - now confirm that it is valid */
		d.getline();
		d.getToken(myr);
		d.getToken(mmo);
		d.getToken(mdy);

		bool lgCOLLv1 = (myr == iyr_v1 && mmo == imo_v1 && mdy == idy_v1);
		bool lgCOLLv2 = (myr == iyr_v2 && mmo == imo_v2 && mdy == idy_v2);

		if( !lgCOLLv1 && !lgCOLLv2 )
			d.errorAbort("Invalid magic number in STOUT .coll file");

	/* Reset per-file state */
		numpoints = 0;
		temps.clear();
		ipCollider = -1;

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"\nStout Species: %s\n", dBaseSpecies[intNS].chLabel);
			fprintf(ioQQQ,"Collision Data File: %s\n", chCOLLFilename.c_str());
			fprintf(ioQQQ,"Species|TEMP|Temperatures (K)\n");
			fprintf(ioQQQ,"Species|Data Type (CS,RATE)|Collider|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Data\n");
		}

	// ---------------- Version 1 COLLISION FORMAT (17 09 05) ----------------
		if( lgCOLLv1 )
		{
			lgSentinelReached = false;
			while( d.getline() )
			{
	/* Stop on *** */
				if( d.lgEODMarker() )
				{
					lgSentinelReached = true;
					break;
				}

				string dataType;
				d.getToken( dataType );

	// Temperature line
				if( dataType == "TEMP" )
				{
					if( DEBUGSTATE )
						fprintf(ioQQQ,"<%s>\tTEMP\t", dBaseSpecies[intNS].chLabel);

					temps.clear();
					double data;
					while( d.getTokenOptional(data) )
					{
						if( data <= 0. )
							d.errorAbort( "invalid temperature" );
						if( temps.size() > 0 && data <= temps.back() )
							d.errorAbort( "temperatures must be monotonically increasing" );
						temps.emplace_back( data );
						if( DEBUGSTATE )
							fprintf(ioQQQ,"%.2e\t", data);
					}
					if( DEBUGSTATE )
						fprintf(ioQQQ,"\n");

					numpoints = (int)temps.size();
					ASSERT( numpoints > 0 );
				}
				else if( dataType == "CS" || dataType == "RATE" )
				{
					bool isRate = ( dataType == "RATE" );
					string colliderType;
					d.getToken( colliderType );

					if( colliderType == "ELECTRON" )
						ipCollider = ipELECTRON;
					else if( colliderType == "PROTON" || colliderType == "H+" )
						ipCollider = ipPROTON;
					else if( colliderType == "HE+2" )
						ipCollider = ipALPHA;
					else if( colliderType == "HE+" )
						ipCollider = ipHE_PLUS;
					else if( colliderType == "H2ORTHO" )
						ipCollider = ipH2_ORTHO;
					else if( colliderType == "H2PARA" )
						ipCollider = ipH2_PARA;
					else if( colliderType == "H2" )
						ipCollider = ipH2;
					else if( colliderType == "HE" )
						ipCollider = ipATOM_HE;
					else if( colliderType == "H" )
						ipCollider = ipATOM_H;
					else
						d.errorAbort( "invalid type of collider for RATE, I know about "
							"ELECTRON PROTON H+ HE+2 HE+ H2ORTHO H2PARA H2 HE H" );

					if( !isRate && ipCollider != ipELECTRON )
						d.errorAbort( "collision strengths are only allowed for electron colliders" );

					if( temps.empty() )
						d.errorAbort( "must specify temperatures before the collision strengths" );

					long ipLoInFile, ipHiInFile, ipLo, ipHi;
					d.getToken( ipLoInFile );
					d.getToken( ipHiInFile );
					processIndices(ipLoInFile, ipHiInFile, lgIsRegular, indexold2new, ipLo, ipHi);

					if( ipHi >= nMolLevs )
						continue;

					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"<%s>\t%s\t%li\t%li:%li\t%li:%li",
							dBaseSpecies[intNS].chLabel, isRate?"RATE":"CS", ipCollider,
							ipLoInFile, ipHiInFile, ipLo+1, ipHi+1);
					}

					StoutCollData[intNS].lgIsRate(ipHi, ipLo, ipCollider) = isRate;

					ASSERT( numpoints > 0 );
					if( StoutCollData[intNS].ntemps(ipHi, ipLo, ipCollider) > 0 )
						d.errorAbort( "duplicate collisional transition found" );

					StoutCollData[intNS].setpoints(ipHi, ipLo, ipCollider, numpoints);

					for( int j = 0; j < numpoints; j++ )
					{
						StoutCollData[intNS].temps(ipHi, ipLo, ipCollider)[j] = temps[j];
						d.getToken( StoutCollData[intNS].collstrs(ipHi, ipLo, ipCollider)[j] );
						if( StoutCollData[intNS].collstrs(ipHi, ipLo, ipCollider)[j] <= 0. )
							d.errorAbort( "invalid collisional data" );
						if( DEBUGSTATE )
							fprintf(ioQQQ,"\t%.2e", StoutCollData[intNS].collstrs(ipHi, ipLo, ipCollider)[j]);
					}
					if( DEBUGSTATE )
						fprintf(ioQQQ,"\n");
				}
				else
				{
					d.errorAbort( "invalid data type" );
				}

				d.checkEOL();
			}
		}
	// ---------------- Version 2 COLLISION FORMAT (25 10 15) ----------------
		else
		{
	/*
	* V2 format: after magic, first non-comment data line is a header:
	*   "CS ELECTRON"  or  "RATE H"  etc.
	* Then TEMP lines + transition lines.
	*/
			d.getline();
			string dataTypeHeader;
			d.getToken( dataTypeHeader );

			bool isRate;
			if( dataTypeHeader == "CS" )
				isRate = false;
			else if( dataTypeHeader == "RATE" )
				isRate = true;
			else
				d.errorAbort( "Invalid header data type in STOUT .coll file" );

			string colliderType;
			d.getToken( colliderType );

			if( colliderType == "ELECTRON" )
				ipCollider = ipELECTRON;
			else if( colliderType == "PROTON" || colliderType == "H+" )
				ipCollider = ipPROTON;
			else if( colliderType == "HE+2" )
				ipCollider = ipALPHA;
			else if( colliderType == "HE+" )
				ipCollider = ipHE_PLUS;
			else if( colliderType == "H2ORTHO" )
				ipCollider = ipH2_ORTHO;
			else if( colliderType == "H2PARA" )
				ipCollider = ipH2_PARA;
			else if( colliderType == "H2" )
				ipCollider = ipH2;
			else if( colliderType == "HE" )
				ipCollider = ipATOM_HE;
			else if( colliderType == "H" )
				ipCollider = ipATOM_H;
			else
				d.errorAbort( "Invalid collider type in STOUT .coll header" );

			if( !isRate && ipCollider != ipELECTRON )
				d.errorAbort( "Collision strengths (CS) are only allowed for electron colliders" );

			lgSentinelReached = false;
			bool lgEarlyStop = false;
			long prevHiInFile = -1;
			long prevLoInFile = -1;
			numpoints = 0; /* require TEMP before transition rows */

			while( d.getline() )
			{
				if( d.lgEODMarker() )
				{
					lgSentinelReached = true;
					break;
				}

				string first;
				d.getToken( first );

				if( first == "TEMP" )
				{
					if( DEBUGSTATE )
						fprintf(ioQQQ,"<%s>\tTEMP\t", dBaseSpecies[intNS].chLabel);

					temps.clear();
					double t;
					while( d.getTokenOptional( t ) )
					{
						if( t <= 0. )
							d.errorAbort( "invalid temperature" );
						if( temps.size() > 0 && t <= temps.back() )
							d.errorAbort( "temperatures must be monotonically increasing" );
						temps.emplace_back( t );
						if( DEBUGSTATE )
							fprintf(ioQQQ,"%.2e\t", t);
					}
					if( DEBUGSTATE )
						fprintf(ioQQQ,"\n");

					numpoints = (int)temps.size();
					ASSERT( numpoints > 0 );
				}
				else
				{
					if( numpoints <= 0 )
						d.errorAbort( "First data line after header must be TEMP in v2 .coll file" );
					char* endptr = nullptr;
					errno = 0;
					long ipLoInFile = strtol( first.c_str(), &endptr, 10 );
					if( errno != 0 || endptr == first.c_str() || *endptr != '\0' )
						d.errorAbort( "failed to read lower level index in STOUT v2 .coll file" );

					long ipHiInFile, ipLo, ipHi;
					d.getToken( ipHiInFile );
					if( ipHiInFile < prevHiInFile ||
					    ( ipHiInFile == prevHiInFile && ipLoInFile < prevLoInFile ) )
					{
						d.errorAbort( "STOUT v2 .coll file is not sorted by (upper, lower) level index" );
					}
					prevHiInFile = ipHiInFile;
					prevLoInFile = ipLoInFile;

					if( ipHiInFile > maxRetainedIndexInFile )
					{
						lgEarlyStop = true;
						break;
					}
					processIndices(ipLoInFile, ipHiInFile, lgIsRegular, indexold2new, ipLo, ipHi);

					/* out-of-range/unmapped transitions: consume values and skip */
					if( ipLo == LONG_MAX || ipHi == LONG_MAX || ipHi >= nMolLevs || ipLo >= nMolLevs )
					{
						double dummy;
						for( int j=0; j<numpoints; ++j )
							d.getToken( dummy );
						d.checkEOL();
						continue;
					}

					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"<%s>\t%s\t%li\t%li:%li\t%li:%li",
							dBaseSpecies[intNS].chLabel, isRate?"RATE":"CS", ipCollider,
							ipLoInFile, ipHiInFile, ipLo+1, ipHi+1);
					}

					StoutCollData[intNS].lgIsRate(ipHi, ipLo, ipCollider) = isRate;

					ASSERT( numpoints > 0 );
					int oldNpts = StoutCollData[intNS].ntemps(ipHi, ipLo, ipCollider);

					if( oldNpts > 0 )
					{
						/* Baseline block: duplicates are an error; higher blocks may override */
						if( !lgOverrideExisting )
							d.errorAbort( "duplicate collisional transition found" );
					}

					/* First time or override: allocate for this transition */
					StoutCollData[intNS].setpoints(ipHi, ipLo, ipCollider, numpoints);

					for( int j = 0; j < numpoints; j++ )
					{
						StoutCollData[intNS].temps(ipHi, ipLo, ipCollider)[j] = temps[j];
						d.getToken( StoutCollData[intNS].collstrs(ipHi, ipLo, ipCollider)[j] );
						if( StoutCollData[intNS].collstrs(ipHi, ipLo, ipCollider)[j] <= 0. )
							d.errorAbort( "invalid collisional data" );
						if( DEBUGSTATE )
							fprintf(ioQQQ,"\t%.2e", StoutCollData[intNS].collstrs(ipHi, ipLo, ipCollider)[j]);
					}
					if( DEBUGSTATE )
						fprintf(ioQQQ,"\n");
				}

				d.checkEOL();
			}
			if( !lgSentinelReached && !lgEarlyStop )
			{
				fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n", chCOLLFilename.c_str() );
				fprintf( ioQQQ, " STOUT v2 .coll expects stars (*****) at end-of-data.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		} // end v2-coll branch
	} // end loop over collToProcess

} // end atmdat_STOUT_readin





/**
 * @brief Reads and processes CHIANTI atomic data files for a given species.
 *
 * This function reads the CHIANTI data files (.elvlc, .wgfa, .splups, .psplups) for a specified species,
 * parses the energy levels, transition probabilities, and collisional data, and initializes the corresponding
 * data structures used by the code. It handles various file formats, checks for data consistency, and supports
 * special cases (e.g., Fe species with more levels). The function also performs error checking and outputs
 * debug information if enabled.
 *
 * @param intNS
 *    The index of the species in the dBaseSpecies array.
 * @param chPrefix
 *    The prefix string used to construct the filenames for the CHIANTI data files.
 *
 * @details
 * The function performs the following steps:
 *   - Reads and validates the energy levels file (.elvlc), sorts levels, and initializes state arrays.
 *   - Reads the transition probability file (.wgfa), processes radiative data, and populates transition arrays.
 *   - Reads the electron and proton collision data files (.splups, .psplups), processes collisional strengths/rates,
 *     and fills collisional data arrays.
 *   - Handles special cases for Fe species and user-specified level limits.
 *   - Performs extensive error checking and outputs debug information if DEBUGSTATE is enabled.
 *
 * @throws
 *   Exits the program with an error message if any file is malformed, missing required data,
 *   or contains invalid entries.
 */
void atmdat_CHIANTI_readin( long intNS, const string& chPrefix )
{
	DEBUG_ENTRY( "atmdat_CHIANTI_readin()" );
	const bool DEBUGSTATE = false;

	int intsplinepts,intTranType,intxs;

	long int nLevelsUsed;
	// total number of chianti levels
	long int nTotalLevels;
	// number of experimental Chianti levels
	long int nExperimentalLevels;
	// number of energy levels with either theoretical or experimental energy
	long int nTheoreticalLevels;
	FILE *ioElecCollData=NULL, *ioProtCollData=NULL;
	realnum  fstatwt,fenergyWN,fWLAng,feinsteina;
	double fScalingParam,fEnergyDiff,EnergyTheory;
	const char chCommentChianti = '#';

	bool lgProtonData=false;

	// this is the largest number of levels allowed by the chianti format, I3
	const int MAX_NUM_LEVELS = 999;

	dBaseSpecies[intNS].lgMolecular = false;
	dBaseSpecies[intNS].lgLTE = false;

	string chUnCaps = chPrefix;
	uncaps(chUnCaps);

	string chEnFilename = chUnCaps;
	string chTraFilename = chUnCaps;
	string chEleColFilename = chUnCaps;
	string chProColFilename = chUnCaps;

	/*For the CHIANTI DATABASE*/
	/*Open the energy levels file*/
	chEnFilename += ".elvlc";

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEnFilename.c_str());

	fstream elvlcstream;
	open_data( elvlcstream, chEnFilename, mode_r );

	/*Open the transition probabilities file*/
	chTraFilename += ".wgfa";

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chTraFilename.c_str());

	fstream wgfastream;
	open_data( wgfastream, chTraFilename, mode_r );

	/*Open the electron collision data*/
	chEleColFilename += ".splups";

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEleColFilename.c_str());

	ioElecCollData = open_data( chEleColFilename, "r" );

	/*Open the proton collision data*/
	chProColFilename += ".psplups";

	/*Open the files*/
	if( DEBUGSTATE )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chProColFilename.c_str());

	/*We will set a flag here to indicate if the proton collision strengths are available */
	if( ( ioProtCollData = open_data( chProColFilename, "r", AS_TRY ) ) != NULL )
	{
		lgProtonData = true;
	}
	else
	{
		lgProtonData = false;
	}

	/*Loop finds how many theoretical and experimental levels are in the elvlc file */
	//eof_col is used get the first 4 charcters per line to find end of file
	const int eof_col = 5;
	//length (+1) of the nrg in the elvlc file
	const int lvl_nrg_col=16;
	//# of columns skipped from the left to get to nrg start
	const int lvl_skipto_nrg = 40;
	/* # of columns to skip from eof check to nrg start */
	const int lvl_eof_to_nrg = lvl_skipto_nrg - eof_col + 1;
	//# of columns to skip over the rydberg energy, we don't use it
	const int lvl_skip_ryd = 15;

	nTotalLevels = 0;
	/* level energy must be positive to be counted since 0 entered when no data.
	 * ground has zero energy so is not counted, as coded, so we start from 1 
	 * to compensate */
	nExperimentalLevels = 1;
	nTheoreticalLevels = 1;

	double EnergyExperimental = 0.;
	if (elvlcstream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		char exptemp[lvl_nrg_col],theotemp[lvl_nrg_col];
		/*This loop counts the number of valid rows within the elvlc file
		  as well as the number of experimental energy levels.*/
		while(nj != -1)
		{
			// count total number of energy level lines of data
			elvlcstream.get(otemp,eof_col);
			nj = atoi(otemp);
			if( nj == -1 )
				break;
			nTotalLevels++;

			// count number of experimental energy levels
			elvlcstream.seekg(lvl_eof_to_nrg,ios::cur);
			elvlcstream.get(exptemp,lvl_nrg_col);
			EnergyExperimental = (double) atof(exptemp);
			if( EnergyExperimental != 0. )
				nExperimentalLevels++;

			// count number of theoretical enerby levels
			elvlcstream.seekg(lvl_skip_ryd,ios::cur);
			elvlcstream.get(theotemp,lvl_nrg_col);
			EnergyTheory = (double) atof(theotemp);
			if( EnergyTheory != 0. || EnergyExperimental != 0. )
				nTheoreticalLevels++;

			elvlcstream.ignore(INT_MAX,'\n');
		}
		elvlcstream.seekg(0,ios::beg);
	}
	if( DEBUGSTATE )
		dprintf(ioQQQ,"CHIANTI in scope nTotalLevels %ld nExperimentalLevels %ld nTheoreticalLevels %ld \n",
			nTotalLevels,	nExperimentalLevels,nTheoreticalLevels);


	/* Sometimes the theoretical chianti level data is incomplete.
	 * If it is bad use experimental instead. 
	 * 25 06 12 This print never happens so perhaps Chianti is now complete? 
	 * 25 06 14 supplemental, text was "warning" so woul not be picked up by our tsuite scripts */

	long HighestIndexInFile = -1;

	/* total number of levels depends on the case; experiment, theory, or mixed.
	 * Previous test passed, so we have theory for every level. Mixed and Theory will
	 * both use all levels. Experiment does not use theory  */
	switch (atmdat.ChiantiType) {
		case t_atmdat::CHIANTI_EXP:
			HighestIndexInFile = nExperimentalLevels;
			if( DEBUGSTATE )
				fprintf(ioQQQ,"DEBUGG CHIANTI case EXP tot=%ld exp=%ld theo=%ld high indx=%ld\n", 	
				nTotalLevels,	nExperimentalLevels,nTheoreticalLevels,HighestIndexInFile);
			break;
		case t_atmdat::CHIANTI_THEO:
			HighestIndexInFile = nTheoreticalLevels;
			if( DEBUGSTATE )
				fprintf(ioQQQ,"DEBUGG CHIANTI case THEO tot=%ld exp=%ld theo=%ld high indx=%ld\n", 	
				nTotalLevels,	nExperimentalLevels,nTheoreticalLevels,HighestIndexInFile);
			break;
		case t_atmdat::CHIANTI_MIXED:
			// mixed option, we use everything we have but prefer experimental. We have theory for all levels
			HighestIndexInFile = nTheoreticalLevels;
			if( DEBUGSTATE )
				fprintf(ioQQQ,"DEBUGG CHIANTI case MIXED tot=%ld exp=%ld theo=%ld high indx=%ld\n", 	
				nTotalLevels,	nExperimentalLevels,nTheoreticalLevels,HighestIndexInFile);
			break;
		default:
			TotalInsanity();
	}
	dBaseSpecies[intNS].numLevels_max = HighestIndexInFile;

	setProperties(dBaseSpecies[intNS]);
	// derive default number of levels, Fe is special case, identify it
	if( tolower(dBaseSpecies[intNS].chLabel[0]) == 'f' && tolower(dBaseSpecies[intNS].chLabel[1]) == 'e')
	{
		// Fe is special case with more levels
		nLevelsUsed = MIN3(HighestIndexInFile, atmdat.nChiantiMaxLevelsFe,MAX_NUM_LEVELS );
	}
	else
	{
		nLevelsUsed = MIN3(HighestIndexInFile, atmdat.nChiantiMaxLevels,MAX_NUM_LEVELS );
	}

	if( nLevelsUsed <= 0 )
	{
		fprintf( ioQQQ, "WARNING: The number of energy levels is non-positive in datafile %s.\n", chEnFilename.c_str() );
		fprintf( ioQQQ, "The file must be corrupted.\n" );
		cdEXIT( EXIT_FAILURE );
	}

	//Consider the number of levels spceified on the masterlist. =1 if not specified
	long numMasterlist = MIN2( dBaseSpecies[intNS].numLevels_masterlist , HighestIndexInFile );
	nLevelsUsed = MAX2(nLevelsUsed,numMasterlist);
		if( DEBUGSTATE )
			fprintf(ioQQQ,"DEBUGG CHIANTI levels post masterlist %ld used of total %ld masterlist=%ld exp=%ld theo=%ld HighestIndexInFile=%ld\n", 	
				nLevelsUsed, nTotalLevels,numMasterlist,	nExperimentalLevels,nTheoreticalLevels,
				HighestIndexInFile);

	if (dBaseSpecies[intNS].setLevels != -1)
	{
		if (dBaseSpecies[intNS].setLevels > HighestIndexInFile)
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel );
			if( dBaseSpecies[intNS].setLevels != LONG_MAX )
				fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li requested,"
					 " only %li energy levels available.\n",
					 dBaseSpecies[intNS].chLabel, chLabelChemical, dBaseSpecies[intNS].setLevels,
					 HighestIndexInFile );
			nLevelsUsed = HighestIndexInFile;
		}
		else
		{
			nLevelsUsed = dBaseSpecies[intNS].setLevels;
		}
	}

	dBaseSpecies[intNS].numLevels_max = nLevelsUsed;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	if( atmdat.lgChiantiPrint )
	{
		char chLabelChemical[CHARS_SPECIES] = "";
		switch( atmdat.ChiantiType )
		{
			case t_atmdat::CHIANTI_EXP:
				spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
				fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li experimental energy levels of %li available.\n",
					dBaseSpecies[intNS].chLabel, chLabelChemical, nLevelsUsed , nExperimentalLevels );
				break;
			case t_atmdat::CHIANTI_MIXED:
				spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
				fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li mixed energy levels of %li available.\n",
					dBaseSpecies[intNS].chLabel, chLabelChemical, nLevelsUsed , nTheoreticalLevels );
				break;
			case t_atmdat::CHIANTI_THEO:
				spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
				fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li theoretical energy levels of %li available.\n",
					dBaseSpecies[intNS].chLabel, chLabelChemical, nLevelsUsed , nTheoreticalLevels );
				break;
			default:
				TotalInsanity();
		}
	}

	/* allocate the States array*/
	dBaseStates[intNS].init(dBaseSpecies[intNS].chLabel,nLevelsUsed);

	/* allocate the Transition array*/
	ipdBaseTrans[intNS].reserve(nLevelsUsed);
	for( long ipHi = 1; ipHi < nLevelsUsed; ipHi++)
		ipdBaseTrans[intNS].reserve(ipHi,ipHi);
	ipdBaseTrans[intNS].alloc();
	dBaseTrans[intNS].resize(ipdBaseTrans[intNS].size());
	dBaseTrans[intNS].states() = &dBaseStates[intNS];
	dBaseTrans[intNS].chLabel() = dBaseSpecies[intNS].chLabel;
	dBaseSpecies[intNS].database = "Chianti";

	int ndBase = 0;
	for( long ipHi = 1; ipHi < nLevelsUsed; ipHi++)
	{
		for( long ipLo = 0; ipLo < ipHi; ipLo++)
		{
			ipdBaseTrans[intNS][ipHi][ipLo] = ndBase;
			dBaseTrans[intNS][ndBase].Junk();
			dBaseTrans[intNS][ndBase].setLo(ipLo);
			dBaseTrans[intNS][ndBase].setHi(ipHi);
			++ndBase;
		}
	}

	/* Keep track of which levels have experimental data and then create a vector
	 * which relates their indices to the default chianti energy indices.
	 */
	long ncounter = 0;

	//Relate Chianti level indices to a set that only include experimental levels
	vector<long> intExperIndex(nTotalLevels,-1);

	DoubleLongPairVector dBaseStatesEnergy;
	vector<double> dBaseStatesStwt(HighestIndexInFile,-1.0);
	for( long ii = 0; ii < HighestIndexInFile; ii++ )
	{
		dBaseStatesEnergy.push_back(make_pair(-1.0,ii));
	}

	//lvl_skipto_statwt is the # of columns to skip to statwt from left
	const int lvl_skipto_statwt = 37;
	//lvl_statwt_col is the length (+1) of statwt
	const int lvl_statwt_col = 4;
	//Read in stat weight and energy

	//lvl_skip_to_exp_nrg is the # of columns to skip to experimental value from lvl_skipto_statwt + lvl_statwt_col
	const int lvl_skip_to_exp_nrg = 3;
	//lvl_skip_to_exp_nrg is the # of columns to skip to experimental value from lvl_skip_to_exp_nrg + lvl_nrg_col
	const int lvl_skip_to_theo_nrg = 13;

	//Read in nrg levels to see if they are in order
	for( long ipLev=0; ipLev<nTotalLevels; ipLev++ )
	{
		if(elvlcstream.is_open())
		{
			char gtemp[lvl_statwt_col],theotemp[lvl_nrg_col],exptemp[lvl_nrg_col];
			//char gtemp[lvl_statwt_col],thtemp[lvl_nrg_col],obtemp[lvl_nrg_col],theotemp[lvl_nrg_col],exptemp[lvl_nrg_col];
			elvlcstream.seekg(lvl_skipto_statwt,ios::cur);
			elvlcstream.get(gtemp,lvl_statwt_col);
			fstatwt = (realnum) atof(gtemp)	;

			//elvlcstream.get(thtemp,lvl_nrg_col);
			//fenergy = (double) atof(thtemp);

			//Reading experimental column
			elvlcstream.seekg(lvl_skip_to_exp_nrg,ios::cur);
			elvlcstream.get(exptemp,lvl_nrg_col);
			EnergyExperimental = (double) atof(exptemp);

			// Reading Theoretical Column
			elvlcstream.seekg(lvl_skip_to_theo_nrg,ios::cur);
			elvlcstream.get(theotemp,lvl_nrg_col);
			EnergyTheory = (double) atof(theotemp);
			//dprintf(ioQQQ,"theory 2nd energy exp theo %f\n",EnergyExperimental,EnergyTheory);

			if(fstatwt <= 0.)
			{
				fprintf( ioQQQ, " WARNING: A positive non zero value is expected for the "
						"statistical weight but was not found in %s"
						" level %li\n", chEnFilename.c_str(),ipLev);
				cdEXIT(EXIT_FAILURE);
			}

			/* Go through the entire level list selectively choosing experimental or theoretical level energies.
			 * Store them, not zeroes, in order using ncounter to count the index.
			 * Any row on the level list where there is no valid energy, put a -1 in the relational vector.
			 * If it is a valid energy level store the new ncounter index.
			 */
			double EnergyUsed = 0.;
			switch( atmdat.ChiantiType )
			{
			case t_atmdat::CHIANTI_EXP:
				// only use experimental energies
				EnergyUsed = EnergyExperimental;
				break;
			case t_atmdat::CHIANTI_THEO:
				// prefer theoretical energy, but if that is absent use experimental energy
				if( EnergyTheory != 0. )
					EnergyUsed = EnergyTheory;
				else
					EnergyUsed = EnergyExperimental;
				break;
			case t_atmdat::CHIANTI_MIXED:
				// prefer experimental energy, but if that is absent use theoretical energy
				if( EnergyExperimental != 0. )
					EnergyUsed = EnergyExperimental;
				else
					EnergyUsed = EnergyTheory;
				break;
			default:
				TotalInsanity();
			}

			if( EnergyUsed != 0. || ipLev == 0 )
			{
				dBaseStatesEnergy.at(ncounter).first = EnergyUsed;
				dBaseStatesEnergy.at(ncounter).second = ncounter;
				dBaseStatesStwt.at(ncounter) = fstatwt;
				intExperIndex.at(ipLev) = ncounter;
				ncounter++;
			}
			else
			{
				intExperIndex.at(ipLev) = -1;
			}

			elvlcstream.ignore(INT_MAX,'\n');
		}
		else
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEnFilename.c_str());
			fclose( ioProtCollData );
			cdEXIT(EXIT_FAILURE);
		}
	}

	elvlcstream.close();

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nintExperIndex Vector:\n");
		fprintf(ioQQQ,"File Index|Exper Index\n");
		for( vector<long>::iterator i = intExperIndex.begin(); i != intExperIndex.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-intExperIndex.begin())+1;
			fprintf(ioQQQ,"%li\t%li\n",iPrt,(*i)+1);
		}

		for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-dBaseStatesEnergy.begin())+1;
			fprintf(ioQQQ,"PreSort:%li\t%li\t%f\t%f\n",iPrt,
					(i->second)+1,i->first,dBaseStatesStwt.at(i->second));
		}
	}

	//Sort energy levels
	sort(dBaseStatesEnergy.begin(),dBaseStatesEnergy.end());

	std::vector<long> indexold2new(dBaseStatesEnergy.size());
	for( DoubleLongPairVector::const_iterator i = dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
	{
		indexold2new[i->second] = i-dBaseStatesEnergy.begin();
	}

	if( DEBUGSTATE )
	{
		for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-dBaseStatesEnergy.begin())+1;
			if( iPrt > nLevelsUsed )
				break;
			fprintf(ioQQQ,"PostSort:%li\t%li\t%f\t%f\n",iPrt,
					(i->second)+1,i->first,dBaseStatesStwt.at(i->second));
		}

		fprintf(ioQQQ,"\nChianti Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Energy Level File: %s\n",chEnFilename.c_str());
		if( atmdat.ChiantiType==t_atmdat::CHIANTI_EXP )
		{
			fprintf(ioQQQ,"Number of Experimental Energy Levels in File: %li\n",nExperimentalLevels);
		}
		else
		{
			fprintf(ioQQQ,"Number of Theoretical Energy Levels in File: %li\n",nTotalLevels);
		}
		fprintf(ioQQQ,"Number of Energy Levels Cloudy is Currently Using: %li\n",nLevelsUsed);
		fprintf(ioQQQ,"Species|File Index|Cloudy Index|StatWT|Energy(WN)\n");
	}

	vector<long> revIntExperIndex;
	revIntExperIndex.resize(dBaseStatesEnergy.size());
	for (size_t i = 0; i<dBaseStatesEnergy.size(); ++i)
		revIntExperIndex[i] = -1;
	for ( vector<long>::const_iterator i= intExperIndex.begin();
		  i != intExperIndex.end(); ++i )
	{
		long ipos = intExperIndex[i-intExperIndex.begin()];
		if (ipos >= 0 && ipos < long(dBaseStatesEnergy.size()))
			revIntExperIndex[ipos] = i-intExperIndex.begin();
	}

	for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
	{

		long ipLevNew = i - dBaseStatesEnergy.begin();
		long ipLevFile = -1;

		if( ipLevNew >= nLevelsUsed )
			break;

		if( atmdat.ChiantiType==t_atmdat::CHIANTI_EXP )
		{
			ipLevFile = revIntExperIndex[ipLevNew];
		}
		else
		{
			ipLevFile = i->second;
		}

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"<%s>\t%li\t%li\t",dBaseSpecies[intNS].chLabel,ipLevFile+1,ipLevNew+1);
		}

		dBaseStates[intNS][ipLevNew].g() = dBaseStatesStwt.at(i->second);
		dBaseStates[intNS][ipLevNew].energy().set(i->first,"cm^-1");
		dBaseStates[intNS][ipLevNew].ipOrg() = ipLevFile+1;

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"%.1f\t",dBaseStatesStwt.at(i->second));
			fprintf(ioQQQ,"%.3f\n",i->first);
		}
	}

	// highest energy transition in chianti
	dBaseSpecies[intNS].maxWN = 0.;
	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for(TransitionList::iterator tr=dBaseTrans[intNS].begin();
		 tr!= dBaseTrans[intNS].end(); ++tr)
	{
		int ipHi = tr->ipHi();
		int ipLo = tr->ipLo();
		fenergyWN = (realnum) dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN();

		tr->EnergyWN() = fenergyWN;

		if( rfield.isEnergyBound( Energy( fenergyWN, "cm^-1" ) ) )
		{
			tr->WLangVac() = wn2angVac( fenergyWN );
			dBaseSpecies[intNS].maxWN = MAX2(dBaseSpecies[intNS].maxWN,fenergyWN);
		}
		else
			tr->WLangVac() = 1e30;
	}

	/************************************************************************/
	/*** Read in the level numbers, Einstein A and transition wavelength  ***/
	/************************************************************************/

	//Count the number of rows first
	long wgfarows = -1;
	if (wgfastream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		while(nj != -1)
		{
			wgfastream.get(otemp,eof_col);
			wgfastream.ignore(INT_MAX,'\n');
			if( otemp[0] == chCommentChianti ) continue;
			nj = atoi(otemp);
			wgfarows++;
		}
		wgfastream.seekg(0,ios::beg);
	}
	else 
		fprintf( ioQQQ, "WARNING The data file %s is corrupted .\n",chTraFilename.c_str());


	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\n\nTransition Probability File: %s\n",chTraFilename.c_str());
		fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Wavelength(A)|Ein A\n");
	}


	//line_index_col is the length(+1) of the level indexes in the WGFA file
	const int line_index_col = 6;
	//line_nrg_to_eina is the # of columns to skip from wavelength to eina in WGFA file
	const int line_nrg_to_eina =15;
	//line_eina_col is the length(+1) of einsteinA in WGFA
	const int line_eina_col = 16;
	char lvltemp[line_index_col];
	//Start reading WGFA file
	if (wgfastream.is_open())
	{
		for (long ii = 0;ii<wgfarows;ii++)
		{
			wgfastream.get(lvltemp,line_index_col);
			if( lvltemp[0] == chCommentChianti )
			{
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			long ipLoInFile = atoi(lvltemp);
			wgfastream.get(lvltemp,line_index_col);
			long ipHiInFile = atoi(lvltemp);

			// ipLo and ipHi will be manipulated below, want to retain original vals for prints
			long ipLo = ipLoInFile - 1;
			long ipHi = ipHiInFile - 1;

			if( atmdat.ChiantiType==t_atmdat::CHIANTI_EXP )
			{
				/* If either upper or lower index is -1 in the relational vector,
				 * skip that line in the wgfa file.
				 * Otherwise translate the level indices.*/
				if( intExperIndex[ipLo] == -1 || intExperIndex[ipHi] == -1 )
				{
					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}
				else
				{
					ipHi = intExperIndex.at(ipHi);
					if (ipHi < long(indexold2new.size()))
					{
						ipHi = indexold2new[ipHi];
					}
					else
					{
						ipHi = -1;
					}
					ipLo = intExperIndex.at(ipLo);
					if (ipLo < long(indexold2new.size()))
					{
						ipLo = indexold2new[ipLo];
					}
					else
					{
						ipLo = -1;
					}
				}
			}
			else
			{
				long testlo = -1, testhi = -1;

				try
				{
					testlo = indexold2new[ipLo];
					testhi = indexold2new[ipHi];
				}
				catch ( out_of_range& /* e */ )
				{
					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"NOTE: An out of range exception has occurred"
								" reading in data from %s\n",chTraFilename.c_str());
						fprintf(ioQQQ," The line in the file containing the unidentifiable"
								" levels has been ignored.\n");
						fprintf(ioQQQ,"There is no reason for alarm."
								" This message is just for documentation.\n");
					}

					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}

				if(  testlo == -1 || testhi == -1 )
				{
					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}
				else
				{
					ipLo = testlo;
					ipHi = testhi;
				}
			}

			if( ipLo >= nLevelsUsed || ipHi >= nLevelsUsed )
			{
				// skip these lines
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			if( ipHi == ipLo )
			{
				fprintf( ioQQQ," WARNING: Upper level = lower for a radiative transition in %s. Ignoring.\n", chTraFilename.c_str() );
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"<%s>\t%li:%li\t%li:%li\t",dBaseSpecies[intNS].chLabel,ipLoInFile,ipHiInFile,ipLo+1,ipHi+1);
			}

			ASSERT( ipHi != ipLo );
			ASSERT( ipHi >= 0 );
			ASSERT( ipLo >= 0 );

			// sometimes the CHIANTI datafiles list the highest index first as in the middle of these five lines in ne_10.wgfa:
			//    ...
			//    8   10       187.5299      0.000e+00      4.127e+05                 3d 2D1.5 -                  4s 2S0.5           E2
			//    9   10       187.6573      0.000e+00      6.197e+05                 3d 2D2.5 -                  4s 2S0.5           E2
			//   11   10   4842624.0000      1.499e-05      9.423e-06                 4p 2P0.5 -                  4s 2S0.5           E1
			//    1   11         9.7085      1.892e-02      6.695e+11                 1s 2S0.5 -                  4p 2P0.5           E1
			//    2   11        48.5157      6.787e-02      9.618e+10                 2s 2S0.5 -                  4p 2P0.5           E1
			//    ...
			// so, just set ipHi (ipLo) equal to the max (min) of the two indices.
			// NB NB NB it looks like this may depend upon whether one uses observed or theoretical energies.

			//Read in wavelengths
			char trantemp[lvl_nrg_col];
			wgfastream.get(trantemp,lvl_nrg_col);
			fWLAng = (realnum)atof(trantemp);
			if( DEBUGSTATE && atmdat.ChiantiType==t_atmdat::CHIANTI_EXP)
			{
				fprintf(ioQQQ,"%.4f\t",fWLAng);
			}

			/* \todo 2 CHIANTI labels the H 1 2-photon transition as z wavelength of zero.
			 * Should we just ignore all of the wavelengths in this file and use the
			 * difference of level energies instead. */

			if( ipHi < ipLo )
			{
				long swap = ipHi;
				ipHi = ipLo;
				ipLo = swap;
			}

			/* If the given wavelength is negative, then theoretical energies are being used.
			 * Take the difference in stored theoretical energies.
			 * It should equal the absolute value of the wavelength in the wgfa file. */
			if( fWLAng <= 0. ) // && !atmdat.lgChiantiExp )
			{
				//if( fWLAng < 0.)
					//fprintf( ioQQQ," WARNING: Negative wavelength for species %6s, indices %3li %3li \n", dBaseSpecies[intNS].chLabel, ipLo, ipHi);
				fWLAng = (realnum)(1e8/abs(dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN()));
			}

			if( DEBUGSTATE && !(atmdat.ChiantiType==t_atmdat::CHIANTI_EXP))
			{
				fprintf(ioQQQ,"%.4f\t",fWLAng);
			}
			//Skip from end of Wavelength to Einstein A and read in
			wgfastream.seekg(line_nrg_to_eina,ios::cur);
			wgfastream.get(trantemp,line_eina_col);
			feinsteina = (realnum)atof(trantemp);
			if( feinsteina == 0. )
			{
				static bool notPrintedYet = true;
				if( notPrintedYet && atmdat.lgChiantiPrint)
				{
					fprintf( ioQQQ," CAUTION: Radiative rate(s) equal to zero in %s.\n", chTraFilename.c_str() );
					notPrintedYet = false;
				}
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}
			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"%.3e\n",feinsteina);
			}

			//Read in the rest of the line and look for auto
			string chLine;
			getline( wgfastream, chLine );
			TransitionList::iterator tr = dBaseTrans[intNS].begin()+ipdBaseTrans[intNS][ipHi][ipLo];
			if( chLine.find("auto") != string::npos )
			{
				if( tr->hasEmis() )
				{
					tr->Emis().AutoIonizFrac() =
						feinsteina/(tr->Emis().Aul() + feinsteina);
					ASSERT( tr->Emis().AutoIonizFrac() >= 0. );
					ASSERT( tr->Emis().AutoIonizFrac() <= 1. );
				}
				continue;
			}

			if( tr->hasEmis() )
			{
				fprintf(ioQQQ," PROBLEM duplicate transition found by atmdat_chianti in %s, "
						"wavelength=%f\n", chTraFilename.c_str(),fWLAng);
				wgfastream.close();
				cdEXIT(EXIT_FAILURE);
			}
			tr->AddLine2Stack();
			tr->Emis().Aul() = feinsteina;

			fenergyWN = (realnum)(1e+8/fWLAng);

			// \todo Check the wavelength in the file with the difference in energy levels

			tr->EnergyWN() = fenergyWN;
			if( rfield.isEnergyBound( Energy( fenergyWN, "cm^-1" ) ) )
			{
				tr->WLangVac() = wn2angVac( fenergyWN );
				tr->Emis().gf() = (realnum)GetGF(tr->Emis().Aul(), tr->EnergyWN(), tr->Hi()->g());
			}
			else
				tr->WLangVac() = 1e30;

			tr->setComment( db_comment_tran_levels() );
		}
	}
	else 
		fprintf( ioQQQ, "WARNING  The data file %s is corrupted .\n",chTraFilename.c_str());
	wgfastream.close();

	/* allocate space for splines */
	AtmolCollSplines[intNS].reserve(nLevelsUsed);
	for( long ipHi=0; ipHi<nLevelsUsed; ipHi++ )
	{
		AtmolCollSplines[intNS].reserve(ipHi,nLevelsUsed);
		for( long ipLo=0; ipLo<nLevelsUsed; ipLo++ )
		{
			AtmolCollSplines[intNS].reserve(ipHi,ipLo,ipNCOLLIDER);
		}
	}
	AtmolCollSplines[intNS].alloc();

	for( long ipHi=0; ipHi<nLevelsUsed; ipHi++ )
	{
		for( long ipLo=0; ipLo<nLevelsUsed; ipLo++ )
		{
			for( long k=0; k<ipNCOLLIDER; k++ )
			{
				/* initialize all spline variables */
				AtmolCollSplines[intNS][ipHi][ipLo][k].nSplinePts = -1; 
				AtmolCollSplines[intNS][ipHi][ipLo][k].intTranType = -1;
				AtmolCollSplines[intNS][ipHi][ipLo][k].EnergyDiff = BIGDOUBLE;
				AtmolCollSplines[intNS][ipHi][ipLo][k].ScalingParam = BIGDOUBLE;
			}
		}
	}

	/************************************/
	/*** Read in the collisional data ***/
	/************************************/

	// ipCollider 0 is electrons, 1 is protons
	for( long ipCollider=0; ipCollider<=1; ipCollider++ )
	{
		string chFilename;

		if( ipCollider == ipELECTRON )
		{
			chFilename = chEleColFilename;
		}
		else if( ipCollider == ipPROTON )
		{
			if( !lgProtonData )
				break;
			chFilename = chProColFilename;
		}
		else
			TotalInsanity();

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"\n\nCollision Data File: %s\n",chTraFilename.c_str());
			fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Spline Points\n");
		}

		fstream splupsstream;
		open_data( splupsstream, chFilename, mode_r );

		//cs_eof_col is the length(+1) of the first column used for finding the end of file
		const int cs_eof_col = 4;
		//cs_index_col is the length(+1) of the indexes in the CS file
		const int cs_index_col = 4;
		//cs_trantype_col is the length(+1) of the transition type in the CS file
		const int cs_trantype_col = 4;
		//cs_values_col is the length(+1) of the other values in the CS file
		//including: GF, nrg diff, scaling parameter, and spline points
		const int cs_values_col = 11;
		//Determine the number of rows in the CS file
		if (splupsstream.is_open())
		{
			int nj = 0;
			//splupslines is -1 since the loop runs 1 extra time
			long splupslines = -1;
			char otemp[cs_eof_col];
			while(nj != -1)
			{
				splupsstream.get(otemp,cs_eof_col);
				splupsstream.ignore(INT_MAX,'\n');
				nj = atoi(otemp);
				splupslines++;
			}
			splupsstream.seekg(0,ios::beg);

			for (int m = 0;m<splupslines;m++)
			{
				if( ipCollider == ipELECTRON )
				{
					splupsstream.seekg(6,ios::cur);
				}

				/* level indices */
				splupsstream.get(otemp,cs_index_col);
				long ipLo = atoi(otemp)-1;
				splupsstream.get(otemp,cs_index_col);
				long ipHi = atoi(otemp)-1;

				long ipLoFile = ipLo;
				long ipHiFile = ipHi;

				/* If either upper or lower index is -1 in the relational vector,
				* skip that line in the splups file.
				* Otherwise translate the level indices.*/
				if( atmdat.ChiantiType==t_atmdat::CHIANTI_EXP )
				{
					if( intExperIndex[ipLo] == - 1 || intExperIndex[ipHi] == -1 )
					{
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}
					else
					{
						ipHi = intExperIndex.at(ipHi);
						if (ipHi < long(indexold2new.size()))
						{
							ipHi = indexold2new[ipHi];
						}
						else
						{
							ipHi = -1;
						}
						ipLo = intExperIndex.at(ipLo);
						if (ipLo < long(indexold2new.size()))
						{
							ipLo = indexold2new[ipLo];
						}
						else
						{
							ipLo = -1;
						}
					}
				}
				else
				{
					long testlo = -1, testhi = -1;

					/* With level trimming on it is possible that there can be rows that
					 * have to be skipped when using theoretical
					 * since the levels no longer exist */
					try
					{
						testlo = indexold2new[ipLo];
						testhi = indexold2new[ipHi];
					}
					catch ( out_of_range& /* e */ )
					{
						if( DEBUGSTATE )
						{
							fprintf(ioQQQ,"NOTE: An out of range exception has occurred"
									" reading in data from %s\n",chEleColFilename.c_str());
							fprintf(ioQQQ," The line in the file containing the unidentifiable"
									" levels has been ignored.\n");
							fprintf(ioQQQ,"There is no reason for alarm."
									" This message is for documentation.\n");
						}
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}

					if( testlo == -1 || testhi == -1 )
					{
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}
					else
					{
						ipLo = testlo;
						ipHi = testhi;
					}
				}

				if( ipLo >= nLevelsUsed || ipHi >= nLevelsUsed )
				{
					// skip these transitions
					splupsstream.ignore(INT_MAX,'\n');
					continue;
				}

				if( ipHi < ipLo )
				{
					long swap = ipHi;
					ipHi = ipLo;
					ipLo = swap;
				}

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"<%s>\t%li:%li\t%li:%li",dBaseSpecies[intNS].chLabel,ipLoFile+1,ipHiFile+1,ipLo+1,ipHi+1);
				}

				/*Transition Type*/
				splupsstream.get(otemp,cs_trantype_col);
				intTranType = atoi(otemp);
				char qtemp[cs_values_col];
				splupsstream.get(qtemp,cs_values_col);
				/*Energy Difference*/
				splupsstream.get(qtemp,cs_values_col);
				fEnergyDiff = atof(qtemp);
				/*Scaling Parameter*/
				splupsstream.get(qtemp,cs_values_col);
				fScalingParam = atof(qtemp);

				ASSERT( ipLo != ipHi );
				ASSERT( ipLo >= 0 && ipLo < nLevelsUsed );
				ASSERT( ipHi >= 0 && ipHi < nLevelsUsed );

				while( splupsstream.peek() != '\n' )
				{
					splupsstream.get(qtemp,cs_values_col);
					if( qtemp[0] == ' ' && qtemp[1] == ' ' )
						break;
					double temp = atof(qtemp);
					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"\t%.3e",temp);
					}
					// intTranType == 6 means log10 of numbers have been fit => allow negative numbers
					// intTranType < 6 means linear numbers have been fit => negative numbers are unphysical
					if( intTranType < 6 )
						temp = max( temp, 0. );
					AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline.push_back( temp );
				}

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"\n");
				}

				intsplinepts = int( AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline.size() );
				ASSERT( intsplinepts > 2 );

				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer.resize( intsplinepts );

				/*The zeroth element contains the number of spline points*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].nSplinePts = intsplinepts;
				/*Transition type*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].intTranType = intTranType;
				/*Energy difference between two levels in Rydbergs*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].EnergyDiff = fEnergyDiff;
				/*Scaling parameter C*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].ScalingParam = fScalingParam;

				/*Once the spline points have been filled,fill the second derivatives structure*/
				/*Creating spline points array*/
				vector<double> xs (intsplinepts),
					spl(intsplinepts),
					spl2(intsplinepts);

				double coeff = (double)1/(intsplinepts-1);
				for(intxs=0;intxs<intsplinepts;intxs++)
				{
					xs[intxs] = coeff*intxs;
					spl[intxs] = AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intxs];
				}

				spline(&xs[0], &spl[0],intsplinepts,2e31,2e31,&spl2[0]);

				/*Filling the second derivative structure*/
				for( long ii=0; ii<intsplinepts; ii++)
				{
					AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer[ii] = spl2[ii];
				}

				splupsstream.ignore(INT_MAX,'\n');
			}
			splupsstream.close();
		}
	}

	// close open file handles
	fclose( ioElecCollData );
	if( lgProtonData )
		fclose( ioProtCollData );

	return;
}

