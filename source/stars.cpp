/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "optimize.h"
#include "continuum.h"
#include "called.h"
#include "rfield.h"
#include "stars.h"
#include "container_classes.h"

/** this is the initial assumed size of the Starburst grid, may be increased during execution if needed */
static const int NSB99 = 1250;
/** maximum number of separate time steps in a Starburst99 model */
static const int MNTS = 200;

/** this is the number of points in each of the stellar continua */
static const int NRAUCH = 19951;
/** The number of models in the original Rauch H-Ca set (first version May 1998, current May 2001) */
static const int NMODS_HCA = 66;
/** The number of models in the new Rauch H-Ni set, Nov 2002 */
static const int NMODS_HNI = 51;
/** The number of models in the new Rauch PG1159 set, Jan 2006 */
static const int NMODS_PG1159 = 71;
/** The number of models in the Rauch Hydrogen only set, Feb 2003 */
static const int NMODS_HYDR = 100;
/** The number of models in the Rauch Helium only set, Jun 2004 */
static const int NMODS_HELIUM = 81;
/** The number of models in the Rauch H+He set, Aug 2004 */
static const int NMODS_HpHE = 117;

/* set to true to turn on debug print statements in these routines */
static const bool DEBUGPRT = false;

static const bool lgSILENT = false;
static const bool lgVERBOSE = true;

static const bool lgLINEAR = false;
static const bool lgTAKELOG = true;

static const bool lgREAD_BIN = false;
static const bool lgREAD_ASCII = true;

typedef enum {
	ISB_COLLECT, ISB_EXECUTE, ISB_FIRST, ISB_SECOND
} IntStageBits;

static const int IS_COLLECT = 2<<ISB_COLLECT;
static const int IS_EXECUTE = 2<<ISB_EXECUTE;
static const int IS_FIRST =   2<<ISB_FIRST;
static const int IS_SECOND =  2<<ISB_SECOND;

/** store the parameters of a single atmosphere model */
struct mpp
{
	double par[MDIM];
	int modid;
	char chGrid;
	mpp() { memset(this, 0, sizeof(mpp)); }
};

/** \todo - check rebinning of Tlusty models
 ** \todo - why was it necessary to change stars_tlusty.in? (change from r43 to r50?)
 ** \todo - check all interpolation modes of CoStar
 ** \todo - compare models with original code, dump atmospheres!
 */

/* this is the structure of the binary atmosphere file (VERSION 20100902[01]):
 *
 *               ============================
 *               * int32 VERSION            *
 *               * int32 MDIM               *
 *               * int32 MNAM               *
 *               * int32 ndim               *
 *               * int32 npar               *
 *               * int32 nmods              *
 *               * int32 ngrid              *
 *               * uint32 nOffset           *
 *               * uint32 nBlocksize        *
 *               * double mesh_elo          *
 *               * double mesh_ehi          *
 *               * double mesh_res_factor   *
 *               * char md5sum[NMD5]        *
 *               * char names[MDIM][MNAM+1] *
 *               * mpp telg[nmods]          *
 *               * realnum anu[ngrid]       *
 *               * realnum mod1[ngrid]      *
 *               *    ...                   *
 *               * realnum modn[ngrid]      *
 *               ============================
 *
 * nOffset == 7*sizeof(int32) + 2*sizeof(uint32) + 3*sizeof(double) +
 *            (NMD5 + MDIM*(MNAM+1))*sizeof(char) + nmods*sizeof(mpp)
 * nBlocksize == ngrid*size(realnum) */

/** store all the relevant information on a binary atmosphere file */
struct stellar_grid
{
	/** the name of the atmosphere file */
	string name;
	/** if true, more relaxed rules for matching log(g) will be used */
	bool lgIsTeffLoggGrid;
	/** where should we search for the binary atmosphere file */
	access_scheme scheme;
	/** the file handle for this file */
	FILE *ioIN;
	/** the identifier for this grid used in the Cloudy output,
	 * this *must* be exactly 12 characters long */
	string ident;
	/** the Cloudy command to recompile the binary atmosphere file */
	string command;
	/** which interpolation mode is requested */
	IntMode imode;
	/** the number of dimensions in the grid */
	int32 ndim;
	/** the number of parameters for each model; npar >= ndim */
	int32 npar;
	/** the number of stellar atmosphere models in this file */
	int32 nmods;
	/** the number of grid points per model, should equal rfield.nflux_with_check */
	int32 ngrid;
	/** the offset to the first data block (the anu grid) */
	uint32 nOffset;
	/** the size of each model block in bytes */
	uint32 nBlocksize;
	/** these are the model parameters in the same
	 * sequence they are stored in the binary file */
	vector<mpp> telg;    /* telg[nmods] */
	/** these are the unique values for each of the model parameters */
	multi_arr<double,2> val; /* val[ndim][nval[n]] */
	/** nval[n] is the number of unique values in val[n][*] */
	vector<long> nval;   /* nval[ndim] */
	/** jlo/jhi will hold indices into the binary model file: jlo/jhi(i,...,n)
	 * will point to the model with parameters val[0][i],...,val[ndim-1][n],
	 * or its closest approximation in log(g) in case the model doesn't exist
	 * and lgIsTeffLoggGrid is true.
	 * jlo will hold the model with the highest log(g) <= than requested
	 * jhi will hold the model with the lowest log(g) >= than requested
	 * in case no suitable model could be found either array will hold -2 */
	vector<long> jlo;   /* jlo(nval[0],...,nval[ndim-1]) */
	vector<long> jhi;   /* jhi(nval[0],...,nval[ndim-1]) */
	/** this array will hold the designations for each dimension of the grid */
	char names[MDIM][MNAM+1];
	/** this array holds the length of each CoStar track */
	vector<long> trackLen; /* trackLen[nTracks] */
	/** this is the number of CoStar tracks */
	long nTracks;
	/** jval will hold indices into the CoStar grid: jval(nModels,nTracks) */
	vector<long> jval;
	/** set to true if we read directly from the ascii file */
	bool lgASCII;
	/** helper variables for reading the ascii atmosphere files */
	double convert_wavl;
	double convert_flux;
	bool lgFreqX;
	bool lgFreqY;
	/** cautions generated by the interpolation routine */
	map<string,int> caution;
	/** the list of SEDs from the grid needed for the interpolation */
	vector<long> index_list, index_list2;
	/** array for holding the rebinned SEDs from the grid */
	multi_arr<realnum,2> CloudyFlux;
	stellar_grid()
	{
		ioIN = NULL;
		memset( names, '\0', MDIM*(MNAM+1) );
	}
	~stellar_grid()
	{
		if( ioIN != NULL )
			fclose( ioIN );
	}
};

/* internal routines */
STATIC bool CoStarInitialize(const char[],const char[]);
STATIC void InterpolateGridCoStar(stellar_grid*,const double[],double*,double*);
STATIC void FindHCoStar(const stellar_grid*,long,double,long,vector<realnum>&,vector<long>&,vector<long>&);
STATIC void FindVCoStar(const stellar_grid*,double,vector<realnum>&,long[]);
STATIC void CoStarListModels(const stellar_grid*);
STATIC bool RauchInitialize(const char[],const char[],const vector<mpp>&,long,long,
				 long,const double[],int);
STATIC void RauchReadMPP(vector<mpp>&,vector<mpp>&,vector<mpp>&,vector<mpp>&,vector<mpp>&,vector<mpp>&);
inline void getdataline(fstream&,string&);
STATIC void WriteASCIIHead(FILE*,long,long,const vector<string>&,long,long,const string&,double,
			   const string&,double,const vector<mpp>&,const char*,int);
STATIC void WriteASCIIData(FILE*,const vector<double>&,long,const char*,int);
STATIC bool lgReadAtmosphereHead(stellar_grid*);
STATIC bool lgReadAtmosphereTail(stellar_grid*,const realnum[],long,const vector<long>&);
STATIC bool lgCompileAtmosphere(const char[],const char[],const realnum[],long,process_counter&);
STATIC void InitGrid(stellar_grid*,bool,bool=lgREAD_BIN);
inline void InitGridASCII(stellar_grid*);
STATIC void InitGridBin(stellar_grid*);
STATIC bool lgValidBinFile(const char*,process_counter&,access_scheme);
STATIC bool lgValidASCIIFile(const char*,access_scheme);
STATIC void InitGridCoStar(stellar_grid*);
STATIC void CheckVal(const stellar_grid*,double[],long*,long*);
STATIC void InterpolateRectGrid(stellar_grid*,const double[],double*,double*,bool=lgTAKELOG,const realnum[]=NULL,
				long=0L);
STATIC void InterpolateModel(stellar_grid*,const double[],vector<double>&,const vector<long>&,const vector<long>&,
			     vector<long>&,long,vector<realnum>&,bool,const realnum[]=NULL,long=0L);
STATIC void InterpolateModel(stellar_grid*,const double[],vector<double>&,const vector<long>&,
			     const vector<long>&,vector<long>&,long,vector<realnum>&,int);
STATIC void InterpolateModelCoStar(const stellar_grid*,const double[],vector<double>&,const long[],
				   const long[],vector<long>&,long,long,vector<realnum>&);
template<class T> void SortUnique(vector<T>&,vector<T>&);
STATIC void GetBins(const stellar_grid*,vector<Energy>&);
STATIC void GetModel(const stellar_grid*,long,realnum*,bool,bool);
STATIC void SetLimits(const stellar_grid*,double,const vector<long>&,const vector<long>&,const long[],
		      const vector<realnum>&,double*,double*);
STATIC void SetLimitsSub(const stellar_grid*,double,const vector<long>&,const vector<long>&,vector<long>&,
			 long,double*,double*);
STATIC void InitIndexArrays(stellar_grid*,bool);
STATIC void FillJ(stellar_grid*,vector<long>&,vector<double>&,long,bool);
STATIC long JIndex(const stellar_grid*,const vector<long>&);
STATIC void SearchModel(const vector<mpp>&,bool,long,const vector<double>&,long,long*,long*);
STATIC void FindIndex(const multi_arr<double,2>&,long,long,double,long*,long*,bool*);
STATIC bool lgFileReadable(const char*, process_counter&,access_scheme);
STATIC void ValidateMesh(const stellar_grid*,const vector<Energy>&);
STATIC bool lgValidMesh(const vector<Energy>&);			   
STATIC void ValidateGrid(const stellar_grid*,double);
STATIC bool lgValidModel(const vector<Energy>&,const vector<realnum>&,double,double);
STATIC void RebinAtmosphere(const vector<realnum>&,const vector<realnum>&,long,const realnum[],long,realnum[]);
STATIC realnum RebinSingleCell(long,const realnum[],const realnum[],const vector<realnum>&,long);
inline long RebinFind(const realnum[],long,realnum);
template<class T> void DumpAtmosphere(const char *fnam,long,long,char[MDIM][MNAM+1],const vector<mpp>&,
				      long,const T[],const realnum[]);


/* the version number for the ascii/binary atmosphere files */
static const long int VERSION_ASCII = 20060612L;
static const long int VERSION_COSTAR = 20160614L; // special ascii format for the CoStar grids
/* binary files are incompatible when floats are converted to doubles */
#ifdef FLT_IS_DBL
static const long int VERSION_BIN = 201009020L;
#else
static const long int VERSION_BIN = 201009021L;
#endif
static const long int VERSION_RAUCH_MPP = 20090324;

/* define the major absorption edges that require special attention during rebinning
 *
 * NB the frequencies should be chosen here such that they are in the middle of
 * the two frequency points that straddle the edge in the atmosphere model. The
 * software in RebinAtmosphere will seek out the exact values of those two points
 * e.g.: in the CoStar models the H I edge is straddled by wavelength points at
 * 911.67 and 911.85 A, so Edges[0] should be chosen at 911.76A.
 *
 * NB beware not to choose edges too close to one another (i.e. on the order of the
 * resolution of the Cloudy frequency grid). E.g. the He II Balmer edge nearly coincides
 * with the H I Ly edge, they should be treated as one edge. Trying to separate them will
 * almost certainly lead to problems in RebinAtmosphere */
static const realnum Edges_CoStar[] = { realnum(RYDLAM/911.76), realnum(RYDLAM/504.26), realnum(RYDLAM/227.84) };

/** List all the available TABLE STAR <grid> commands by checking installed *.mod files */
void AtmospheresAvail()
{
	DEBUG_ENTRY( "AtmospheresAvail()" );

	/* This routine makes a list of all the stellar atmosphere grids that are valid,
	 * giving the parameters for use in the input script as well. It is simply a long
	 * list of if-statements, so if any grid is added to Cloudy, it should be added in
	 * this routine as well.
	 *
	 * NB NB NB -- test this routine regularly to see if the list is still complete! */

	fprintf( ioQQQ, "\n I will now list all stellar atmosphere grids that are ready to be used (if any).\n" );
	fprintf( ioQQQ, " User-defined stellar atmosphere grids will not be included in this list.\n\n" );

	process_counter dum;

	/* we always look in the data directory regardless of where we are,
	 * it would be very confusing to the user if we did otherwise... */
	access_scheme as = AS_DATA_ONLY_TRY;

	if( lgValidBinFile( "atlas_fp10k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp05k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp03k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp02k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.2 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp01k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.1 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp00k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm01k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-0.1 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm02k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-0.2 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm03k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm05k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-0.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm10k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm15k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-1.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm20k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-2.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm25k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-2.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm30k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-3.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm35k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-3.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm40k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-4.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm45k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-4.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm50k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-5.0 <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "atlas_fp05k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z+0.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp02k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z+0.2 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp00k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm05k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-0.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm10k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm15k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-1.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm20k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-2.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm25k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-2.5 <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "atlas_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "atlas_3d_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "Sc1_costar_solar.mod", dum, as ) )
		fprintf( ioQQQ, "   table star costar solar (see Hazy for parameters)\n" );
	if( lgValidBinFile( "Sc1_costar_halo.mod", dum, as ) )
		fprintf( ioQQQ, "   table star costar halo (see Hazy for parameters)\n" );

	if( lgValidBinFile( "kurucz79.mod", dum, as ) )
		fprintf( ioQQQ, "   table star kurucz79 <Teff>\n" );

	if( lgValidBinFile( "mihalas.mod", dum, as ) )
		fprintf( ioQQQ, "   table star mihalas <Teff>\n" );

	if( lgValidBinFile( "rauch_h-ca_solar.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ca solar <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_h-ca_halo.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ca halo <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_h-ca_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ca 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "rauch_h-ni_solar.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ni solar <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_h-ni_halo.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ni halo <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_h-ni_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ni 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "rauch_pg1159.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch pg1159 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_cowd.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch co wd <Teff>\n" );

	if( lgValidBinFile( "rauch_hydr.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch hydrogen <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "rauch_helium.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch helium <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "rauch_h+he_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H+He <Teff> <log(g)> <frac(He)>\n" );

	if( lgValidBinFile( "starburst99.mod", dum, as ) )
		fprintf( ioQQQ, "   table star \"starburst99.mod\" <age>\n" );
	if( lgValidBinFile( "starburst99_2d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star \"starburst99_2d.mod\" <age> <Z>\n" );

	if( lgValidBinFile( "obstar_merged_p03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z+0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_p00.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_m03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z-0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_m07.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z-0.7 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_m10.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_m99.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z-inf <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "obstar_merged_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "bstar2006_p03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z+0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_p00.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_m03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z-0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_m07.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z-0.7 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_m10.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_m99.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z-inf <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "bstar2006_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "ostar2002_p03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z+0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_p00.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m07.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-0.7 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m10.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m15.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-1.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m17.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-1.7 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m20.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-2.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m30.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-3.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m99.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-inf <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "ostar2002_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "kwerner.mod", dum, as ) )
		fprintf( ioQQQ, "   table star werner <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "wmbasic.mod", dum, as ) )
		fprintf( ioQQQ, "   table star wmbasic <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidASCIIFile( "hm05_galaxy.ascii", as ) )
		fprintf( ioQQQ, "   table HM05 <z> [ <factor> ]\n" );
	if( lgValidASCIIFile( "hm05_quasar.ascii", as ) )
		fprintf( ioQQQ, "   table HM05 quasar <z> [ <factor> ]\n" );
	if( lgValidASCIIFile( "hm12_galaxy.ascii", as ) )
		fprintf( ioQQQ, "   table HM12 <z> [ <factor> ]\n" );
}

/* AtlasCompile rebin Kurucz stellar models to match energy grid of code */
/* >>chng 05 nov 16, added return value to indicate success (0) or failure (1) */
int AtlasCompile(process_counter& pc)
{
	DEBUG_ENTRY( "AtlasCompile()" );

	/* This is a program to re-bin the Kurucz stellar models spectrum to match the 
	 * CLOUDY grid.  For wavelengths shorter than supplied in the Kurucz files,
	 * the flux will be set to zero.  At long wavelengths a Rayleigh-Jeans
	 * extrapolation will be used. */

	/* This version uses power-law interpolation between the points of the stellar
	 * model.*/

	fprintf( ioQQQ, " AtlasCompile on the job.\n" );

	access_scheme as = AS_LOCAL_ONLY_TRY;

	/* >>chng 05 nov 19, add support for non-solar metalicities as well as odfnew models, PvH */
	bool lgFail = false;
	if( lgFileReadable( "atlas_fp10k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp10k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp10k2.ascii", "atlas_fp10k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fp05k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp05k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp05k2.ascii", "atlas_fp05k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fp03k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp03k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp03k2.ascii", "atlas_fp03k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fp02k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp02k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp02k2.ascii", "atlas_fp02k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fp01k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp01k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp01k2.ascii", "atlas_fp01k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fp00k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp00k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp00k2.ascii", "atlas_fp00k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm01k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm01k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm01k2.ascii", "atlas_fm01k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm02k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm02k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm02k2.ascii", "atlas_fm02k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm03k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm03k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm03k2.ascii", "atlas_fm03k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm05k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm05k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm05k2.ascii", "atlas_fm05k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm10k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm10k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm10k2.ascii", "atlas_fm10k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm15k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm15k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm15k2.ascii", "atlas_fm15k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm20k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm20k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm20k2.ascii", "atlas_fm20k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm25k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm25k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm25k2.ascii", "atlas_fm25k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm30k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm30k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm30k2.ascii", "atlas_fm30k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm35k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm35k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm35k2.ascii", "atlas_fm35k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm40k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm40k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm40k2.ascii", "atlas_fm40k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm45k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm45k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm45k2.ascii", "atlas_fm45k2.mod", NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm50k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm50k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm50k2.ascii", "atlas_fm50k2.mod", NULL, 0L, pc );

	if( lgFileReadable( "atlas_fp05k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fp05k2_odfnew.mod", pc, as ) )
	    
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp05k2_odfnew.ascii", "atlas_fp05k2_odfnew.mod",
							NULL, 0L, pc );
	if( lgFileReadable( "atlas_fp02k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fp02k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp02k2_odfnew.ascii", "atlas_fp02k2_odfnew.mod",
							NULL, 0L, pc );
	if( lgFileReadable( "atlas_fp00k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fp00k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp00k2_odfnew.ascii", "atlas_fp00k2_odfnew.mod",
							NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm05k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm05k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm05k2_odfnew.ascii", "atlas_fm05k2_odfnew.mod",
							NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm10k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm10k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm10k2_odfnew.ascii", "atlas_fm10k2_odfnew.mod",
							NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm15k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm15k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm15k2_odfnew.ascii", "atlas_fm15k2_odfnew.mod",
							NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm20k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm20k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm20k2_odfnew.ascii", "atlas_fm20k2_odfnew.mod",
							NULL, 0L, pc );
	if( lgFileReadable( "atlas_fm25k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm25k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm25k2_odfnew.ascii", "atlas_fm25k2_odfnew.mod",
							NULL, 0L, pc );

	if( lgFileReadable( "atlas_3d.ascii", pc, as ) && !lgValidBinFile( "atlas_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_3d.ascii", "atlas_3d.mod", NULL, 0L, pc );

	if( lgFileReadable( "atlas_3d_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_3d_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_3d_odfnew.ascii", "atlas_3d_odfnew.mod", NULL, 0L, pc );
	return lgFail;
}

/* AtlasInterpolate read in and interpolate on Kurucz grid of atmospheres, originally by K Volk */
long AtlasInterpolate(double val[], /* val[nval] */
		      long *nval,
		      long *ndim,
		      const char *chMetalicity,
		      const char *chODFNew,
		      bool lgList,
		      double *Tlow,
		      double *Thigh)
{
	DEBUG_ENTRY( "AtlasInterpolate()" );

	stellar_grid grid;
	grid.name = "atlas_";
	if( *ndim == 3 )
		grid.name += "3d";
	else
	{
		grid.name += "f";
		grid.name += chMetalicity;
		grid.name += "k2";
	}
	grid.name += chODFNew;
	grid.name += ".mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	char chIdent[13];
	if( *ndim == 3 )
	{
		strcpy( chIdent, "3-dim" );
	}
	else
	{
		strcpy( chIdent, "Z " );
		strcat( chIdent, chMetalicity );
	}
	strcat( chIdent, ( strlen(chODFNew) == 0 ? " Kurucz" : " ODFNew" ) );
	grid.ident = chIdent;
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	/* Note on the interpolation (solar abundance grid): 26 October 2000 (Peter van Hoof)
	 *
	 * I computed the effective temperature for a random sample of interpolated
	 * atmospheres by integrating the flux as shown above and compared the results
	 * with the expected effective temperature using DELTA = (COMP-EXPEC)/EXPEC.
	 *
	 * I found that the average discrepancy was:
	 *
	 *     DELTA = -0.10% +/- 0.06% (sample size 5000)
	 *
	 * The most extreme discrepancies were
	 *     -0.30% <= DELTA <= 0.21%
	 *
	 * The most negative discrepancies were for Teff =  36 -  39 kK, log(g) = 4.5 - 5
	 * The most positive discrepancies were for Teff = 3.5 - 4.0 kK, log(g) = 0 - 1
	 *
	 * The interpolation in the ATLAS grid is clearly very accurate */

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* CoStarCompile rebin costar stellar models to match energy grid of code*/
int CoStarCompile(process_counter& pc)
{
	DEBUG_ENTRY( "CoStarCompile()" );

	fprintf( ioQQQ, " CoStarCompile on the job.\n" );

	access_scheme as = AS_LOCAL_ONLY_TRY;

	bool lgFail = false;
	if( lgFileReadable( "Sc1_costar_z020_lb.fluxes", pc, as ) && !lgValidASCIIFile( "Sc1_costar_solar.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating Sc1_costar_solar.ascii....\n" );
		lgFail = lgFail || CoStarInitialize( "Sc1_costar_z020_lb.fluxes", "Sc1_costar_solar.ascii" );
	}
	if( lgFileReadable( "Sc1_costar_z004_lb.fluxes", pc, as ) && !lgValidASCIIFile( "Sc1_costar_halo.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating Sc1_costar_halo.ascii....\n" );
		lgFail = lgFail || CoStarInitialize( "Sc1_costar_z004_lb.fluxes", "Sc1_costar_halo.ascii" );
	}

	if( lgFileReadable( "Sc1_costar_solar.ascii", pc, as ) && !lgValidBinFile( "Sc1_costar_solar.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "Sc1_costar_solar.ascii", "Sc1_costar_solar.mod",
							Edges_CoStar, 3L, pc );
	if( lgFileReadable( "Sc1_costar_halo.ascii", pc, as ) && !lgValidBinFile( "Sc1_costar_halo.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "Sc1_costar_halo.ascii", "Sc1_costar_halo.mod",
							Edges_CoStar, 3L, pc );

	return lgFail;
}

/* CoStarInterpolate read in and interpolate on CoStar grid of atmospheres */
long CoStarInterpolate(double val[], /* requested model parameters */
		       long *nval,
		       long *ndim,
		       IntMode imode, /* which interpolation mode is requested */
		       bool lgHalo,  /* flag indicating whether solar (==0) or halo (==1) abundances */
		       bool lgList,
		       double *val0_lo,
		       double *val0_hi)
{
	DEBUG_ENTRY( "CoStarInterpolate()" );

	stellar_grid grid;
	grid.name = ( lgHalo ? "Sc1_costar_halo.mod" : "Sc1_costar_solar.mod" );
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "      costar";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	/* listing the models in the grid is implemented in CoStarListModels() */
	InitGrid( &grid, false );
	/* now sort the models according to track */
	InitGridCoStar( &grid );
	/* override default interpolation mode */
	grid.imode = imode;

	if( lgList )
	{
		CoStarListModels( &grid );
		cdEXIT(EXIT_SUCCESS);
	}

	CheckVal( &grid, val, nval, ndim );

	/* Note on the interpolation: 26 October 2000 (Peter van Hoof)
	 *
	 * I computed the effective temperature for a random sample of interpolated
	 * atmospheres by integrating the flux as shown above and compared the results
	 * with the expected effective temperature using DELTA = (COMP-EXPEC)/EXPEC.
	 *
	 * I found that the average discrepancy was:
	 *
	 *     DELTA = -1.16% +/- 0.69% (SOLAR models, sample size 4590)
	 *     DELTA = -1.17% +/- 0.70% (HALO models, sample size 4828)
	 *
	 * The most extreme discrepancies for the SOLAR models were
	 *     -3.18% <= DELTA <= -0.16%
	 *
	 * The most negative discrepancies were for  Teff = 35 kK, log(g) = 3.5
	 * The least negative discrepancies were for Teff = 50 kK, log(g) = 4.1
	 *
	 * The most extreme discrepancies for the HALO models were
	 *     -2.90% <= DELTA <= -0.13%
	 *
	 * The most negative discrepancies were for  Teff = 35 kK, log(g) = 3.5
	 * The least negative discrepancies were for Teff = 50 kK, log(g) = 4.1
	 *
	 * Since Cloudy checks the scaling elsewhere there is no need to re-scale 
	 * things here, but this inaccuracy should be kept in mind since it could
	 * indicate problems with the flux distribution */

	InterpolateGridCoStar( &grid, val, val0_lo, val0_hi );

	return rfield.nflux_with_check;
}

/* GridCompile rebin user supplied stellar models to match energy grid of code */
bool GridCompile(const char *InName)
{
	DEBUG_ENTRY( "GridCompile()" );

	fprintf( ioQQQ, " GridCompile on the job.\n" );

	// replace filename extension with ".mod"
	string OutName( InName );
	string::size_type ptr = OutName.rfind( '.' );
	ASSERT( ptr != string::npos );
	OutName.replace( ptr, string::npos, ".mod" );

	process_counter dum;
	bool lgFail = lgCompileAtmosphere( InName, OutName.c_str(), NULL, 0L, dum );

	if( !lgFail )
	{
		stellar_grid grid;

		/* the file must be local */
		grid.name = OutName;
		grid.scheme = AS_LOCAL_ONLY;
		grid.ident = "bogus ident.";
		grid.command = "bogus command.";

		InitGrid( &grid, false );

		/* check whether the models in the grid have the correct effective temperature */

		if( strcmp( grid.names[0], "Teff" ) == 0 )
		{
			fprintf( ioQQQ, " GridCompile: checking effective temperatures...\n" );
			ValidateGrid( &grid, 0.02 );
		}
	}
	return lgFail;
}

/* GridInterpolate read in and interpolate on user supplied grid of atmospheres */
long GridInterpolate(double val[], /* val[nval] */
		     long *nval,
		     long *ndim,
		     const char *FileName,
		     bool lgList,
		     double *Tlow,
		     double *Thigh)
{
	DEBUG_ENTRY( "GridInterpolate()" );

	// make filename without extension
	string chTruncName( FileName );
	string::size_type ptr = chTruncName.rfind( '.' );
	if( ptr != string::npos )
		chTruncName.replace( ptr, string::npos, "" );

	stellar_grid grid;
	grid.name = FileName;
	grid.scheme = AS_DATA_OPTIONAL;
	bool lgASCII = lgValidASCIIFile( grid.name.c_str(), AS_DATA_ONLY_TRY );
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = chTruncName.substr(0,12);
	if( grid.ident.length() < 12 )
		grid.ident.append(12-grid.ident.length(),' ');
	/* the Cloudy command needed to recompile the binary model file */
	if( !lgASCII )
		grid.command = "COMPILE STARS \"" + chTruncName + ".ascii\"";

	InitGrid( &grid, lgList, lgASCII );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* HaardtMadauInterpolate read in and interpolate on Haardt & Madau SEDs */
long HaardtMadauInterpolate(double val,
			    int version,
			    bool lgQuasar,
			    double *zlow,
			    double *zhigh)
{
	DEBUG_ENTRY( "HaardtMadauInterpolate()" );

	stellar_grid grid;
	string name, ident;
	if( version == 2005 )
	{
		grid.name = "hm05_";
		/* identification of this atmosphere set, used in
		 * the Cloudy output, *must* be 12 characters long */
		grid.ident = " HM05 ";
	}
	else if( version == 2012 )
	{
		grid.name = "hm12_";
		grid.ident = " HM12 ";
	}
	else
		TotalInsanity();
	if( lgQuasar )
	{
		grid.name += "quasar.ascii";
		grid.ident += "QUASAR";
	}
	else
	{
		grid.name += "galaxy.ascii";
		grid.ident += "GALAXY";
	}
	grid.scheme = AS_DEFAULT;

	/* define the major absorption edges that require special attention during rebinning
	 * see the routine CoStarCompile() for a more detailed discussion of this array */
	/* the CUBA code has unusual edges coinciding with the H I and He II Lyman lines. See
	 * section 2.2 of Haardt & Madau (2012) for a detailed discussion. We omit the edges
	 * from the highest Lyman lines since they are spaced too closely (see CoStarCompile).
	 * Omitted are H I 930.7, 926.2, 923.1, 920.9 and He II 232.7, 231.5, 230.8, 230.2 */
	vector<realnum> Edges;
	Edges.push_back( realnum(RYDLAM/1216.0) );
	if( version == 2012 )
	{
		Edges.push_back( realnum(RYDLAM/1026.0) );
		Edges.push_back( realnum(RYDLAM/972.5) );
		Edges.push_back( realnum(RYDLAM/949.7) );
		Edges.push_back( realnum(RYDLAM/937.8) );
	}

	Edges.push_back( realnum(RYDLAM/304.0) );
	if( version == 2012 )
	{
		Edges.push_back( realnum(RYDLAM/256.4) );
		Edges.push_back( realnum(RYDLAM/243.1) );
		Edges.push_back( realnum(RYDLAM/237.4) );
		Edges.push_back( realnum(RYDLAM/234.4) );		
	}
	Edges.push_back( realnum(RYDLAM/228.0) );

	InitGrid( &grid, false, lgREAD_ASCII );

	long nval = 1, ndim = 1;
	CheckVal( &grid, &val, &nval, &ndim );

	InterpolateRectGrid( &grid, &val, zlow, zhigh, lgLINEAR, get_ptr(Edges), Edges.size() );

	return rfield.nflux_with_check;
}

/* KhaireSrianandInterpolate read in and interpolate on Khaire & Srianand SEDs */
long KhaireSrianandInterpolate(double val,
			       int Q,
			       double *zlow,
			       double *zhigh)
{
	DEBUG_ENTRY( "KhaireSrianandInterpolate()" );

	stellar_grid grid;
	ostringstream oss;
	oss << "ks18_q" << Q << ".ascii";
	grid.name = oss.str();
	ostringstream oss2;
	oss2 << " KS18, Q=" << Q << " ";
	grid.ident = oss2.str();
	/* these files are part of the distribution in the data directory */
	grid.scheme = AS_DEFAULT;

	/* there are edges at 1215.67, 911.78, 503.98, 303.78, and 227.84 A
	 * after discussion with Vikram Khaire it was decided not to protect these */

	InitGrid( &grid, false, lgREAD_ASCII );

	long nval = 1, ndim = 1;
	CheckVal( &grid, &val, &nval, &ndim );

	InterpolateRectGrid( &grid, &val, zlow, zhigh, lgLINEAR, NULL, 0L );

	return rfield.nflux_with_check;
}

/* Kurucz79Compile rebin Kurucz 1979 stellar models to match energy grid of code */
int Kurucz79Compile(process_counter& pc)
{
	DEBUG_ENTRY( "Kurucz79Compile()" );

	fprintf( ioQQQ, " Kurucz79Compile on the job.\n" );

	/* following atmospheres LTE from Kurucz 1979, Ap.J. Sup 40, 1. and
	 * Kurucz (1989) private communication, newer opacities */

	access_scheme as = AS_LOCAL_ONLY_TRY;

	bool lgFail = false;
	if( lgFileReadable( "kurucz79.ascii", pc, as ) && !lgValidBinFile( "kurucz79.mod", pc, as ) )
		lgFail = lgCompileAtmosphere( "kurucz79.ascii", "kurucz79.mod", NULL, 0L, pc );
	return lgFail;
}

/* Kurucz79Interpolate read in and interpolate on Kurucz79 grid of atmospheres */
long Kurucz79Interpolate(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgList,
			 double *Tlow,
			 double *Thigh)
{
	DEBUG_ENTRY( "Kurucz79Interpolate()" );

	stellar_grid grid;
	grid.name = "kurucz79.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = " Kurucz 1979";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* MihalasCompile rebin Mihalas stellar models to match energy grid of code */
int MihalasCompile(process_counter& pc)
{
	DEBUG_ENTRY( "MihalasCompile()" );

	fprintf( ioQQQ, " MihalasCompile on the job.\n" );

	/* following atmospheres NLTE from Mihalas, NCAR-TN/STR-76 */

	/* define the major absorption edges that require special attention during rebinning
	 * see the routine CoStarCompile() for a more detailed discussion of this array */
	realnum Edges[3];
	Edges[0] = (realnum)(RYDLAM/911.204);
	Edges[1] = (realnum)(RYDLAM/503.968);
	Edges[2] = (realnum)(RYDLAM/227.806);

	access_scheme as = AS_LOCAL_ONLY_TRY;

	bool lgFail = false;
	if( lgFileReadable( "mihalas.ascii", pc, as ) && !lgValidBinFile( "mihalas.mod", pc, as ) )
		lgFail = lgCompileAtmosphere( "mihalas.ascii", "mihalas.mod", Edges, 3L, pc );
	return lgFail;
}

/* MihalasInterpolate read in and interpolate on Mihalas grid of atmospheres */
long MihalasInterpolate(double val[], /* val[nval] */
			long *nval,
			long *ndim,
			bool lgList,
			double *Tlow,
			double *Thigh)
{
	DEBUG_ENTRY( "MihalasInterpolate()" );

	stellar_grid grid;
	grid.name = "mihalas.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "     Mihalas";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* RauchCompile create ascii and mod files for Rauch atmospheres */
int RauchCompile(process_counter& pc)
{
	DEBUG_ENTRY( "RauchCompile()" );

	/* metalicities of the solar and halo grid */
	static const double par2[2] = { 0., -1. };

	/* Helium fraction by mass */
	static const double par3[11] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

	/* Before running this program issue the following command where the Rauch
	 * model atmosphere files are kept (0050000_50_solar_bin_0.1 and so on)
	 *
	 *   ls *solar_bin_0.1 > rauchmods.list
	 *
	 * and check to see that there are 66 lines in the file.
	 */

	fprintf( ioQQQ, " RauchCompile on the job.\n" );

	vector<mpp> telg1(NMODS_HCA);
	vector<mpp> telg2(NMODS_HNI);
	vector<mpp> telg3(NMODS_PG1159);
	vector<mpp> telg4(NMODS_HYDR);
	vector<mpp> telg5(NMODS_HELIUM);
	vector<mpp> telg6(NMODS_HpHE);

	RauchReadMPP( telg1, telg2, telg3, telg4, telg5, telg6 );

	process_counter dum;
	access_scheme as = AS_LOCAL_ONLY_TRY;
	bool lgFail = false;

	/* this is the H-Ca grid */
	if( lgFileReadable( "0050000_50_solar_bin_0.1", dum, as ) && !lgValidASCIIFile( "rauch_h-ca_solar.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ca_solar.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_h-ca_solar.ascii", "_solar_bin_0.1", 
						    telg1, NMODS_HCA, 1, 1, par2, 1 );
	}

	if( lgFileReadable( "0050000_50_halo__bin_0.1", dum, as ) && !lgValidASCIIFile( "rauch_h-ca_halo.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ca_halo.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_h-ca_halo.ascii", "_halo__bin_0.1",
						    telg1, NMODS_HCA, 1, 1, par2, 1 );
	}

	if( lgFileReadable( "0050000_50_solar_bin_0.1", dum, as ) &&
	    lgFileReadable( "0050000_50_halo__bin_0.1", dum, as ) &&
	    !lgValidASCIIFile( "rauch_h-ca_3d.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ca_3d.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_h-ca_3d.ascii", "_solar_bin_0.1",
						    telg1, NMODS_HCA, 1, 2, par2, 1 );
		lgFail = lgFail || RauchInitialize( "rauch_h-ca_3d.ascii", "_halo__bin_0.1",
						    telg1, NMODS_HCA, 2, 2, par2, 1 );
	}

	/* this is the H-Ni grid */
	if( lgFileReadable( "0050000_50_solar_iron.bin_0.1", dum, as ) &&
	    !lgValidASCIIFile( "rauch_h-ni_solar.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ni_solar.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_h-ni_solar.ascii", "_solar_iron.bin_0.1",
						    telg2, NMODS_HNI, 1, 1, par2, 1 );
	}

	if( lgFileReadable( "0050000_50_halo__iron.bin_0.1", dum, as ) &&
	    !lgValidASCIIFile( "rauch_h-ni_halo.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ni_halo.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_h-ni_halo.ascii", "_halo__iron.bin_0.1",
						    telg2, NMODS_HNI, 1, 1, par2, 1 );
	}

	if( lgFileReadable( "0050000_50_solar_iron.bin_0.1", dum, as ) &&
	    lgFileReadable( "0050000_50_halo__iron.bin_0.1", dum, as ) &&
	    !lgValidASCIIFile( "rauch_h-ni_3d.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ni_3d.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_h-ni_3d.ascii", "_solar_iron.bin_0.1",
						    telg2, NMODS_HNI, 1, 2, par2, 1 );
		lgFail = lgFail || RauchInitialize( "rauch_h-ni_3d.ascii", "_halo__iron.bin_0.1",
						    telg2, NMODS_HNI, 2, 2, par2, 1 );
	}

	/* this is the hydrogen deficient PG1159 grid */
	if( lgFileReadable( "0040000_5.00_33_50_02_15.bin_0.1", dum, as ) &&
	    !lgValidASCIIFile( "rauch_pg1159.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_pg1159.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_pg1159.ascii", "_33_50_02_15.bin_0.1",
						    telg3, NMODS_PG1159, 1, 1, par2, 2 );
	}

	/* this is the pure hydrogen grid */
	if( lgFileReadable( "0020000_4.00_H_00005-02000A.bin_0.1", dum, as ) &&
	    !lgValidASCIIFile( "rauch_hydr.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_hydr.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_hydr.ascii", "_H_00005-02000A.bin_0.1",
						    telg4, NMODS_HYDR, 1, 1, par2, 2 );
	}

	/* this is the pure helium grid */
	if( lgFileReadable( "0050000_5.00_He_00005-02000A.bin_0.1", dum, as ) &&
	    !lgValidASCIIFile( "rauch_helium.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_helium.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_helium.ascii", "_He_00005-02000A.bin_0.1",
						    telg5, NMODS_HELIUM, 1, 1, par2, 2 );
	}

	/* this is the 3D grid for arbitrary H+He mixtures */
	if( lgFileReadable( "0050000_5.00_H+He_1.000_0.000_00005-02000A.bin_0.1", dum, as ) &&
	    !lgValidASCIIFile( "rauch_h+he_3d.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h+he_3d.ascii....\n" );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_1.000_0.000_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  1, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.900_0.100_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  2, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.800_0.200_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  3, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.700_0.300_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  4, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.600_0.400_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  5, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.500_0.500_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  6, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.400_0.600_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  7, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.300_0.700_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  8, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.200_0.800_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE,  9, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.100_0.900_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE, 10, 11, par3, 2 );
		lgFail = lgFail || RauchInitialize( "rauch_h+he_3d.ascii", "_H+He_0.000_1.000_00005-02000A.bin_0.1",
						    telg6, NMODS_HpHE, 11, 11, par3, 2 );
	}

	if( lgFileReadable( "rauch_h-ca_solar.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ca_solar.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ca_solar.ascii", "rauch_h-ca_solar.mod", NULL,0L, pc );
	if( lgFileReadable( "rauch_h-ca_halo.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ca_halo.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ca_halo.ascii", "rauch_h-ca_halo.mod", NULL, 0L, pc );
	if( lgFileReadable( "rauch_h-ca_3d.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ca_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ca_3d.ascii", "rauch_h-ca_3d.mod", NULL, 0L, pc );

	if( lgFileReadable( "rauch_h-ni_solar.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ni_solar.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ni_solar.ascii", "rauch_h-ni_solar.mod", NULL,0L, pc );
	if( lgFileReadable( "rauch_h-ni_halo.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ni_halo.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ni_halo.ascii", "rauch_h-ni_halo.mod", NULL, 0L, pc );
	if( lgFileReadable( "rauch_h-ni_3d.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ni_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ni_3d.ascii", "rauch_h-ni_3d.mod", NULL, 0L, pc );

	if( lgFileReadable( "rauch_pg1159.ascii", pc, as ) && !lgValidBinFile( "rauch_pg1159.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_pg1159.ascii", "rauch_pg1159.mod", NULL, 0L, pc );
	if( lgFileReadable( "rauch_cowd.ascii", pc, as ) && !lgValidBinFile( "rauch_cowd.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_cowd.ascii", "rauch_cowd.mod", NULL, 0L, pc );

	if( lgFileReadable( "rauch_hydr.ascii", pc, as ) && !lgValidBinFile( "rauch_hydr.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_hydr.ascii", "rauch_hydr.mod", NULL, 0L, pc );

	if( lgFileReadable( "rauch_helium.ascii", pc, as ) && !lgValidBinFile( "rauch_helium.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_helium.ascii", "rauch_helium.mod", NULL, 0L, pc );

	if( lgFileReadable( "rauch_h+he_3d.ascii", pc, as ) && !lgValidBinFile( "rauch_h+he_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h+he_3d.ascii", "rauch_h+he_3d.mod", NULL, 0L, pc );
	return lgFail;
}

/* RauchInterpolateHCa get one of the Rauch H-Ca model atmospheres, originally by K. Volk */
long RauchInterpolateHCa(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgHalo,
			 bool lgList,
			 double *Tlow,
			 double *Thigh)
{
	DEBUG_ENTRY( "RauchInterpolateHCa()" );

	stellar_grid grid;
	if( *ndim == 3 )
		grid.name = "rauch_h-ca_3d.mod";
	else
		grid.name = ( lgHalo ? "rauch_h-ca_halo.mod" : "rauch_h-ca_solar.mod" );
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "  H-Ca Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* RauchInterpolateHNi get one of the Rauch H-Ni model atmospheres */
long RauchInterpolateHNi(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgHalo,
			 bool lgList,
			 double *Tlow,
			 double *Thigh)
{
	DEBUG_ENTRY( "RauchInterpolateHNi()" );

	stellar_grid grid;
	if( *ndim == 3 )
		grid.name = "rauch_h-ni_3d.mod";
	else
		grid.name = ( lgHalo ? "rauch_h-ni_halo.mod" : "rauch_h-ni_solar.mod" );
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "  H-Ni Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* RauchInterpolatePG1159 get one of the Rauch PG1159 model atmospheres */
long RauchInterpolatePG1159(double val[], /* val[nval] */
			    long *nval,
			    long *ndim,
			    bool lgList,
			    double *Tlow,
			    double *Thigh)
{
	DEBUG_ENTRY( "RauchInterpolatePG1159()" );

	stellar_grid grid;
	grid.name = "rauch_pg1159.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "PG1159 Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* RauchInterpolateCOWD get one of the Rauch C/O white dwarf model atmospheres */
long RauchInterpolateCOWD(double val[], /* val[nval] */
			  long *nval,
			  long *ndim,
			  bool lgList,
			  double *Tlow,
			  double *Thigh)
{
	DEBUG_ENTRY( "RauchInterpolateCOWD()" );

	stellar_grid grid;
	grid.name = "rauch_cowd.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "C/O WD Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* RauchInterpolateHydr get one of the Rauch pure hydrogen model atmospheres */
long RauchInterpolateHydr(double val[], /* val[nval] */
			  long *nval,
			  long *ndim,
			  bool lgList,
			  double *Tlow,
			  double *Thigh)
{
	DEBUG_ENTRY( "RauchInterpolateHydr()" );

	stellar_grid grid;
	grid.name = "rauch_hydr.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "  Hydr Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* RauchInterpolateHelium get one of the Rauch pure helium model atmospheres */
long RauchInterpolateHelium(double val[], /* val[nval] */
			    long *nval,
			    long *ndim,
			    bool lgList,
			    double *Tlow,
			    double *Thigh)
{
	DEBUG_ENTRY( "RauchInterpolateHelium()" );

	stellar_grid grid;
	grid.name = "rauch_helium.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "Helium Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* RauchInterpolateHpHe get one of the Rauch hydrogen plus helium model atmospheres */
long RauchInterpolateHpHe(double val[], /* val[nval] */
			  long *nval,
			  long *ndim,
			  bool lgList,
			  double *Tlow,
			  double *Thigh)
{
	DEBUG_ENTRY( "RauchInterpolateHpHe()" );

	stellar_grid grid;
	grid.name = "rauch_h+he_3d.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "  H+He Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* StarburstInitialize does the actual work of preparing the ascii file */
bool StarburstInitialize(const char chInName[],
			 const char chOutName[],
			 sb_mode mode)
{
	DEBUG_ENTRY( "StarburstInitialize()" );

	/* grab some space for the wavelengths and fluxes */
	vector<mpp> telg(MNTS);
	vector<double> wavl, fluxes[MNTS];
	wavl.reserve(NSB99);

	FILE *ioIn = open_data( chInName, "r", AS_LOCAL_ONLY );

	double lwavl = 0.;
	long nmods = 0;
	long ngp = 0;

	bool lgHeader = true;
	char chLine[INPUT_LINE_LENGTH];
	while( read_whole_line( chLine, INPUT_LINE_LENGTH, ioIn ) != NULL )
	{
		if( !lgHeader )
		{
			/* format: age/yr wavl/Angstrom log10(flux_total) log10(flux_stellar) log10(flux_neb) */
			/* we are only interested in the total flux, so we ignore the remaining numbers */
			double cage, cwavl, cfl, cfl1, cfl2, cfl3;
			if( sscanf( chLine, " %le %le %le %le %le", &cage, &cwavl, &cfl1, &cfl2, &cfl3 ) != 5 )
			{
				fprintf( ioQQQ, "syntax error in data of Starburst grid.\n" );
				return true;
			}

			if( mode == SB_TOTAL )
				cfl = cfl1;
			else if( mode == SB_STELLAR )
				cfl = cfl2;
			else if( mode == SB_NEBULAR )
				cfl = cfl3;
			else
				TotalInsanity();

			if( cwavl < lwavl )
			{
				++nmods;
				ngp = 0;

				if( nmods >= MNTS )
				{
					fprintf( ioQQQ, "too many time steps in Starburst grid.\n" );
					fprintf( ioQQQ, "please increase MNTS and recompile.\n" );
					return true;
				}
			}

			if( ngp == 0 )
			{
				if( nmods == 0 )
					fluxes[nmods].reserve(NSB99);
				else
					fluxes[nmods].reserve(fluxes[nmods-1].size());
				telg[nmods].par[0] = cage;
			}

			if( !fp_equal(telg[nmods].par[0],cage,10) )
			{
				fprintf( ioQQQ, "age error in Starburst grid.\n" );
				return true;
			}

			if( nmods == 0 )
				wavl.push_back(cwavl);
			else
			{
				if( !fp_equal(wavl[ngp],cwavl,10) )
				{
					fprintf( ioQQQ, "wavelength error in Starburst grid.\n" );
					return true;
				}
			}

			/* arbitrarily renormalize to flux in erg/cm^2/s/A at 1kpc */
			/* constant is log10( 4*pi*(kpc/cm)^2 ) */
			fluxes[nmods].push_back( exp10(cfl - 44.077911) );

			lwavl = cwavl;
			++ngp;
		}

		if( lgHeader && strncmp( &chLine[1], "TIME [YR]", 9 ) == 0 )
			lgHeader = false;
	}

	if( lgHeader )
	{
		/* this happens when the "TIME [YR]" string was not found in column 1 of the file */
		fprintf( ioQQQ, "syntax error in header of Starburst grid.\n" );
		return true;
	}

	++nmods;

	/* finished - close the unit */
	fclose(ioIn);

	/* now write the ascii file */
	FILE *ioOut = open_data( chOutName, "w" );

	vector<string> names;
	names.push_back( "Age" );
	WriteASCIIHead(ioOut, VERSION_ASCII, 1, names, nmods, ngp, "lambda", 1., "F_lambda", 1., telg, " %.3e", 4);
	WriteASCIIData(ioOut, wavl, ngp, "  %.4e", 5);
	for( long i=0; i < nmods; i++ )
		WriteASCIIData(ioOut, fluxes[i], ngp, "  %.4e", 5);

	fclose(ioOut);
	return false;
}

/* StarburstCompile, rebin Starburst99 model output to match energy grid of code */
bool StarburstCompile(process_counter& pc)
{
	DEBUG_ENTRY( "StarburstCompile()" );

	fprintf( ioQQQ, " StarburstCompile on the job.\n" );

	process_counter dum;
	access_scheme as = AS_LOCAL_ONLY_TRY;

	bool lgFail = false;
	if( lgFileReadable( "starburst99.stb99", dum, as ) && !lgValidASCIIFile( "starburst99.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating starburst99.ascii....\n" );		
		lgFail = lgFail || StarburstInitialize( "starburst99.stb99", "starburst99.ascii", SB_TOTAL );
	}
	if( lgFileReadable( "starburst99.ascii", pc, as ) && !lgValidBinFile( "starburst99.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "starburst99.ascii", "starburst99.mod", NULL, 0L, pc );

	if( lgFileReadable( "starburst99_2d.ascii", pc, as ) && !lgValidBinFile( "starburst99_2d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "starburst99_2d.ascii", "starburst99_2d.mod", NULL, 0L, pc );
	return lgFail;
}

/* TlustyCompile rebin Tlusty BSTAR2006/OSTAR2002 stellar models to match energy grid of code */
int TlustyCompile(process_counter& pc)
{
	DEBUG_ENTRY( "TlustyCompile()" );

	fprintf( ioQQQ, " TlustyCompile on the job.\n" );

	access_scheme as = AS_LOCAL_ONLY_TRY;

	bool lgFail = false;
	if( lgFileReadable( "obstar_merged_p03.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_p03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere("obstar_merged_p03.ascii","obstar_merged_p03.mod",NULL,0L,pc);
	if( lgFileReadable( "obstar_merged_p00.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_p00.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere("obstar_merged_p00.ascii","obstar_merged_p00.mod",NULL,0L,pc);
	if( lgFileReadable( "obstar_merged_m03.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_m03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere("obstar_merged_m03.ascii","obstar_merged_m03.mod",NULL,0L,pc);
	if( lgFileReadable( "obstar_merged_m07.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_m07.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere("obstar_merged_m07.ascii","obstar_merged_m07.mod",NULL,0L,pc);
	if( lgFileReadable( "obstar_merged_m10.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_m10.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere("obstar_merged_m10.ascii","obstar_merged_m10.mod",NULL,0L,pc);
	if( lgFileReadable( "obstar_merged_m99.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_m99.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere("obstar_merged_m99.ascii","obstar_merged_m99.mod",NULL,0L,pc);

	if( lgFileReadable( "obstar_merged_3d.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere("obstar_merged_3d.ascii", "obstar_merged_3d.mod", NULL, 0L, pc);

	if( lgFileReadable( "bstar2006_p03.ascii", pc, as ) && !lgValidBinFile( "bstar2006_p03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_p03.ascii", "bstar2006_p03.mod", NULL, 0L, pc );
	if( lgFileReadable( "bstar2006_p00.ascii", pc, as ) && !lgValidBinFile( "bstar2006_p00.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_p00.ascii", "bstar2006_p00.mod", NULL, 0L, pc );
	if( lgFileReadable( "bstar2006_m03.ascii", pc, as ) && !lgValidBinFile( "bstar2006_m03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_m03.ascii", "bstar2006_m03.mod", NULL, 0L, pc );
	if( lgFileReadable( "bstar2006_m07.ascii", pc, as ) && !lgValidBinFile( "bstar2006_m07.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_m07.ascii", "bstar2006_m07.mod", NULL, 0L, pc );
	if( lgFileReadable( "bstar2006_m10.ascii", pc, as ) && !lgValidBinFile( "bstar2006_m10.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_m10.ascii", "bstar2006_m10.mod", NULL, 0L, pc );
	if( lgFileReadable( "bstar2006_m99.ascii", pc, as ) && !lgValidBinFile( "bstar2006_m99.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_m99.ascii", "bstar2006_m99.mod", NULL, 0L, pc );

	if( lgFileReadable( "bstar2006_3d.ascii", pc, as ) && !lgValidBinFile( "bstar2006_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_3d.ascii", "bstar2006_3d.mod", NULL, 0L, pc );

	if( lgFileReadable( "ostar2002_p03.ascii", pc, as ) && !lgValidBinFile( "ostar2002_p03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_p03.ascii", "ostar2002_p03.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_p00.ascii", pc, as ) && !lgValidBinFile( "ostar2002_p00.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_p00.ascii", "ostar2002_p00.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_m03.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m03.ascii", "ostar2002_m03.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_m07.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m07.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m07.ascii", "ostar2002_m07.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_m10.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m10.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m10.ascii", "ostar2002_m10.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_m15.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m15.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m15.ascii", "ostar2002_m15.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_m17.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m17.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m17.ascii", "ostar2002_m17.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_m20.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m20.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m20.ascii", "ostar2002_m20.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_m30.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m30.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m30.ascii", "ostar2002_m30.mod", NULL, 0L, pc );
	if( lgFileReadable( "ostar2002_m99.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m99.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m99.ascii", "ostar2002_m99.mod", NULL, 0L, pc );

	if( lgFileReadable( "ostar2002_3d.ascii", pc, as ) && !lgValidBinFile( "ostar2002_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_3d.ascii", "ostar2002_3d.mod", NULL, 0L, pc );
	return lgFail;
}

/* TlustyInterpolate get one of the Tlusty OBSTAR_MERGED/BSTAR2006/OSTAR2002 model atmospheres */
long TlustyInterpolate(double val[], /* val[nval] */
		       long *nval,
		       long *ndim,
		       tl_grid tlg,
		       const char *chMetalicity,
		       bool lgList,
		       double *Tlow,
		       double *Thigh)
{
	DEBUG_ENTRY( "TlustyInterpolate()" );

	stellar_grid grid;
	if( tlg == TL_OBSTAR )
		grid.name = "obstar_merged_";
	else if( tlg == TL_BSTAR )
		grid.name = "bstar2006_";
	else if( tlg == TL_OSTAR )
		grid.name = "ostar2002_";
	else
		TotalInsanity();
	if( *ndim == 3 )
		grid.name += "3d";
	else
		grid.name += chMetalicity;
	grid.name += ".mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	char chIdent[13];
	if( *ndim == 3 )
	{
		strcpy( chIdent, "3-dim" );
	}
	else
	{
		strcpy( chIdent, "Z " );
		strcat( chIdent, chMetalicity );
	}
	if( tlg == TL_OBSTAR )
		strcat( chIdent, " OBstar" );
	else if( tlg == TL_BSTAR )
		strcat( chIdent, " Bstr06" );
	else if( tlg == TL_OSTAR )
		strcat( chIdent, " Ostr02" );
	else
		TotalInsanity();
	grid.ident = chIdent;
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* WernerCompile rebin Werner stellar models to match energy grid of code */
/* >>chng 05 nov 16, added return value to indicate success (0) or failure (1) */
int WernerCompile(process_counter& pc)
{
	DEBUG_ENTRY( "WernerCompile()" );

	fprintf( ioQQQ, " WernerCompile on the job.\n" );

	/* define the major absorption edges that require special attention during rebinning
	 * see the routine CoStarCompile() for a more detailed discussion of this array */
	realnum Edges[3] = { 0.99946789f, 1.8071406f, 3.9996377f };

	/* The "kwerner.ascii" file is a modified ascii dump of the Klaus Werner 
	 * stellar model files which he gave to me in 1992.  The first set of values 
	 * is the frequency grid (in Ryd) followed by the atmosphere models in order
	 * of increasing temperature and log(g). The following comments are already
	 * incorporated in the modified kwerner.ascii file that is supplied with Cloudy.
	 *
	 * >>chng 00 oct 18, The frequency grid was slightly tweaked compared to the
	 * original values supplied by Klaus Werner to make it monotonically increasing;
	 * this is due to there being fluxes above and below certain wavelengths where
	 * the opacity changes (i.e. the Lyman and Balmer limits for example) which are 
	 * assigned the same wavelength in the original Klaus Werner files. PvH
	 *
	 * >>chng 00 oct 20, StarEner[172] is out of sequence. As per the Klaus Werner comment,
	 * it should be omitted. The energy grid is very dense in this region and was most likely
	 * intended to sample an absorption line which was not included in this particular grid.
	 * StarFlux[172] is therefore always equal to the flux in neighbouring points (at least
	 * those with slightly smaller energies). It is therefore safe to ignore this point. PvH
	 *
	 * >>chng 00 oct 20, As per the same comment, StarFlux[172] is also deleted. PvH */

	access_scheme as = AS_LOCAL_ONLY_TRY;

	bool lgFail = false;
	if( lgFileReadable( "kwerner.ascii", pc, as ) && !lgValidBinFile( "kwerner.mod", pc, as ) )
		lgFail = lgCompileAtmosphere( "kwerner.ascii", "kwerner.mod", Edges, 3L, pc );
	return lgFail;
}

/* WernerInterpolate read in and interpolate on Werner grid of PN atmospheres, originally by K Volk */
long WernerInterpolate(double val[], /* val[nval] */
		       long *nval,
		       long *ndim,
		       bool lgList,
		       double *Tlow,
		       double *Thigh)
{
	DEBUG_ENTRY( "WernerInterpolate()" );

	/* This subroutine was added (28 dec 1992) to read from the set of
	 * hot white dwarf model atmospheres from Klaus Werner at Kiel. The 
	 * values are read in (energy in Rydberg units, f_nu in cgs units)
	 * for any of the 20 models. Each model had 513 points before rebinning.
	 * The Rayleigh-Jeans tail was extrapolated. */

	stellar_grid grid;
	grid.name = "kwerner.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "Klaus Werner";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	/* Note on the interpolation: 26 October 2000 (Peter van Hoof)
	 *
	 * I computed the effective temperature for a random sample of interpolated
	 * atmospheres by integrating the flux as shown above and compared the results
	 * with the expected effective temperature using DELTA = (COMP-EXPEC)/EXPEC.
	 *
	 * I found that the average discrepancy was:
	 *
	 *     DELTA = -0.71% +/- 0.71% (sample size 5000)
	 *
	 * The most extreme discrepancies were
	 *     -4.37% <= DELTA <= 0.24%
	 *
	 * The most negative discrepancies were for Teff =  95 kK, log(g) = 5
	 * The most positive discrepancies were for Teff = 160 kK, log(g) = 8
	 *
	 * Since Cloudy checks the scaling elsewhere there is no need to re-scale 
	 * things here, but this inaccuracy should be kept in mind since it could
	 * indicate problems with the flux distribution */

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* WMBASICCompile rebin WMBASIC stellar models to match energy grid of code */
int WMBASICCompile(process_counter& pc)
{
	DEBUG_ENTRY( "WMBASICCompile()" );

	fprintf( ioQQQ, " WMBASICCompile on the job.\n" );

	/* define the major absorption edges that require special attention during rebinning
	 * see the routine CoStarCompile() for a more detailed discussion of this array */
	realnum Edges[3] = { 0.9994665f, 1.807140f, 3.999632f };

	access_scheme as = AS_LOCAL_ONLY_TRY;

	bool lgFail = false;
	if( lgFileReadable( "wmbasic.ascii", pc, as ) && !lgValidBinFile( "wmbasic.mod", pc, as ) )
		lgFail = lgCompileAtmosphere( "wmbasic.ascii", "wmbasic.mod", Edges, 3L, pc );
	return lgFail;
}

/* WMBASICInterpolate read in and interpolate on WMBASIC grid of hot star atmospheres */
long WMBASICInterpolate(double val[], /* val[nval] */
			long *nval,
			long *ndim,
			bool lgList,
			double *Tlow,
			double *Thigh)
{
	DEBUG_ENTRY( "WMBASICInterpolate()" );

	stellar_grid grid;
	grid.name = "wmbasic.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "     WMBASIC";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	return rfield.nflux_with_check;
}

/* CoStarInitialize create ascii file in Cloudy format */
STATIC bool CoStarInitialize(const char chFNameIn[],
			     const char chFNameOut[])
{
	DEBUG_ENTRY( "CoStarInitialize()" );

	/* read the original data file obtained off the web, 
	 * open as read only */
	FILE *ioIN;
	try
	{
		ioIN = open_data( chFNameIn, "r", AS_LOCAL_ONLY );
	}
	catch( cloudy_exit& )
	{
		return true;
	}

	/* get first line and see how many more to skip */
	long nskip;
	char chLine[INPUT_LINE_LENGTH];
	if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
	{
		fprintf( ioQQQ, " CoStarInitialize fails reading nskip.\n" );
		return true;
	}
	sscanf( chLine, "%li", &nskip );

	/* now skip the header information */
	for( long i=0; i < nskip; ++i )
	{
		if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
		{
			fprintf( ioQQQ, " CoStarInitialize fails skipping header.\n" );
			return true;
		}
	}

	/* now get number of models and number of wavelengths */
	long nModels, nWL;
	if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
	{
		fprintf( ioQQQ, " CoStarInitialize fails reading nModels, nWL.\n" );
		return true;
	}
	sscanf( chLine, "%li%li", &nModels, &nWL );

	if( nModels <= 0 || nWL <= 0 )
	{
		fprintf( ioQQQ, " CoStarInitialize scanned off impossible values for nModels=%li or nWL=%li\n",
			 nModels, nWL );
		return true;
	}

	/* this will hold all the model parameters */
	vector<mpp> telg(nModels);

	/* get all model parameters for the atmospheres */
	for( long i=0; i < nModels; ++i )
	{
		if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
		{
			fprintf( ioQQQ, " CoStarInitialize fails reading model parameters.\n" );
			return true;
		}
		/* first letter on line is indicator of grid */
		telg[i].chGrid = chLine[0];
		/* get the model id number */
		sscanf( chLine+1, "%i", &telg[i].modid );
		/* get the temperature */
		sscanf( chLine+23, "%lg", &telg[i].par[0] );
		/* get the surface gravity */
		sscanf( chLine+31, "%lg", &telg[i].par[1] );
		/* get the ZAMS mass */
		sscanf( chLine+7, "%lg", &telg[i].par[2] );
		/* get the model age */
		sscanf( chLine+15, "%lg", &telg[i].par[3] );

		/* the code in parse_table.cpp implicitly depends on this! */
		ASSERT( telg[i].par[2] > 10. );
		ASSERT( telg[i].par[3] > 10. );

		/* convert ZAMS masses to logarithms */
		telg[i].par[2] = log10(telg[i].par[2]);
	}

	/* this will be the file we create, that will be read to compute models, 
	 * open to write binary */
	FILE* ioOUT;
	try
	{
		ioOUT = open_data( chFNameOut, "w" );
	}
	catch( cloudy_exit& )
	{
		return true;
	}

	vector<string> names;
	names.push_back( "Teff" );
	names.push_back( "log(g)" );
	names.push_back( "log(M)" );
	names.push_back( "Age" );
	WriteASCIIHead(ioOUT, VERSION_COSTAR, 2, names, nModels, nWL-1, "lambda", 1., "F_nu", PI, telg, " %.6e", 1);

	/* get some workspace */
	vector<double> StarWavl(nWL);
	vector<double> StarFlux(nWL);

	/* get the star data */
	for( long i=0; i < nModels; ++i )
	{
		/* get number to skip */
		if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
		{
			fprintf( ioQQQ, " CoStarInitialize fails reading the skip to next spectrum.\n" );
			return true;
		}
		sscanf( chLine, "%li", &nskip );

		for( long j=0; j < nskip; ++j )
		{
			if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
			{
				fprintf( ioQQQ, " CoStarInitialize fails doing the skip.\n" );
				return true;
			}
		}

		/* now read in the wavelength and flux for this star */
		for( long j=0; j < nWL; ++j )
		{
			if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
			{
				fprintf( ioQQQ, " CoStarInitialize fails reading the spectral data.\n" );
				return true;
			}
			double help1, help2;
			sscanf( chLine, "%lg %lg", &help1, &help2 );

			/* lambda in Angstrom */
			if( i == 0 )
				StarWavl[j] = help1;
			else
				ASSERT( fp_equal(StarWavl[j], help1) );
			/* continuum flux is log of "astrophysical" flux in erg cm^-2 s^-1 Hz^-1 */
			StarFlux[j] = exp10(help2);

			/* sanity check */
			if( j > 0 )
				ASSERT( StarWavl[j] > StarWavl[j-1] );
		}

		// skip the last point as it appends a bogus RJ tail
		if( i == 0 )
			WriteASCIIData(ioOUT, StarWavl, nWL-1, "  %.6e", 4);
		WriteASCIIData(ioOUT, StarFlux, nWL-1, "  %.6e", 4);
	}

	fclose( ioIN );
	fclose( ioOUT );

	return false;
}

/* InterpolateGridCoStar read in and interpolate on costar grid of windy O atmospheres */
STATIC void InterpolateGridCoStar(stellar_grid *grid, /* struct with all the grid parameters */
				  const double val[], /* val[0]: Teff for imode = 1,2; M_ZAMS for imode = 3;
						       * age for imode = 4 */
				                      /* val[1]: nmodid for imode = 1; log(g) for imode = 2;
						       * age for imode = 3; M_ZAMS for imode = 4 */
				  double *val0_lo,
				  double *val0_hi)
{
	DEBUG_ENTRY( "InterpolateGridCoStar()" );

	long off;
	double lval[2];
	switch( grid->imode )
	{
	case IM_COSTAR_TEFF_MODID:
	case IM_COSTAR_TEFF_LOGG:
		lval[0] = val[0];
		lval[1] = val[1];
		off = 0;
		break;
	case IM_COSTAR_MZAMS_AGE:
		lval[0] = log10(val[0]); /* use log10(M_ZAMS) internally */
		lval[1] = val[1];
		off = 2;
		break;
	case IM_COSTAR_AGE_MZAMS:
		/* swap parameters, hence mimic IM_COSTAR_MZAMS_AGE */
		lval[0] = log10(val[1]); /* use log10(M_ZAMS) internally */
		lval[1] = val[0];
		off = 2;
		break;
	default:
		fprintf( ioQQQ, " InterpolateGridCoStar called with insane value for imode: %d.\n", grid->imode );
		cdEXIT(EXIT_FAILURE);
	}

	long nmodid = (long)(lval[1]+0.5);

	rfield.tNu[rfield.nShape].resize(rfield.nflux_with_check);
	rfield.tslop[rfield.nShape].resize(rfield.nflux_with_check);

	if( grid->lgASCII )
	{
		for( long i=0; i < rfield.nflux_with_check; ++i )
			rfield.tNu[rfield.nShape][i] = (realnum)rfield.anu(i);
	}
	else
	{
		/* read in the saved cloudy energy scale so we can confirm this is a good image */
		GetBins( grid, rfield.tNu[rfield.nShape] );

		/* check that the stored frequency mesh matches what is used in Cloudy */
		ValidateMesh( grid, rfield.tNu[rfield.nShape] );
	}

	if( DEBUGPRT )
	{
		/* check whether the models in the grid have the correct effective temperature */
		ValidateGrid( grid, 0.005 );
	}

	/* now allocate some temp workspace */
	vector<realnum> ValTr(grid->nTracks);
	vector<long> indloTr(grid->nTracks);
	vector<long> indhiTr(grid->nTracks);
	vector<long> index(2);

	/* first do horizontal search, i.e. search along individual tracks */
	for( long j=0; j < grid->nTracks; j++ )
	{
		if( grid->imode == IM_COSTAR_TEFF_MODID )
		{
			if( grid->trackLen[j] >= nmodid ) {
				index[0] = nmodid - 1;
				index[1] = j;
				long ptr = grid->jval[JIndex(grid,index)];
				indloTr[j] = ptr;
				indhiTr[j] = ptr;
				ValTr[j] = (realnum)grid->telg[ptr].par[off];
			}
			else
			{
				indloTr[j] = -2;
				indhiTr[j] = -2;
				ValTr[j] = -FLT_MAX;
			}
		}
		else
		{
			FindHCoStar( grid, j, lval[1], off, ValTr, indloTr, indhiTr );
		}
	}

	if( DEBUGPRT )
	{
		for( long j=0; j < grid->nTracks; j++ ) 
		{
			if( indloTr[j] >= 0 ) 
				printf( "track %c: models %c%d, %c%d, val %g\n",
					(char)('A'+j), grid->telg[indloTr[j]].chGrid, grid->telg[indloTr[j]].modid,
					grid->telg[indhiTr[j]].chGrid, grid->telg[indhiTr[j]].modid, ValTr[j]);
		}
	}

	long useTr[2], indlo[2], indhi[2];

	/* now do vertical search, i.e. interpolate between tracks */
	FindVCoStar( grid, lval[0], ValTr, useTr );

	/* This should only happen when InterpolateGridCoStar is called in non-optimizing mode,
	 * when optimizing InterpolateGridCoStar should report back to optimize_func()...
	 * The fact that FindVCoStar allows interpolation between non-adjoining tracks
	 * should guarantee that this will not happen. */
	if( useTr[0] < 0 )
	{
		fprintf( ioQQQ, " The parameters for the requested CoStar model are out of range.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	ASSERT( useTr[0] >= 0 && useTr[0] < grid->nTracks );
	ASSERT( useTr[1] >= 0 && useTr[1] < grid->nTracks );
	ASSERT( indloTr[useTr[0]] >= 0 && indloTr[useTr[0]] < (int)grid->nmods );
	ASSERT( indhiTr[useTr[0]] >= 0 && indhiTr[useTr[0]] < (int)grid->nmods );
	ASSERT( indloTr[useTr[1]] >= 0 && indloTr[useTr[1]] < (int)grid->nmods );
	ASSERT( indhiTr[useTr[1]] >= 0 && indhiTr[useTr[1]] < (int)grid->nmods );

	if( DEBUGPRT )
		printf( "interpolate between tracks %c and %c\n", (char)('A'+useTr[0]), (char)('A'+useTr[1]) );

	indlo[0] = indloTr[useTr[0]];
	indhi[0] = indhiTr[useTr[0]];
	indlo[1] = indloTr[useTr[1]];
	indhi[1] = indhiTr[useTr[1]];

	grid->index_list.push_back( indlo[0] );
	grid->index_list.push_back( indhi[0] );
	grid->index_list.push_back( indlo[1] );
	grid->index_list.push_back( indhi[1] );

	// read the necessary models
	SortUnique( grid->index_list, grid->index_list2 );
	grid->CloudyFlux.alloc(grid->index_list2.size(), rfield.nflux_with_check);
	if( grid->lgASCII )
	{
		if( lgReadAtmosphereTail(grid, Edges_CoStar, 3L, grid->index_list2) )
		{
			fprintf( ioQQQ, "Failed to read atmosphere models.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	for( size_t i=0; i < grid->index_list2.size(); ++i )
		GetModel( grid, grid->index_list2[i], &grid->CloudyFlux[i][0], lgVERBOSE, true );

	vector<double> aval(4);
	InterpolateModelCoStar( grid, lval, aval, indlo, indhi, index, 0, off, rfield.tslop[rfield.nShape] );

	for( long i=0; i < rfield.nflux_with_check; i++ )
	{
		rfield.tslop[rfield.nShape][i] = exp10(rfield.tslop[rfield.nShape][i]);
		if( rfield.tslop[rfield.nShape][i] < 1e-37 )
			rfield.tslop[rfield.nShape][i] = 0.;
	}

	if( false )
	{
		FILE *ioBUG = open_data( "interpolated.txt", "w" );
		for( long k=0; k < rfield.nflux_with_check; ++k )
			fprintf( ioBUG, "%e %e\n", rfield.tNu[rfield.nShape][k].Ryd(), rfield.tslop[rfield.nShape][k] );
		fclose( ioBUG );
	}

	/* sanity check: see whether this model has the correct effective temperature */
	if( ! lgValidModel( rfield.tNu[rfield.nShape], rfield.tslop[rfield.nShape], aval[0], 0.05 ) )
		TotalInsanity();

	/* set limits for optimizer */
	vector<long> dum;
	SetLimits( grid, lval[0], dum, dum, useTr, ValTr, val0_lo, val0_hi );

	/* now write some final info */
	if( called.lgTalk )
	{
		fprintf( ioQQQ, "                       * #<< FINAL: T_eff = %7.1f, ", aval[0] );
		fprintf( ioQQQ, "log(g) = %4.2f, M(ZAMS) = %5.1f, age = ", aval[1], exp10(aval[2]) );
		fprintf( ioQQQ, PrintEfmt("%8.2e",aval[3]) );
		fprintf( ioQQQ, "  >>> *\n" );
	}
}

/* find which models to use for interpolation along a given evolutionary track */
STATIC void FindHCoStar(const stellar_grid *grid,
			long track,
			double par2,           /* requested log(g) or age */
			long off,              /* determines which parameter to match 0 -> log(g), 2 -> age */
			vector<realnum>& ValTr,/* ValTr[track]: Teff/log(M) value for interpolated model along track */
			vector<long>& indloTr, /* indloTr[track]: model number for first model used in interpolation */
			vector<long>& indhiTr) /* indhiTr[track]: model number for second model used in interpolation */
{
	DEBUG_ENTRY( "FindHCoStar()" );

	indloTr[track] = -2;
	indhiTr[track] = -2;
	ValTr[track] = -FLT_MAX;

	long mod1, mod2;
	vector<long> index(2); 
	index[1] = track;

	for( long j=0; j < grid->trackLen[track]; j++ )
	{
		index[0] = j;
		mod1 = grid->jval[JIndex(grid,index)];

		/* do we have an exact match ? */
		if( fabs(par2-grid->telg[mod1].par[off+1]) <= 10.*FLT_EPSILON*fabs(grid->telg[mod1].par[off+1]) )
		{
			indloTr[track] = mod1;
			indhiTr[track] = mod1;
			ValTr[track] = (realnum)grid->telg[mod1].par[off];
			return;
		}
	}

	for( long j=0; j < grid->trackLen[track]-1; j++ )
	{
		index[0] = j;
		mod1 = grid->jval[JIndex(grid,index)];
		index[0] = j+1;
		mod2 = grid->jval[JIndex(grid,index)];

		/* do we interpolate ? */
		if( (par2 - grid->telg[mod1].par[off+1])*(par2 - grid->telg[mod2].par[off+1]) < 0. )
		{
			double frac;

			indloTr[track] = mod1;
			indhiTr[track] = mod2;
			frac = (par2 - grid->telg[mod2].par[off+1])/
				(grid->telg[mod1].par[off+1] - grid->telg[mod2].par[off+1]);
			ValTr[track] = (realnum)(frac*grid->telg[mod1].par[off] + 
				(1.-frac)*grid->telg[mod2].par[off] );
			break;
		}
	}
}

/* find which tracks to use for interpolation in between tracks */
STATIC void FindVCoStar(const stellar_grid *grid,
			double par1,  /* requested Teff or ZAMS mass */
			vector<realnum>& ValTr, /* internal workspace */
			long useTr[]) /* useTr[0]: track number for first track to be used in interpolation
				       *            (i.e., 0 means 'A', etc.)
				       * useTr[1]: track number for second track to be used in interpolation
				       * NOTE: FindVCoStar raises a flag when interpolating between non-adjoining
				       *       tracks, i.e. when (useTr[1]-useTr[0]) > 1 */
{
	DEBUG_ENTRY( "FindVCoStar()" );

	useTr[0] = -1;
	useTr[1] = -1;

	for( long j=0; j < grid->nTracks; j++ )
	{
		/* do we have an exact match ? */
		if( ValTr[j] != -FLT_MAX && fabs(par1-(double)ValTr[j]) <= 10.*FLT_EPSILON*fabs(ValTr[j]) )
		{
			useTr[0] = j;
			useTr[1] = j;
			break;
		}
	}

	if( useTr[0] >= 0 )
		return;

	for( long j=0; j < grid->nTracks-1; j++ )
	{
		if( ValTr[j] != -FLT_MAX )
		{
			/* find next valid track */
			long j2 = 0;
			for( long i = j+1; i < grid->nTracks; i++ )
			{
				if( ValTr[i] != -FLT_MAX )
				{
					j2 = i;
					break;
				}
			}

			/* do we interpolate ? */
			if( j2 > 0 && ((realnum)par1-ValTr[j])*((realnum)par1-ValTr[j2]) < 0.f )
			{
				useTr[0] = j;
				useTr[1] = j2;
				break;
			}
		}
	}

	/* raise caution when we interpolate between non-adjoining tracks */
	continuum.lgCoStarInterpolationCaution = ( useTr[1]-useTr[0] > 1 );
}

/* Make a listing of all the models in the CoStar grid */
STATIC void CoStarListModels(const stellar_grid *grid)
{
	DEBUG_ENTRY( "CoStarListModels()" );

	long maxlen = 0;
	for( long n=0; n < grid->nTracks; n++ )
		maxlen = MAX2( maxlen, grid->trackLen[n] );

	fprintf( ioQQQ, "\n" );
	fprintf( ioQQQ, "  Track\\Index |" );
	for( long n = 0; n < maxlen; n++ )
		fprintf( ioQQQ, "     %5ld      ", n+1 );
	fprintf( ioQQQ, "\n" );
	fprintf( ioQQQ, "--------------|" );
	for( long n = 0; n < maxlen; n++ )
		fprintf( ioQQQ, "----------------" );
	fprintf( ioQQQ, "\n" );

	vector<long> index(2);
	for( index[1]=0; index[1] < grid->nTracks; ++index[1] )
	{
		long ptr;
		double Teff, alogg, Mass;

		fprintf( ioQQQ, " %c", (char)('A'+index[1]) );
		index[0] = 0;
		ptr = grid->jval[JIndex(grid,index)];
		Mass = exp10(grid->telg[ptr].par[2]);
		fprintf( ioQQQ, " (%3.0f Msol) |", Mass );

		for( index[0]=0; index[0] < grid->trackLen[index[1]]; ++index[0] )
		{
			ptr = grid->jval[JIndex(grid,index)];
			Teff = grid->telg[ptr].par[0];
			alogg = grid->telg[ptr].par[1];
			fprintf( ioQQQ, "  (%6.1f,%4.2f)", Teff, alogg );
		}
		fprintf( ioQQQ, "\n" );
	}
}

/*  RauchInitialize does the actual work of preparing the ascii file */
STATIC bool RauchInitialize(const char chFName[],
			    const char chSuff[],
			    const vector<mpp>& telg,
			    long nmods,
			    long n,
			    long ngrids,
			    const double par2[], /* par2[ngrids] */
			    int format)
{
	DEBUG_ENTRY( "RauchInitialize()" );

	/* grab some space for the wavelengths and fluxes */
	vector<double> wavl(NRAUCH);
	vector<double> fluxes(NRAUCH);

	FILE *ioOut;
	try
	{
		if( n == 1 )
			ioOut = open_data( chFName, "w" );
		else
			ioOut = open_data( chFName, "a" );
	}
	catch( cloudy_exit& )
	{
		return true;
	}

	if( n == 1 )
	{
		long ndim = ( ngrids == 1 ) ? 2 : 3;
		vector<string> names;
		names.push_back( "Teff" );
		names.push_back( "log(g)" );
		if( ngrids == 2 )
			names.push_back( "log(Z)" );
		else if( ngrids == 11 )
			names.push_back( "f(He)" );
		else if( ngrids != 1 )
			TotalInsanity();
		// make local copy so that we can fill in par2[]...
		vector<mpp> mytelg(ngrids*telg.size());
		size_t k = 0;
		/* NB - this is based on the assumption that each of the planes in the cubic grid is the same */
		for( long j=0; j < ngrids; j++ )
			for( size_t i=0; i < telg.size(); i++ )
			{
				mytelg[k] = telg[i];
				if( ngrids > 1 )
					mytelg[k].par[2] = par2[j];
				++k;
			}
		/* Rauch models give the "Astrophysical" flux F_lambda in erg/cm^2/s/cm */
		/* the factor PI*1e-8 is needed to convert to "regular" flux in erg/cm^2/s/Angstrom */
		WriteASCIIHead(ioOut, VERSION_ASCII, ndim, names, nmods*ngrids, NRAUCH, "lambda", 1.,
			       "F_lambda", PI*1.e-8, mytelg, " %.1f", 4);
	}

	for( long i=0; i < nmods; i++ )
	{
		char chLine[INPUT_LINE_LENGTH];
		/* must create name of next stellar atmosphere */
		if( format == 1 )
			sprintf( chLine, "%7.7ld_%2ld", (long)(telg[i].par[0]+0.5), (long)(10.*telg[i].par[1]+0.5) );
		else if( format == 2 )
			sprintf( chLine, "%7.7ld_%.2f", (long)(telg[i].par[0]+0.5), telg[i].par[1] );
		else
			TotalInsanity();
		string chFileName( chLine );
		chFileName += chSuff;
		/* now open next stellar atmosphere for reading*/
		FILE *ioIn;
		try
		{
			ioIn = open_data( chFileName.c_str(), "r", AS_LOCAL_ONLY );
		}
		catch( cloudy_exit& )
		{
			return true;
		}

		/* get first line */
		long j = 0;
		if( read_whole_line( chLine, (int)sizeof(chLine), ioIn ) == NULL )
		{
			fprintf( ioQQQ, " RauchInitialize error in atmosphere file %ld %ld\n", i, j );
			return true;
		}
		/* >>chng 02 nov 20, now keep reading them until don't hit the *
		 * since number of comments may change */
		while( chLine[0] == '*' )
		{
			if( read_whole_line( chLine, (int)sizeof(chLine), ioIn ) == NULL )
			{
				fprintf( ioQQQ, " RauchInitialize error in atmosphere file %ld %ld\n", i, j );
				return true;
			}
			++j;
		}

		for( j=0; j < NRAUCH; j++ )
		{
			double ttemp, wl;
			/* get the input line */
			/* >>chng 02 nov 20, don't reread very first line image since we got it above */
			if( j > 0 )
			{
				if(read_whole_line( chLine, (int)sizeof(chLine), ioIn )==NULL )
				{
					fprintf( ioQQQ, " RauchInitialize error in atmosphere file %ld %ld\n", i, j );
					return true;
				}
			}

			/* scan off wavelength and flux)*/
			if( sscanf( chLine, "%lf %le", &wl, &ttemp ) != 2 )
			{
				fprintf( ioQQQ, " RauchInitialize error in atmosphere file %ld %ld\n", i, j );
				return true;
			}

			if( i == 0 )
				wavl[j] = wl;
			else
			{
				/* check if this model is on the same wavelength grid as the first */
				if( !fp_equal(wavl[j],wl,10) )
				{
					fprintf( ioQQQ, " RauchInitialize error in atmosphere file %ld %ld\n", i, j );
					return true;
				}
			}
			fluxes[j] = ttemp; 
		}

		/* finished - close the unit */
		fclose(ioIn);

		/* now write to output file */
		if( i == 0 && n == 1 )
			WriteASCIIData(ioOut, wavl, NRAUCH, "  %.4e", 5);
		WriteASCIIData(ioOut, fluxes, NRAUCH, "  %.4e", 5);
	}

	fclose(ioOut);

	return false;
}

STATIC void RauchReadMPP(vector<mpp>& telg1,
			 vector<mpp>& telg2,
			 vector<mpp>& telg3,
			 vector<mpp>& telg4,
			 vector<mpp>& telg5,
			 vector<mpp>& telg6)
{
	DEBUG_ENTRY( "RauchReadMPP()" );

	const char fnam[] = "rauch_models.dat";
	fstream ioDATA;
	open_data( ioDATA, fnam, mode_r );

	string line;
	getdataline( ioDATA, line );
	long version;
	istringstream iss( line );
	iss >> version;
	if( version != VERSION_RAUCH_MPP )
	{
		fprintf( ioQQQ, " RauchReadMPP: the version of %s is not the current version.\n", fnam );
		fprintf( ioQQQ, " Please obtain the current version from the Cloudy web site.\n" );
		fprintf( ioQQQ, " I expected to find version %ld and got %ld instead.\n",
			 VERSION_RAUCH_MPP, version );
		cdEXIT(EXIT_FAILURE);
	}

	getdataline( ioDATA, line );
	unsigned long ndata;
	istringstream iss2( line );
	iss2 >> ndata;
	ASSERT( ndata == telg1.size() );
	// this implicitly assumes there is exactly one comment line between
	// the number of data points and the start of the data
	getline( ioDATA, line );
	// read data for H-Ca grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg1[i].par[0] >> telg1[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss3( line );
	iss3 >> ndata;
	ASSERT( ndata == telg2.size() );
	getline( ioDATA, line );
	// read data for H-Ni grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg2[i].par[0] >> telg2[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss4( line );
	iss4 >> ndata;
	ASSERT( ndata == telg3.size() );
	getline( ioDATA, line );
	// read data for PG1159 grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg3[i].par[0] >> telg3[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss5( line );
	iss5 >> ndata;
	ASSERT( ndata == telg4.size() );
	getline( ioDATA, line );
	// read data for pure H grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg4[i].par[0] >> telg4[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss6( line );
	iss6 >> ndata;
	ASSERT( ndata == telg5.size() );
	getline( ioDATA, line );
	// read data for pure He grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg5[i].par[0] >> telg5[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss7( line );
	iss7 >> ndata;
	ASSERT( ndata == telg6.size() );
	getline( ioDATA, line );
	// read data for pure H+He grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg6[i].par[0] >> telg6[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss8( line );
	iss8 >> version;
	ASSERT( version == VERSION_RAUCH_MPP );
}

inline void getdataline(fstream& ioDATA,
			string& line)
{
	do
	{
		getline( ioDATA, line );
	}
	while( line[0] == '#' );
}

STATIC void WriteASCIIHead(FILE* ioOut,
			   long version,
			   long ndim,
			   const vector<string>& names,
			   long nmods,
			   long ngrid,
			   const string& wtype,
			   double wfac,
			   const string& ftype,
			   double ffac,
			   const vector<mpp>& telg,
			   const char* format,
			   int nmult)
{
	DEBUG_ENTRY( "WriteASCIIHead()" );

	long npar = (long)names.size();

	ASSERT( ndim <= npar && npar <= MDIM );
	ASSERT( wtype == "lambda" || wtype == "nu" );
	ASSERT( ftype == "F_lambda" || ftype == "F_nu" ||
		ftype == "H_lambda" || ftype == "H_nu" );
	ASSERT( version == VERSION_ASCII || version == VERSION_COSTAR );

	fprintf( ioOut, "  %ld\n", version );
	fprintf( ioOut, "  %ld\n", ndim );
	fprintf( ioOut, "  %ld\n", npar );
	for( long i=0; i < npar; ++i )
		fprintf( ioOut, "  %s\n", names[i].c_str() );
	fprintf( ioOut, "  %ld\n", nmods );
	fprintf( ioOut, "  %ld\n", ngrid );
	fprintf( ioOut, "  %s\n", wtype.c_str() );
	fprintf( ioOut, "  %.8e\n", wfac );
	fprintf( ioOut, "  %s\n", ftype.c_str() );
	fprintf( ioOut, "  %.8e\n", ffac );
	/* write out the parameter grid */
	long i;
	for( i=0; i < nmods; i++ )
	{
		fprintf( ioOut, " " );
		for( long j=0; j < npar; ++j )
			fprintf( ioOut, format, telg[i].par[j] );
		if( version == VERSION_COSTAR )
			fprintf( ioOut, " %c%d", telg[i].chGrid, telg[i].modid );
		if( ((i+1)%nmult) == 0 )
			fprintf( ioOut, "\n" );
	}
	if( (i%nmult) != 0 )
		fprintf( ioOut, "\n" );
}

STATIC void WriteASCIIData(FILE* ioOut,
			   const vector<double>& data, // data[ngrid]
			   long ngrid,
			   const char* format,
			   int nmult)
{
	DEBUG_ENTRY( "WriteASCIIData()" );

	long i;
	for( i=0; i < ngrid; i++ )
	{
		fprintf( ioOut, format, data[i] );
		if( ((i+1)%nmult) == 0 )
			fprintf( ioOut, "\n" );
	}
	if( (i%nmult) != 0 )
		fprintf( ioOut, "\n" );
}

STATIC bool lgReadAtmosphereHead(stellar_grid* grid)
{
	DEBUG_ENTRY( "lgReadAtmosphereHead()" );

	// skip leading comments
	char chLine[INPUT_LINE_LENGTH];
	do
	{
		if( read_whole_line( chLine, (int)sizeof(chLine), grid->ioIN ) == NULL )
		{
			fprintf( ioQQQ, " ReadAtmosphere fails reading line.\n" );
			return true;
		}
	}
	while( chLine[0] == '#' );

	/* read version number */
	long version;
	if( sscanf( chLine, "%ld", &version ) != 1 )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading VERSION.\n" );
		return true;
	}

	if( version != VERSION_ASCII && version != VERSION_COSTAR )
	{
		fprintf( ioQQQ, " ReadAtmosphere: there is a version number mismatch in"
			 " the ascii atmosphere file: %s.\n", grid->name.c_str() );
		fprintf( ioQQQ, " ReadAtmosphere: Please recreate this file or download the"
			 " latest version following the instructions on the Cloudy website.\n" );
		return true;
	}

	/* >>chng 06 jun 10, read the dimension of the grid, PvH */
	if( fscanf( grid->ioIN, "%d", &grid->ndim ) != 1 || grid->ndim <= 0 || grid->ndim > MDIM )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading valid dimension of grid.\n" );
		return true;
	}

	/* >>chng 06 jun 12, read the number of model parameters, PvH */
	if( fscanf( grid->ioIN, "%d", &grid->npar ) != 1 || grid->npar <= 0 ||
	    grid->npar < grid->ndim || grid->npar > MDIM )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading valid no. of model parameters.\n" );
		return true;
	}

	for( long nd=0; nd < grid->npar; nd++ )
	{
		if( fscanf( grid->ioIN, "%6s", grid->names[nd] ) != 1 )
		{
			fprintf( ioQQQ, " ReadAtmosphere failed reading parameter label.\n" );
			return true;
		}
	}

	/* >>chng 05 nov 18, read the following extra parameters from the ascii file, PvH */
	if( fscanf( grid->ioIN, "%d", &grid->nmods ) != 1 || grid->nmods <= 0 )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading valid number of models.\n" );
		return true;
	}

	if( fscanf( grid->ioIN, "%d", &grid->ngrid ) != 1 || grid->ngrid <= 1 )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading valid number of grid points.\n" );
		return true;
	}

	char chDataType[11];
	/* read data type for wavelengths, allowed values are lambda, nu */
	if( fscanf( grid->ioIN, "%10s", chDataType ) != 1 )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading wavl DataType string.\n" );
		return true;
	}

	if( strcmp( chDataType, "lambda" ) == 0 )
		grid->lgFreqX = false;
	else if( strcmp( chDataType, "nu" ) == 0 )
		grid->lgFreqX = true;
	else {
		fprintf( ioQQQ, " ReadAtmosphere found illegal wavl DataType: %s.\n", chDataType );
		return true;
	}

	if( fscanf( grid->ioIN, "%le", &grid->convert_wavl ) != 1 || grid->convert_wavl <= 0. )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading valid wavl conversion factor.\n" );
		return true;
	}

	/* read data type for flux, allowed values F_lambda, H_lambda, F_nu, H_nu */
	if( fscanf( grid->ioIN, "%10s", chDataType ) != 1 )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading flux DataType string.\n" );
		return true;
	}

	if( strcmp( chDataType, "F_lambda" ) == 0 || strcmp( chDataType, "H_lambda" ) == 0 )
		grid->lgFreqY = false;
	else if( strcmp( chDataType, "F_nu" ) == 0 || strcmp( chDataType, "H_nu" ) == 0 )
		grid->lgFreqY = true;
	else {
		fprintf( ioQQQ, " ReadAtmosphere found illegal flux DataType: %s.\n", chDataType );
		return true;
	}

	if( fscanf( grid->ioIN, "%le", &grid->convert_flux ) != 1 || grid->convert_flux <= 0. )
	{
		fprintf( ioQQQ, " ReadAtmosphere failed reading valid flux conversion factor.\n" );
		return true;
	}

	grid->telg.resize(grid->nmods);

	for( long i=0; i < grid->nmods; i++ )
	{
		for( long nd=0; nd < grid->npar; nd++ )
		{
			if( fscanf( grid->ioIN, "%le", &grid->telg[i].par[nd] ) != 1 )
			{
				fprintf( ioQQQ, " ReadAtmosphere failed reading valid model parameter.\n" );
				return true;
			}
		}
		if( version == VERSION_COSTAR )
		{
			if( fscanf( grid->ioIN, " %c%d", &grid->telg[i].chGrid, &grid->telg[i].modid ) != 2 )
			{
				fprintf( ioQQQ, " ReadAtmosphere failed reading valid track ID.\n" );
				return true;
			}
		}
	}
	return false;
}

STATIC bool lgReadAtmosphereTail(stellar_grid* grid,
				 const realnum Edges[], /* Edges[nEdges] */
				 long nEdges,
				 const vector<long>& index_list)
{
	DEBUG_ENTRY( "lgReadAtmosphereTail()" );

	/* get some workspace */
	vector<realnum> StarEner(grid->ngrid);
	vector<realnum> scratch(grid->ngrid);
	vector<realnum> StarFlux(grid->ngrid);

	/* read wavelength grid */
	for( long i=0; i < grid->ngrid; i++ )
	{
		double help;
		if( fscanf( grid->ioIN, "%lg", &help ) != 1 )
		{
			fprintf( ioQQQ, " ReadAtmosphere failed reading wavelength.\n" );
			return true;
		}
		/* this conversion makes sure that scratch[i] is
		 * either wavelength in Angstrom or frequency in Hz */
		scratch[i] = (realnum)(help*grid->convert_wavl);

		if( scratch[i] <= 0.f )
		{
			fprintf( ioQQQ, " PROBLEM: a non-positive %s was found, value: %e\n",
				 grid->lgFreqX ? "frequency" : "wavelength", scratch[i] );
			cdEXIT(EXIT_FAILURE);
		}
	}

	bool lgFlip = ( !grid->lgFreqX && scratch[0] < scratch[1] ) || ( grid->lgFreqX && scratch[0] > scratch[1] );

	/* convert continuum over to increasing frequency in Ryd */
	for( long i=0; i < grid->ngrid; i++ )
	{
		/* convert scratch[i] to frequency in Ryd */
		if( grid->lgFreqX )
			scratch[i] /= (realnum)FR1RYD;
		else
			scratch[i] = (realnum)(RYDLAM/scratch[i]);

		if( lgFlip )
			StarEner[grid->ngrid-i-1] = scratch[i];
		else
			StarEner[i] = scratch[i];
	}

	ASSERT( StarEner[0] > 0.f );
	/* make sure the array is in ascending order */
	for( long i=1; i < grid->ngrid; i++ )
	{
		if( StarEner[i] <= StarEner[i-1] )
		{
			fprintf( ioQQQ, " PROBLEM: the %s grid is not strictly monotonically increasing/decreasing\n",
				 grid->lgFreqX ? "frequency" : "wavelength" );
			cdEXIT(EXIT_FAILURE);	
		}
	}

	size_t ni = 0;
	size_t ni_end = ( index_list.size() == 0 ) ? grid->nmods : index_list.size();

	for( long imod=0; imod < grid->nmods; imod++ )
	{
		const realnum CONVERT_FNU = (realnum)(1.e8*SPEEDLIGHT/POW2(FR1RYD));

		/* now read the stellar fluxes */
		for( long i=0; i < grid->ngrid; i++ )
		{
			double help;
			if( fscanf( grid->ioIN, "%lg", &help ) != 1 )
			{
				fprintf( ioQQQ, " ReadAtmosphere failed reading star flux.\n" );
				return true;
			}
			/* this conversion makes sure that scratch[i] is either
			 * F_nu in erg/cm^2/s/Hz or F_lambda in erg/cm^2/s/A */
			scratch[i] = (realnum)(help*grid->convert_flux);

			/* this can underflow on the Wien tail */
			if( scratch[i] < 0.f )
			{
				fprintf( ioQQQ, "\n PROBLEM: a negative flux was found, model number %ld, value: %e\n",
					 imod+1, help );
				cdEXIT(EXIT_FAILURE);
			}
		}

		for( long i=0; i < grid->ngrid; i++ )
		{
			if( lgFlip )
				StarFlux[grid->ngrid-i-1] = scratch[i];
			else
				StarFlux[i] = scratch[i];
		}

		for( long i=0; i < grid->ngrid; i++ )
		{
			/* this converts to F_nu in erg/cm^2/s/Hz */
			if( !grid->lgFreqY )
				StarFlux[i] *= CONVERT_FNU/POW2(StarEner[i]);
			ASSERT( StarFlux[i] >= 0.f );
		}

		if( false )
		{
			DumpAtmosphere( "atmosphere_input_dump.txt", imod, grid->npar, grid->names, grid->telg,
					grid->ngrid, get_ptr(StarEner), get_ptr(StarFlux) );
		}

		if( index_list.size() == 0 || imod == index_list[ni] )
		{
			/* the re-binned values are returned in the "CloudyFlux" array */
			RebinAtmosphere(StarEner, StarFlux, grid->ngrid, Edges, nEdges, &grid->CloudyFlux[ni][0]);

			if( false )
			{
				DumpAtmosphere( "atmosphere_output_dump.txt", imod, grid->npar, grid->names, grid->telg,
						rfield.nflux_with_check, rfield.anuptr(), &grid->CloudyFlux[ni][0] );
			}

			++ni;
		}

		if( ni == ni_end )
			break;
	}
	return false;
}
			   
/* lgCompileAtmosphere does the actual rebinning onto the Cloudy grid and writes the binary file */
/* >>chng 01 feb 12, added return value to indicate success (0) or failure (1) */
STATIC bool lgCompileAtmosphere(const char chFNameIn[],
				const char chFNameOut[],
				const realnum Edges[], /* Edges[nEdges] */
				long nEdges,
				process_counter& pc)
{
	DEBUG_ENTRY( "lgCompileAtmosphere()" );

	++pc.nFail; // claim failure here in case we exit early
	stellar_grid grid;
	grid.name = chFNameIn;
	try
	{
		grid.ioIN = open_data( chFNameIn, "r", AS_LOCAL_ONLY );
	}
	catch( cloudy_exit& )
	{
		return true;
	}
	fprintf( ioQQQ, " lgCompileAtmosphere got %s.\n", chFNameIn );

	if( lgReadAtmosphereHead(&grid) )
		return true;

	FILE *ioOUT;
	try
	{
		ioOUT = open_data( chFNameOut, "wb" );
	}
	catch( cloudy_exit& )
	{
		return true;
	}

	int32 val[7];
	uint32 uval[2];
	double dval[3];
	char md5sum[NMD5];

	val[0] = (int32)VERSION_BIN;
	val[1] = (int32)MDIM;
	val[2] = (int32)MNAM;
	val[3] = (int32)grid.ndim;
	val[4] = (int32)grid.npar;
	val[5] = (int32)grid.nmods;
	val[6] = (int32)rfield.nflux_with_check;
	uval[0] = sizeof(val) + sizeof(uval) + sizeof(dval) + sizeof(md5sum) +
		sizeof(grid.names) + grid.nmods*sizeof(mpp); /* nOffset */
	uval[1] = rfield.nflux_with_check*sizeof(realnum); /* nBlocksize */
	dval[0] = rfield.emm();
	dval[1] = rfield.egamry();
	dval[2] = rfield.getResolutionScaleFactor();

	for( unsigned int i=0; i < NMD5; i++ )
		md5sum[i] = rfield.mesh_md5sum()[i];

	vector<realnum> SaveAnu(rfield.nflux_with_check);
	for( long i=0; i < rfield.nflux_with_check; ++i )
		SaveAnu[i] = (realnum)rfield.anu(i);

	if( fwrite( val, sizeof(val), 1, ioOUT ) != 1 ||
	    fwrite( uval, sizeof(uval), 1, ioOUT ) != 1 ||
	    /* write out the lower, upper bound of the energy mesh, and the res scale factor */
	    fwrite( dval, sizeof(dval), 1, ioOUT ) != 1 ||
	    /* write out the (modified) md5 checksum of continuum_mesh.ini */
	    fwrite( md5sum, sizeof(md5sum), 1, ioOUT ) != 1 ||
	    fwrite( grid.names, sizeof(grid.names), 1, ioOUT ) != 1 ||
	    /* write out the array of {Teff,log(g)} pairs */
	    fwrite( get_ptr(grid.telg), sizeof(mpp), (size_t)grid.nmods, ioOUT ) != (size_t)grid.nmods ||
	    /* write out the cloudy energy grid for later sanity checks */
	    fwrite( get_ptr(SaveAnu), (size_t)uval[1], 1, ioOUT ) != 1 )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed writing header of output file.\n" );
		return true;
	}

	vector<long> empty;
	grid.CloudyFlux.alloc(grid.nmods, rfield.nflux_with_check);
	if( lgReadAtmosphereTail(&grid, Edges, nEdges, empty) )
		return true;

	for( long imod=0; imod < grid.nmods; imod++ )
	{
		/* write the continuum out as a binary file */
		if( fwrite( &grid.CloudyFlux[imod][0], (size_t)uval[1], 1, ioOUT ) != 1 )
		{
			fprintf( ioQQQ, " lgCompileAtmosphere failed writing star flux.\n" );
			return true;
		}
	}

	fclose(ioOUT);

	fprintf( ioQQQ, " lgCompileAtmosphere completed ok.\n\n" );

	--pc.nFail; // already registered failure at the start -> revert this
	++pc.nOK;
	return false;
}

STATIC void InitGrid(stellar_grid *grid,
		     bool lgList,
		     bool lgASCII)
{
	DEBUG_ENTRY( "InitGrid()" );

	try
	{
		const char* mode = lgASCII ? "r" : "rb";
		grid->ioIN = open_data( grid->name.c_str(), mode, grid->scheme );
	}
	catch( cloudy_exit& )
	{
		/* something went wrong */
		/* NB NB - DO NOT CHANGE THE FOLLOWING ERROR MESSAGE! checkall.pl picks it up */
		const char* compilemsg = lgASCII ? "" : " and compiled with the COMPILE STARS command";
		fprintf( ioQQQ, " Error: stellar atmosphere file not found.\n" );
		fprintf( ioQQQ, "\n\n If the path is set then it is possible that the stellar"
			 " atmosphere data files do not exist.\n");
		fprintf( ioQQQ, " Have the stellar data files been downloaded%s?\n", compilemsg );
		fprintf( ioQQQ, " If you are simply running the test suite and do not need the"
			 " stellar continua then you should simply ignore this failure\n");
		cdEXIT(EXIT_FAILURE);
	}

	if( lgASCII )
		InitGridASCII(grid);
	else
		InitGridBin(grid);

	grid->val.alloc(grid->ndim,grid->nmods);
	grid->nval.resize(grid->ndim);

	grid->lgIsTeffLoggGrid = ( grid->ndim >= 2 &&
				   strcmp( grid->names[0], "Teff" ) == 0 &&
				   strcmp( grid->names[1], "log(g)" ) == 0 );

	InitIndexArrays( grid, lgList );

	grid->lgASCII = lgASCII;
	/* set default interpolation mode */
	grid->imode = IM_RECT_GRID;
	/* these are only used by CoStar grids */
	grid->nTracks = 0;
}

inline void InitGridASCII(stellar_grid *grid)
{
	if( lgReadAtmosphereHead(grid) )
	{
		fprintf( ioQQQ, " InitGrid failed reading header.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	rfield.RSFCheck[rfield.nShape] = rfield.getResolutionScaleFactor();
}

STATIC void InitGridBin(stellar_grid *grid)
{
	DEBUG_ENTRY( "InitGridBin()" );

	int32 version, mdim, mnam;
	double mesh_elo, mesh_ehi;
	char md5sum[NMD5];
	/* >>chng 01 oct 17, add version and size to this array */
	if( fread( &version, sizeof(version), 1, grid->ioIN ) != 1 ||
	    fread( &mdim, sizeof(mdim), 1, grid->ioIN ) != 1 ||
	    fread( &mnam, sizeof(mnam), 1, grid->ioIN ) != 1 ||
	    fread( &grid->ndim, sizeof(grid->ndim), 1, grid->ioIN ) != 1 ||
	    fread( &grid->npar, sizeof(grid->npar), 1, grid->ioIN ) != 1 ||
	    fread( &grid->nmods, sizeof(grid->nmods), 1, grid->ioIN ) != 1 ||
	    fread( &grid->ngrid, sizeof(grid->ngrid), 1, grid->ioIN ) != 1 ||
	    fread( &grid->nOffset, sizeof(grid->nOffset), 1, grid->ioIN ) != 1 ||
	    fread( &grid->nBlocksize, sizeof(grid->nBlocksize), 1, grid->ioIN ) != 1 ||
	    fread( &mesh_elo, sizeof(mesh_elo), 1, grid->ioIN ) != 1 ||
	    fread( &mesh_ehi, sizeof(mesh_ehi), 1, grid->ioIN ) != 1 ||
	    fread( &rfield.RSFCheck[rfield.nShape], sizeof(rfield.RSFCheck[rfield.nShape]), 1, grid->ioIN ) != 1 ||
	    fread( md5sum, sizeof(md5sum), 1, grid->ioIN ) != 1 )
	{
		fprintf( ioQQQ, " InitGrid failed reading header.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* do some sanity checks */
	if( version != VERSION_BIN )
	{
		fprintf( ioQQQ, " InitGrid: there is a version mismatch between"
			 " the compiled atmospheres file I expected and the one I found.\n" );
		fprintf( ioQQQ, " InitGrid: Please recompile the stellar"
			 " atmospheres file with the command: %s.\n", grid->command.c_str() );
		cdEXIT(EXIT_FAILURE);
	}

	if( mdim != MDIM || mnam != MNAM )
	{
		fprintf( ioQQQ, " InitGrid: the compiled atmospheres file is produced"
			 " with an incompatible version of Cloudy.\n" );
		fprintf( ioQQQ, " InitGrid: Please recompile the stellar"
			 " atmospheres file with the command: %s.\n", grid->command.c_str() );
		cdEXIT(EXIT_FAILURE);
	}

	if( !fp_equal_tol( rfield.emm(), mesh_elo, 1.e-11*rfield.emm() ) ||
	    !fp_equal_tol( rfield.egamry(), mesh_ehi, 1.e-7*rfield.egamry() ) ||
	    strncmp( rfield.mesh_md5sum().c_str(), md5sum, NMD5 ) != 0 ||
	    rfield.nflux_with_check != grid->ngrid )
	{
		fprintf( ioQQQ, " InitGrid: the compiled atmospheres file is produced"
			 " with an incompatible frequency grid.\n" );
		fprintf( ioQQQ, " InitGrid: Please recompile the stellar"
			 " atmospheres file with the command: %s.\n", grid->command.c_str() );
		cdEXIT(EXIT_FAILURE);
	}

	ASSERT( grid->ndim > 0 && grid->ndim <= MDIM );
	ASSERT( grid->npar >= grid->ndim && grid->npar <= MDIM );
	ASSERT( grid->nmods > 0 );
	ASSERT( grid->ngrid > 0 );
	ASSERT( grid->nOffset > 0 );
	ASSERT( grid->nBlocksize > 0 );

	if( fread( &grid->names, sizeof(grid->names), 1, grid->ioIN ) != 1 )
	{
		fprintf( ioQQQ, " InitGrid failed reading names array.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	grid->telg.resize(grid->nmods);

	if( fread( get_ptr(grid->telg), sizeof(mpp), grid->nmods, grid->ioIN ) != (size_t)grid->nmods )
	{
		fprintf( ioQQQ, " InitGrid failed reading model parameter block.\n" );
		cdEXIT(EXIT_FAILURE);
	}

#	ifdef SEEK_END
	/* sanity check: does the file have the correct length ? */
	/* NOTE: this operation is not necessarily supported by all operating systems
	 * but if the preprocessor symbol SEEK_END exists it is assumed to be supported */
	int res = fseek( grid->ioIN, 0, SEEK_END );
	if( res == 0 )
	{
		long End = ftell( grid->ioIN );
		long Expected = grid->nOffset + (grid->nmods+1)*grid->nBlocksize;
		if( End != Expected )
		{
			fprintf( ioQQQ, " InitGrid: Problem performing sanity check for size of binary file.\n" );
			fprintf( ioQQQ, " InitGrid: I expected to find %ld bytes, but actually found %ld bytes.\n",
				 Expected, End );
			fprintf( ioQQQ, " InitGrid: Please recompile the stellar"
				 " atmospheres file with the command: %s.\n", grid->command.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
	}
#	endif
}

/* check whether a binary atmosphere exists and is up-to-date */
STATIC bool lgValidBinFile(const char *binName, process_counter& pc, access_scheme scheme)
{
	DEBUG_ENTRY( "lgValidBinFile()" );

	//
	// this routine is called when either of these two commands is issued:
	// 
	// TABLE STAR AVAIL
	// COMPILE STAR [ additional parameters ]
	//

	stellar_grid grid;
	grid.name = binName;

	if( (grid.ioIN = open_data( grid.name.c_str(), "rb", scheme )) == NULL )
		return false;

	int32 version, mdim, mnam;
	double mesh_elo, mesh_ehi, mesh_res_factor;
	char md5sum[NMD5];
	if( fread( &version, sizeof(version), 1, grid.ioIN ) != 1 ||
	    fread( &mdim, sizeof(mdim), 1, grid.ioIN ) != 1 ||
	    fread( &mnam, sizeof(mnam), 1, grid.ioIN ) != 1 ||
	    fread( &grid.ndim, sizeof(grid.ndim), 1, grid.ioIN ) != 1 ||
	    fread( &grid.npar, sizeof(grid.npar), 1, grid.ioIN ) != 1 ||
	    fread( &grid.nmods, sizeof(grid.nmods), 1, grid.ioIN ) != 1 ||
	    fread( &grid.ngrid, sizeof(grid.ngrid), 1, grid.ioIN ) != 1 ||
	    fread( &grid.nOffset, sizeof(grid.nOffset), 1, grid.ioIN ) != 1 ||
	    fread( &grid.nBlocksize, sizeof(grid.nBlocksize), 1, grid.ioIN ) != 1 ||
	    fread( &mesh_elo, sizeof(mesh_elo), 1, grid.ioIN ) != 1 ||
	    fread( &mesh_ehi, sizeof(mesh_ehi), 1, grid.ioIN ) != 1 ||
	    fread( &mesh_res_factor, sizeof(mesh_res_factor), 1, grid.ioIN ) != 1 ||
	    fread( md5sum, sizeof(md5sum), 1, grid.ioIN ) != 1 )
	{
		return false;
	}

	/* do some sanity checks */
	if( version != VERSION_BIN || mdim != MDIM || mnam != MNAM ||
	    !fp_equal_tol( rfield.emm(), mesh_elo, 1.e-11*rfield.emm() ) ||
	    !fp_equal_tol( rfield.egamry(), mesh_ehi, 1.e-7*rfield.egamry() ) ||
	    !fp_equal( rfield.getResolutionScaleFactor(), mesh_res_factor ) ||
	    strncmp( rfield.mesh_md5sum().c_str(), md5sum, NMD5 ) != 0 )
	{
		return false;
	}

	/* now check the full frequency mesh */
	vector<Energy> anu(rfield.nflux_with_check);
	GetBins( &grid, anu );
	if( !lgValidMesh( anu ) )
	{
		return false;
	}

#	ifdef SEEK_END
	/* sanity check: does the file have the correct length ? */
	/* NOTE: this operation is not necessarily supported by all operating systems
	 * but if the preprocessor symbol SEEK_END exists it is assumed to be supported */
	int res = fseek( grid.ioIN, 0, SEEK_END );
	if( res == 0 )
	{
		long End = ftell( grid.ioIN );
		long Expected = grid.nOffset + (grid.nmods+1)*grid.nBlocksize;
		if( End != Expected )
		{
			return false;
		}
	}
#	endif

	++pc.notProcessed; // the file is up-to-date -> no processing 
	return true;
}

/* check whether a ascii atmosphere file exists and is up-to-date */
STATIC bool lgValidASCIIFile(const char *ascName, access_scheme scheme)
{
	DEBUG_ENTRY( "lgValidASCIIFile()" );

	/* can we read the file? */
	fstream ioIN;
	open_data( ioIN, ascName, mode_r, scheme );
	if( !ioIN.is_open() )
		return false;

	// skip leading comments
	string chLine;
	while( getline( ioIN, chLine ) )
	{
		if( chLine[0] != '#' )
			break;
	}
	if( !ioIN.good() )
		return false;

	/* check version number */
	long version;
	if( sscanf( chLine.c_str(), "%ld", &version ) != 1 ||
	    ( version != VERSION_ASCII && version != VERSION_COSTAR ) )
		return false;

	return true;
}

/* sort CoStar models according to track and index number, store indices in grid->jval[] */
STATIC void InitGridCoStar(stellar_grid *grid) /* the grid parameters */
{
	DEBUG_ENTRY( "InitGridCoStar()" );

	ASSERT( grid->ndim == 2 );
	ASSERT( grid->jlo.size() > 0 );

	swap(grid->jval, grid->jlo);
	grid->jlo.clear();
	grid->jhi.clear();

	/* invalidate contents set by InitGrid first */
	invalidate_array( get_ptr(grid->jval), grid->jval.size()*sizeof(grid->jval[0]) );

	grid->trackLen.resize(grid->nmods);

	vector<long> index(2);
	index[1] = 0;
	while( true )
	{
		bool lgFound;
		index[0] = 0;
		char track = (char)('A'+index[1]);
		do
		{
			lgFound = false;
			for( long i=0; i < grid->nmods; i++ )
			{
				if( grid->telg[i].chGrid == track && grid->telg[i].modid == index[0]+1 )
				{
					grid->jval[JIndex(grid,index)] = i;
					++index[0];
					lgFound = true;
					break;
				}
			}
		}
		while( lgFound );

		if( index[0] == 0 )
			break;

		grid->trackLen[index[1]] = index[0];
		++index[1];
	}

	grid->nTracks = index[1];
}

STATIC void CheckVal(const stellar_grid *grid,
		     double val[], /* val[ndim] */
		     long *nval,
		     long *ndim)
{
	DEBUG_ENTRY( "CheckVal()" );

	if( *ndim == 0 )
		*ndim = (long)grid->ndim;
	if( *ndim == 2 && *nval == 1 && grid->lgIsTeffLoggGrid )
	{
		/* default gravity is maximum gravity */
		val[*nval] = grid->val[1][grid->nval[1]-1];
		++(*nval);
	}
	if( *ndim != (long)grid->ndim )
	{
		fprintf( ioQQQ, " A %ld-dim grid was requested, but a %ld-dim grid was found.\n",
			 *ndim, (long)grid->ndim );
		cdEXIT(EXIT_FAILURE);
	}
	if( *nval < *ndim )
	{
		fprintf( ioQQQ, " A %ld-dim grid was requested, but only %ld parameters were entered.\n",
			 *ndim, *nval );
		cdEXIT(EXIT_FAILURE);
	}
}

STATIC void InterpolateRectGrid(stellar_grid *grid,
				const double val[], /* val[ndim] */
				double *Tlow,
				double *Thigh,
				bool lgTakeLog,
				const realnum Edges[],
				long nEdges)
{
	DEBUG_ENTRY( "InterpolateRectGrid()" );

	/* create some space */
	vector<long> indlo(grid->ndim);
	vector<long> indhi(grid->ndim);
	vector<long> index(grid->ndim);
	vector<double> aval(grid->npar);

	rfield.tNu[rfield.nShape].resize(rfield.nflux_with_check);
	rfield.tslop[rfield.nShape].resize(rfield.nflux_with_check);

	if( grid->lgASCII )
	{
		for( long i=0; i < rfield.nflux_with_check; ++i )
			rfield.tNu[rfield.nShape][i] = (realnum)rfield.anu(i);
	}
	else
	{
		ASSERT( grid->nBlocksize == rfield.nflux_with_check*sizeof(realnum) );

		/* save energy scale for check against code's in conorm */
		GetBins( grid, rfield.tNu[rfield.nShape] );

		/* check that the stored frequency mesh matches what is used in Cloudy */
		ValidateMesh( grid, rfield.tNu[rfield.nShape] );
	}

	if( DEBUGPRT )
	{
		/* check whether the models have the correct effective temperature, for debugging only */
		ValidateGrid( grid, 0.02 );
	}

	/* now generate pointers for models to use */
	for( long nd=0; nd < grid->ndim; nd++ )
	{
		bool lgInvalid;
		FindIndex( grid->val, nd, grid->nval[nd], val[nd], &indlo[nd], &indhi[nd], &lgInvalid );
		if( lgInvalid )
		{
			fprintf( ioQQQ, 
				 " Requested parameter %s = %.2f is not within the range %.2f to %.2f\n",
				 grid->names[nd], val[nd], grid->val[nd][0], grid->val[nd][grid->nval[nd]-1] );
			cdEXIT(EXIT_FAILURE);
		}
	}

	InterpolateModel( grid, val, aval, indlo, indhi, index, grid->ndim, rfield.tslop[rfield.nShape],
			  lgTakeLog, Edges, nEdges );

	/* print the parameters of the interpolated model */
	if( called.lgTalk )
	{
		if( grid->npar == 1 )
			fprintf( ioQQQ, 
				 "                       * #<< FINAL:  %6s = %13.2f"
				 "                                          >>> *\n", 
				 grid->names[0], aval[0] );
		else if( grid->npar == 2 )
			fprintf( ioQQQ, 
				 "                       * #<< FINAL:  %6s = %10.2f"
				 "   %6s = %8.5f                         >>> *\n", 
				 grid->names[0], aval[0], grid->names[1], aval[1] );
		else if( grid->npar == 3 )
			fprintf( ioQQQ, 
				 "                       * #<< FINAL:  %6s = %7.0f"
				 "   %6s = %5.2f   %6s = %5.2f              >>> *\n", 
				 grid->names[0], aval[0], grid->names[1], aval[1],
				 grid->names[2], aval[2] );
		else if( grid->npar >= 4 )
		{
			fprintf( ioQQQ, 
				 "                       * #<< FINAL:  %4s = %7.0f"
				 " %6s = %4.2f %6s = %5.2f %6s = ",
				 grid->names[0], aval[0], grid->names[1], aval[1],
				 grid->names[2], aval[2], grid->names[3] );
			fprintf( ioQQQ, PrintEfmt( "%9.2e", aval[3] ) );
			fprintf( ioQQQ, "  >>> *\n" );
		}
	}	

	for( long i=0; i < rfield.nflux_with_check; i++ )
	{
		if( lgTakeLog )
			rfield.tslop[rfield.nShape][i] = exp10(rfield.tslop[rfield.nShape][i]);
		if( rfield.tslop[rfield.nShape][i] < 1e-37 )
			rfield.tslop[rfield.nShape][i] = 0.;
	}

	if( false )
	{
		FILE *ioBUG = open_data( "interpolated.txt", "w" );
		for( long k=0; k < rfield.nflux_with_check; ++k )
			fprintf( ioBUG, "%e %e\n", rfield.tNu[rfield.nShape][k].Ryd(), rfield.tslop[rfield.nShape][k] );
		fclose( ioBUG );
	}

	if( strcmp( grid->names[0], "Teff" ) == 0 )
	{
		if( ! lgValidModel( rfield.tNu[rfield.nShape], rfield.tslop[rfield.nShape], val[0], 0.10 ) )
			TotalInsanity();
	}

	/* set limits for optimizer */
	vector<realnum> dum;
	SetLimits( grid, val[0], indlo, indhi, NULL, dum, Tlow, Thigh );
}

STATIC void InterpolateModel(stellar_grid *grid,
			     const double val[],
			     vector<double>& aval,
			     const vector<long>& indlo,
			     const vector<long>& indhi,
			     vector<long>& index,
			     long nd,
			     vector<realnum>& flux1,
			     bool lgTakeLog,
			     const realnum Edges[],
			     long nEdges)
{
	DEBUG_ENTRY( "InterpolateModel()" );

	// first determine which models we need to do the interpolation
	InterpolateModel(grid, val, aval, indlo, indhi, index, nd, flux1, IS_COLLECT);

	// emit cautions
	for( map<string,int>::const_iterator p=grid->caution.begin(); p != grid->caution.end(); ++p )
		fprintf( ioQQQ, "%s\n", p->first.c_str() );

	// read the necessary models
	SortUnique( grid->index_list, grid->index_list2 );
	grid->CloudyFlux.alloc(grid->index_list2.size(), rfield.nflux_with_check);
	if( grid->lgASCII )
	{
		if( lgReadAtmosphereTail(grid, Edges, nEdges, grid->index_list2) )
		{
			fprintf( ioQQQ, "Failed to read atmosphere models.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	for( size_t i=0; i < grid->index_list2.size(); ++i )
		GetModel( grid, grid->index_list2[i], &grid->CloudyFlux[i][0], lgVERBOSE, lgTakeLog );

	// and finally carry out the interpolation...
	InterpolateModel(grid, val, aval, indlo, indhi, index, nd, flux1, IS_EXECUTE);
}

STATIC void InterpolateModel(stellar_grid *grid,
			     const double val[],
			     vector<double>& aval,
			     const vector<long>& indlo,
			     const vector<long>& indhi,
			     vector<long>& index,
			     long nd,
			     vector<realnum>& flux1,
			     int stage)
{
	DEBUG_ENTRY( "InterpolateModel()" );

	--nd;

	if( nd < 0 )
	{
		long ind, n = JIndex(grid,index);
		if( stage&IS_FIRST )
			ind = ( grid->jlo[n] >= 0 ) ? grid->jlo[n] : grid->jhi[n];
		else if( stage&IS_SECOND ) 
			ind = ( grid->jhi[n] >= 0 ) ? grid->jhi[n] : grid->jlo[n];
		else if( grid->ndim == 1 )
			/* in this case grid->jlo[n] and grid->jhi[n] should be identical */
			ind = grid->jlo[n];
		else
			TotalInsanity();

		if( ind < 0 )
		{
			fprintf( ioQQQ, " The requested interpolation could not be completed, sorry.\n" );
			fprintf( ioQQQ, " No suitable match was found for a model with" );
			for( long i=0; i < grid->ndim; i++ )
				fprintf( ioQQQ, " %s=%.6g ", grid->names[i], grid->val[i][index[i]] );
			fprintf( ioQQQ, "\n" );
			cdEXIT(EXIT_FAILURE);
		}

		for( long i=0; i < grid->npar; i++ )
			aval[i] = grid->telg[ind].par[i];

		if( stage&IS_COLLECT )
		{
			for( long i=0; i < grid->ndim && called.lgTalk; i++ )
			{
				if( !fp_equal(grid->val[i][index[i]],aval[i],10) )
				{
					ostringstream oss;
					oss << " No exact match was found for a model with";
					for( long j=0; j < grid->ndim; j++ )
						oss << " " << grid->names[j] << "=" << grid->val[j][index[j]] << " ";
					oss << "- using model " << ind+1 << " instead.";
					// use a map to weed out duplicate cautions
					grid->caution[oss.str()] = 1;
					break;
				}
			}

			grid->index_list.push_back(ind);
		}
		else if( stage&IS_EXECUTE )
		{
			size_t i = 0;
			while( i < grid->index_list2.size() && grid->index_list2[i] != ind )
				++i;
			ASSERT( i < grid->index_list2.size() && grid->CloudyFlux.size() > 0 );

			for( long j=0; j < rfield.nflux_with_check; ++j )
				flux1[j] = grid->CloudyFlux[i][j];
		}
		else
			TotalInsanity();
	}
	else
	{
#		if !defined NDEBUG
		const realnum SECURE = 10.f*FLT_EPSILON;
#		endif

		/* Interpolation is carried out first in the parameter with nd == 0 (usually
		 * Teff), then the parameter with nd == 1 (usually log(g)), etc. One or two
		 * atmosphere models are read depending on whether the parameter was matched
		 * exactly or not. If needed, logarithmic interpolation is done.
		 */

		index[nd] = indlo[nd];
		int next_stage = ( nd == 1 ) ? stage|IS_FIRST : stage;
		InterpolateModel( grid, val, aval, indlo, indhi, index, nd, flux1, next_stage );

		vector<realnum> flux2(rfield.nflux_with_check);
		vector<double> aval2(grid->npar);

		index[nd] = indhi[nd];
		next_stage = ( nd == 1 ) ? stage|IS_SECOND : stage;
		InterpolateModel( grid, val, aval2, indlo, indhi, index, nd, flux2, next_stage );

		if( (stage&IS_EXECUTE) && !fp_equal(aval2[nd],aval[nd],10) )
		{
			double fr1 = (aval2[nd]-val[nd])/(aval2[nd]-aval[nd]);
			/* when interpolating in log(g) it can happen that fr1 is outside the range 0 .. 1
			 * this can be the case when the requested log(g) was not present in the grid
			 * and it had to be approximated by another model. In this case do not extrapolate */
			if( nd == 1 )
				fr1 = MIN2( MAX2( fr1, 0. ), 1. );
			double fr2 = 1. - fr1;

			ASSERT( 0.-SECURE <= fr1 && fr1 <= 1.+SECURE );

			if( DEBUGPRT )
				fprintf( ioQQQ, "interpolation nd=%ld fr1=%g\n", nd, fr1 );

			/* special treatment for high-temperature Rauch models */
			double fc1 = 0., fc2 = 0.;
			if( nd == 0 && strcmp( grid->names[nd], "Teff" ) == 0 )
			{
				/* The following is an approximate scaling to use for the range of 
				 * temperatures above 200000 K in the H-Ca Rauch models where the
				 * temperature steps are large and thus the interpolations are over
				 * large ranges.  For the lower temperatures I assume that there is
				 * no need for this.
				 *
				 * It should be remembered that this interpolation is not exact, and 
				 * the possible error at high temperatures might be large enough to 
				 * matter. (Kevin Volk)
				 */
				fc1 = ( val[nd] > 200000. ) ? log10(val[nd]/grid->val[nd][indlo[nd]])*4. : 0.;
				fc2 = ( val[nd] > 200000. ) ? log10(val[nd]/grid->val[nd][indhi[nd]])*4. : 0.;
			}

			for( long i=0; i < rfield.nflux_with_check; ++i )
				flux1[i] = (realnum)(fr1*(flux1[i]+fc1) + fr2*(flux2[i]+fc2));

			for( long i=0; i < grid->npar; i++ )
				aval[i] = fr1*aval[i] + fr2*aval2[i];
		}
	}
}

STATIC void InterpolateModelCoStar(const stellar_grid *grid,
				   const double val[],
				   vector<double>& aval,
				   const long indlo[],
				   const long indhi[],
				   vector<long>& index,
				   long nd,
				   long off,
				   vector<realnum>& flux1)
{
	DEBUG_ENTRY( "InterpolateModelCoStar()" );

	if( nd == 2 )
	{
		long ind = ( index[1] == 0 ) ? indlo[index[0]] : indhi[index[0]];

		size_t i = 0;
		while( i < grid->index_list2.size() && grid->index_list2[i] != ind )
			++i;
		ASSERT( i < grid->index_list2.size() && grid->CloudyFlux.size() > 0 );

		for( long j=0; j < rfield.nflux_with_check; ++j )
			flux1[j] = grid->CloudyFlux[i][j];

		for( long j=0; j < grid->npar; j++ )
			aval[j] = grid->telg[ind].par[j];
	}
	else
	{
#		if !defined NDEBUG
		const realnum SECURE = 10.f*FLT_EPSILON;
#		endif

		/* Interpolation is carried out first along evolutionary tracks, then
		 * in between evolutionary tracks. Between 1 and 4 atmosphere models are read
		 * depending on whether the parameter/track was matched exactly or not.
		 */

		index[nd] = 0;
		InterpolateModelCoStar( grid, val, aval, indlo, indhi, index, nd+1, off, flux1 );

		bool lgSkip = ( nd == 1 ) ?  ( indhi[index[0]] == indlo[index[0]] ) :
			( indlo[0] == indlo[1] && indhi[0] == indhi[1] );

		if( ! lgSkip )
		{
			vector<realnum> flux2(rfield.nflux_with_check);
			vector<double> aval2(grid->npar);

			index[nd] = 1;
			InterpolateModelCoStar( grid, val, aval2, indlo, indhi, index, nd+1, off, flux2 );

			double fr1 = (aval2[nd+off]-val[nd])/(aval2[nd+off]-aval[nd+off]);
			double fr2 = 1. - fr1;

			if( DEBUGPRT )
				fprintf( ioQQQ, "interpolation nd=%ld fr1=%g\n", nd, fr1 );

			ASSERT( 0.-SECURE <= fr1 && fr1 <= 1.+SECURE );

			for( long i=0; i < rfield.nflux_with_check; ++i )
				flux1[i] = (realnum)(fr1*flux1[i] + fr2*flux2[i]);

			for( long i=0; i < grid->npar; i++ )
				aval[i] = fr1*aval[i] + fr2*aval2[i];
		}
	}
}

template<class T>
void SortUnique(vector<T>& in,
		vector<T>& out)
{
	DEBUG_ENTRY( "SortUnique()" );

	// first sort array "in" and then create a copy in "out" removing duplicate entries
	sort( in.begin(), in.end() );
	out.clear();
	T lastval = in[0];
	out.push_back( lastval );
	for( size_t i=1; i < in.size(); ++i )
		if( in[i] != lastval )
		{
			lastval = in[i];
			out.push_back( lastval );
		}
}

STATIC void GetBins(const stellar_grid *grid,
		    vector<Energy>& ener)
{
	DEBUG_ENTRY( "GetBins()" );

	ASSERT( grid->nBlocksize == rfield.nflux_with_check*sizeof(realnum) );

	/* skip over ind stars */
	/* >>chng 01 oct 18, add nOffset */
	if( fseek( grid->ioIN, (long)(grid->nOffset), SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " Error finding atmosphere frequency bins\n");
		cdEXIT(EXIT_FAILURE);
	}

	vector<realnum> data(rfield.nflux_with_check);
	if( fread( get_ptr(data), 1, grid->nBlocksize, grid->ioIN ) != grid->nBlocksize )
	{
		fprintf( ioQQQ, " Error reading atmosphere frequency bins\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( ener.size() != size_t(rfield.nflux_with_check) )
		ener.resize(rfield.nflux_with_check);

	for( long i=0; i < rfield.nflux_with_check; ++i )
		ener[i].set(data[i]);
}

STATIC void GetModel(const stellar_grid *grid,
		     long ind,
		     realnum *flux,
		     bool lgTalk,
		     bool lgTakeLog)
{
	DEBUG_ENTRY( "GetModel()" );

	/* add 1 to account for frequency grid that is stored in front of all the atmospheres */
	ind++;

	/* make sure ident is exactly 12 characters long, otherwise output won't fit */
	ASSERT( grid->ident.length() == 12 );
	/* ind == 0 is the frequency grid, ind == 1 .. nmods are the atmosphere models */
	ASSERT( ind >= 0 && ind <= grid->nmods );

	if( !grid->lgASCII )
	{
		/* skip over ind stars */
		/* >>chng 01 oct 18, add nOffset */
		if( fseek( grid->ioIN, (long)(ind*grid->nBlocksize+grid->nOffset), SEEK_SET ) != 0 )
		{
			fprintf( ioQQQ, " Error seeking atmosphere %ld\n", ind );
			cdEXIT(EXIT_FAILURE);
		}

		if( fread( get_ptr(flux), 1, grid->nBlocksize, grid->ioIN ) != grid->nBlocksize )
		{
			fprintf( ioQQQ, " Error trying to read atmosphere %ld\n", ind );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* print the parameters of the atmosphere model */
	if( called.lgTalk && lgTalk )
	{
		/* ind-1 below since telg doesn't have the entry for the frequency grid */
		if( grid->npar == 1 )
			fprintf( ioQQQ, 
				 "                       * #<< %s model%5ld read.  "
				 "  %6s = %13.2f                 >>> *\n", 
				 grid->ident.c_str(), ind, grid->names[0], grid->telg[ind-1].par[0] );
		else if( grid->npar == 2 )
			fprintf( ioQQQ, 
				 "                       * #<< %s model%5ld read.  "
				 "  %6s = %10.2f %6s = %8.5f  >>> *\n", 
				 grid->ident.c_str(), ind, grid->names[0], grid->telg[ind-1].par[0],
				 grid->names[1], grid->telg[ind-1].par[1] );
		else if( grid->npar == 3 )
			fprintf( ioQQQ, 
				 "                       * #<< %s model%5ld read. "
				 " %6s=%7.0f %6s=%5.2f %6s=%5.2f >>> *\n", 
				 grid->ident.c_str(), ind, grid->names[0], grid->telg[ind-1].par[0],
				 grid->names[1], grid->telg[ind-1].par[1],
				 grid->names[2], grid->telg[ind-1].par[2] );
		else if( grid->npar >= 4 )
		{
			fprintf( ioQQQ, 
				 "                       * #< %s mdl%4ld"
				 " %4s=%5.0f %6s=%4.2f %6s=%5.2f %6s=",
				 grid->ident.c_str(), ind, grid->names[0], grid->telg[ind-1].par[0],
				 grid->names[1], grid->telg[ind-1].par[1],
				 grid->names[2], grid->telg[ind-1].par[2], grid->names[3] );
			fprintf( ioQQQ, PrintEfmt( "%9.2e", grid->telg[ind-1].par[3] ) );
			fprintf( ioQQQ, " >> *\n" ); 
		}
	}	

	if( lgTakeLog )
	{
		/* convert to logs since we will interpolate in log flux */
		for( long i=0; i < rfield.nflux_with_check; ++i )
		{
			// the keyword volatile is needed to work around a
			// compiler bug in g++ versions 4.7.0 and later
			// see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=65425
		        volatile double help = flux[i];
			if( help > 0. )
				help = log10(help);
			else
				help = -99999.;
			flux[i] = realnum(help);
		}
	}
}

STATIC void SetLimits(const stellar_grid *grid,
		      double val,
		      const vector<long>& indlo,
		      const vector<long>& indhi,
		      const long useTr[],
		      const vector<realnum>& ValTr,
		      double *loLim,
		      double *hiLim)
{
	DEBUG_ENTRY( "SetLimits()" );

	if( optimize.lgVarOn )
	{
		const double SECURE = (1. + 20.*(double)FLT_EPSILON);

		int ptr0, ptr1;
		vector<long> index(MDIM);
		*loLim = +DBL_MAX;
		*hiLim = -DBL_MAX;

		switch( grid->imode )
		{
		case IM_RECT_GRID:
			*loLim = -DBL_MAX;
			*hiLim = +DBL_MAX;
			SetLimitsSub( grid, val, indlo, indhi, index, grid->ndim, loLim, hiLim );
			break;
		case IM_COSTAR_TEFF_MODID:
		case IM_COSTAR_TEFF_LOGG:
		case IM_COSTAR_MZAMS_AGE:
			for( long j=0; j < grid->nTracks; j++ )
			{
				if( ValTr[j] != -FLT_MAX )
				{
					/* M_ZAMS is already logarithm, Teff is linear */
					double temp = ( grid->imode == IM_COSTAR_MZAMS_AGE ) ?
						exp10((double)ValTr[j]) : ValTr[j];
					*loLim = MIN2(*loLim,temp);
					*hiLim = MAX2(*hiLim,temp);
				}
			}
			break;
		case IM_COSTAR_AGE_MZAMS:
			index[0] = 0;
			index[1] = useTr[0];
			ptr0 = grid->jval[JIndex(grid,index)];
			index[1] = useTr[1];
			ptr1 = grid->jval[JIndex(grid,index)];
			*loLim = MAX2(grid->telg[ptr0].par[3],grid->telg[ptr1].par[3]);
			if( DEBUGPRT )
			{
				printf( "set limit 0: (models %d, %d) %f %f\n",
					ptr0+1, ptr1+1, grid->telg[ptr0].par[3], grid->telg[ptr1].par[3] );
			}
			index[0] = grid->trackLen[useTr[0]]-1;
			index[1] = useTr[0];
			ptr0 = grid->jval[JIndex(grid,index)];
			index[0] = grid->trackLen[useTr[1]]-1;
			index[1] = useTr[1];
			ptr1 = grid->jval[JIndex(grid,index)];
			*hiLim = MIN2(grid->telg[ptr0].par[3],grid->telg[ptr1].par[3]);
			if( DEBUGPRT )
			{
				printf( "set limit 1: (models %d, %d) %f %f\n",
					ptr0+1, ptr1+1, grid->telg[ptr0].par[3], grid->telg[ptr1].par[3] );
			}
			break;
		default:
			fprintf( ioQQQ, " SetLimits called with insane value for imode: %d.\n", grid->imode );
			cdEXIT(EXIT_FAILURE);
		}

		ASSERT( fabs(*loLim) < DBL_MAX && fabs(*hiLim) < DBL_MAX );

		/* check sanity of optimization limits */
		if( *hiLim <= *loLim )
		{
			fprintf( ioQQQ, " no room to optimize: lower limit %.4f, upper limit %.4f.\n",
				 *loLim,*hiLim );
			cdEXIT(EXIT_FAILURE);
		}

		/* make a bit of room for round-off errors */
		*loLim *= SECURE;
		*hiLim /= SECURE;

		if( DEBUGPRT )
			printf("set limits: %g %g\n",*loLim,*hiLim);
	}
	else
	{
		*loLim = 0.;
		*hiLim = 0.;
	}
}

STATIC void SetLimitsSub(const stellar_grid *grid,
			 double val,
			 const vector<long>& indlo,
			 const vector<long>& indhi,
			 vector<long>& index,
			 long nd,
			 double *loLim,
			 double *hiLim)
{
	DEBUG_ENTRY( "SetLimitsSub()" );

	--nd;

	if( nd < 1 )
	{
		double loLoc = +DBL_MAX;
		double hiLoc = -DBL_MAX;

		for( index[0]=0; index[0] < grid->nval[0]; ++index[0] )
		{
			/* grid->val[0][i] is the array of Par0 values (Teff/Age/...) in the
			 * grid, which it is sorted in strict monotonically increasing order.
			 * This routine searches for the largest range [loLoc,hiLoc] in Par0
			 * such that loLoc <= val <= hiLoc, and at least one model exists for
			 * each Par0 value in this range. This assures that interpolation is
			 * safe and the optimizer will not trip... */
			long n = JIndex(grid,index);
			if( grid->jlo[n] < 0 && grid->jhi[n] < 0 )
			{
				/* there are no models with this value of Par0 */
				/* this value of Par0 should be outside of allowed range */
				if( grid->val[0][index[0]] < val )
					loLoc = DBL_MAX;
				/* this is beyond the legal range, so terminate the search */
				if( grid->val[0][index[0]] > val )
					break;
			}
			else
			{
				/* there are models with this value of Par0 */
				/* update range to include this value of Par0 */
				if( grid->val[0][index[0]] <= val )
				{
					/* remember lowest legal value of loLoc
					 * -> only update if previous value was illegal */
					if( loLoc == DBL_MAX )
						loLoc = grid->val[0][index[0]];
				}
				if( grid->val[0][index[0]] >= val )
				{
					/* remember highest legal value of hiLoc
					 * -> always update */
					hiLoc = grid->val[0][index[0]];
				}
			}
		}

		ASSERT( fabs(loLoc) < DBL_MAX && fabs(hiLoc) < DBL_MAX && loLoc <= hiLoc );

		*loLim = MAX2(*loLim,loLoc);
		*hiLim = MIN2(*hiLim,hiLoc);
	}
	else
	{
		index[nd] = indlo[nd];
		SetLimitsSub( grid, val, indlo, indhi, index, nd, loLim, hiLim );

		if( indhi[nd] != indlo[nd] )
		{
			index[nd] = indhi[nd];
			SetLimitsSub( grid, val, indlo, indhi, index, nd, loLim, hiLim );
		}
	}
}

STATIC void InitIndexArrays(stellar_grid *grid,
			    bool lgList)
{
	DEBUG_ENTRY( "InitIndexArrays()" );

	ASSERT( grid->telg.size() > 0 );
	ASSERT( grid->nmods > 0 );

	long jsize = 1;

	/* this loop creates a list of all unique model parameter values in increasing order */
	for( long nd=0; nd < grid->ndim; nd++ )
	{
		double pval = grid->telg[0].par[nd];
		grid->val[nd][0] = pval;
		grid->nval[nd] = 1;

		for( long i=1; i < grid->nmods; i++ )
		{
			bool lgOutOfRange;
			long i1, i2;

			pval = grid->telg[i].par[nd];
			FindIndex( grid->val, nd, grid->nval[nd], pval, &i1, &i2, &lgOutOfRange );
			/* if i1 < i2, the new parameter value was not present yet and
			 * it needs to be inserted in between i1 and i2 --> first move
			 * all entries from i2 to grid->nval[nd]-1 one slot upward and
			 * then insert the new value at i2; this also works correctly
			 * if lgOutOfRange is set, hence no special check is needed */ 
			if( i1 < i2 )
			{
				/* val[nd] has grid->nmods entries, so cannot overflow */
				for( long j = grid->nval[nd]-1; j >= i2; j-- )
					grid->val[nd][j+1] = grid->val[nd][j];
				grid->val[nd][i2] = pval;
				grid->nval[nd]++;
			}
		}

		jsize *= grid->nval[nd];

		if( DEBUGPRT )
		{
			printf( "%s[%ld]:", grid->names[nd], grid->nval[nd] );
			for( long i=0; i < grid->nval[nd]; i++ )
				printf( " %g", grid->val[nd][i] );
			printf( "\n" );
		}
	}

	vector<long> index(grid->ndim);
	vector<double> val(grid->ndim);

	grid->jlo.resize(jsize);
	grid->jhi.resize(jsize);

	/* set up square array of model indices; this will be used to
	 * choose the correct models for the interpolation process */
	FillJ( grid, index, val, grid->ndim, lgList );

	if( lgList )
		cdEXIT(EXIT_SUCCESS);
}

STATIC void FillJ(stellar_grid *grid,
		  vector<long>& index, /* index[grid->ndim] */
		  vector<double>& val, /* val[grid->ndim] */
		  long nd,
		  bool lgList)
{
	DEBUG_ENTRY( "FillJ()" );

	--nd;

	if( nd < 0 )
	{
		long n = JIndex(grid,index);
		SearchModel( grid->telg, grid->lgIsTeffLoggGrid, grid->nmods, val, grid->ndim,
			     &grid->jlo[n], &grid->jhi[n] );
	}
	else
	{
		for( index[nd]=0; index[nd] < grid->nval[nd]; index[nd]++ )
		{
			val[nd] = grid->val[nd][index[nd]];
			FillJ( grid, index, val, nd, lgList );
		}
	}

	if( lgList && nd == MIN2(grid->ndim-1,1) )
	{
		fprintf( ioQQQ, "\n" );
		if( grid->ndim > 2 )
		{
			fprintf( ioQQQ, "subgrid for" );
			for( long n = nd+1; n < grid->ndim; n++ )
				fprintf( ioQQQ, " %s=%g", grid->names[n], val[n] );
			fprintf( ioQQQ, ":\n\n" );
		}
		if( grid->ndim > 1 )
		{
			fprintf( ioQQQ, "%6.6s\\%6.6s |", grid->names[0], grid->names[1] );
			for( long n = 0; n < grid->nval[1]; n++ )
				fprintf( ioQQQ, " %9.3g", grid->val[1][n] );
			fprintf( ioQQQ, "\n" );
			fprintf( ioQQQ, "--------------|" );
			for( long n = 0; n < grid->nval[1]; n++ )
				fprintf( ioQQQ, "----------" );
		}
		else
		{
			fprintf( ioQQQ, "%13.13s |\n", grid->names[0] );
			fprintf( ioQQQ, "--------------|----------" );
		}
		fprintf( ioQQQ, "\n" );
		for( index[0]=0; index[0] < grid->nval[0]; index[0]++ )
		{
			fprintf( ioQQQ, "%13.7g |", grid->val[0][index[0]] );
			if( grid->ndim > 1 )
			{
				for( index[1]=0; index[1] < grid->nval[1]; index[1]++ )
					if( grid->jlo[JIndex(grid,index)] == grid->jhi[JIndex(grid,index)] &&
					    grid->jlo[JIndex(grid,index)] >= 0 )
						fprintf( ioQQQ, " %9ld", grid->jlo[JIndex(grid,index)]+1 );
					else
						fprintf( ioQQQ, "        --" );
			}
			else
			{
				fprintf( ioQQQ, " %9ld", grid->jlo[JIndex(grid,index)]+1 );
			}
			fprintf( ioQQQ, "\n" );
		}
		fprintf( ioQQQ, "\n" );
	}
}

STATIC long JIndex(const stellar_grid *grid,
		   const vector<long>& index) /* index[grid->ndim] */
{
	DEBUG_ENTRY( "JIndex()" );

	long ind = 0;
	long mul = 1;
	for( long i=0; i < grid->ndim; i++ )
	{
		ind += index[i]*mul;
		mul *= grid->nval[i];
	}
	return ind;
}

STATIC void SearchModel(const vector<mpp>& telg, /* telg[nmods] */
			bool lgIsTeffLoggGrid,
			long nmods,
			const vector<double>& val, /* val[ndim] */
			long ndim,
			long *index_low,
			long *index_high)
{
	DEBUG_ENTRY( "SearchModel()" );

	/* given values for the model parameters, this routine searches for the atmosphere
	 * model that is the best match. If all parameters can be matched simultaneously the
	 * choice is obvious, but this cannot always be achieved (typically for high Teff, the
	 * low log(g) models will be missing). If lgIsTeffLoggGrid is true, the rule is that
	 * all parameters except log(g) must always be matched (such a model is not always
	 * guaranteed to exist). If all requested parameters can be matched exactly, both
	 * index_low and index_high will point to that model. If all parameters except log(g)
	 * can be matched exactly, it will return the model with the lowest log(g) value larger
	 * than the requested value in index_high, and the model with the highest log(g) value
	 * lower than the requested value in index_low. If either requirement cannot be
	 * fulfilled, -2 will be returned. When lgIsTeffLoggGrid is false, all parameters must
	 * be matched and both index_low and index_high will point to that model. If no such
	 * model can be found, -2 will be returned. */

	*index_low = *index_high = -2;
	double alogg_low = -DBL_MAX, alogg_high = DBL_MAX;
	for( long i=0; i < nmods; i++ )
	{
		bool lgNext = false;
		/* ignore models with different parameters */
		for( long nd=0; nd < ndim; nd++ )
		{
			if( nd != 1 && !fp_equal(telg[i].par[nd],val[nd],10) )
			{
				lgNext = true;
				break;
			}
		}
		if( lgNext )
			continue;

		/* an exact match is found */
		if( ndim == 1 || fp_equal(telg[i].par[1],val[1],10) )
		{
			*index_low = i;
			*index_high = i;
			return;
		}
		if( lgIsTeffLoggGrid )
		{
			/* keep a record of the highest log(g) model smaller than alogg */
			if( telg[i].par[1] < val[1] && telg[i].par[1] > alogg_low )
			{
				*index_low = i;
				alogg_low = telg[i].par[1];
			}
			/* also keep a record of the lowest log(g) model greater than alogg */
			if( telg[i].par[1] > val[1] && telg[i].par[1] < alogg_high )
			{
				*index_high = i;
				alogg_high = telg[i].par[1];
			}
		}
	}
}

STATIC void FindIndex(const multi_arr<double,2>& xval, /* xval[NVAL] */
		      long nd,
		      long NVAL,
		      double x,
		      long *ind1,
		      long *ind2,
		      bool *lgInvalid)
{
	DEBUG_ENTRY( "FindIndex()" );

	/* this routine searches for indices ind1, ind2 such that
	 *   xval[nd][ind1] < x < xval[nd][ind2]
	 * if x is equal to one of the values in xval, then
	 *   ind1 == ind2  and  xval[nd][ind1] == x
	 *
	 * if x is outside the range xval[nd][0] ... xval[nd][NVAL-1]
	 * then lgInvalid will be set to true
	 *
	 * NB NB -- this routine implicitly assumes that xval is
	 *          strictly monotonically increasing!
	 */

	ASSERT( NVAL > 0 );

	/* is x outside of range xval[nd][0] ... xval[nd][NVAL-1]? */
	bool lgOutLo = ( x-xval[nd][0] < -10.*DBL_EPSILON*fabs(xval[nd][0]) );
	bool lgOutHi = ( x-xval[nd][NVAL-1] > 10.*DBL_EPSILON*fabs(xval[nd][NVAL-1]) );

	if( lgOutLo || lgOutHi )
	{
		/* pretend there are two fictitious array elements
		 *   xval[nd][-1] = -Inf  and  xval[nd][NVAL] = +Inf,
		 * and return ind1 and ind2 accordingly. This behavior
		 * is needed for InitIndexArrays() to work correctly */
		*ind1 = lgOutLo ? -1 : NVAL-1;
		*ind2 = lgOutLo ?  0 : NVAL;
		*lgInvalid = true;
		return;
	}

	*lgInvalid = false;

	/* there are more efficient ways of doing this, e.g. a binary search.
	 * However, the xval arrays typically only have 1 or 2 dozen elements,
	 * so the overhead is negligible and the clarity of this code is preferred */

	/* first look for an "exact" match */
	for( long i=0; i < NVAL; i++ )
	{
		if( fp_equal(xval[nd][i],x,10) )
		{
			*ind1 = i;
			*ind2 = i;
			return;
		}
	}

	/* no match was found -> bracket the x value */
	for( long i=0; i < NVAL-1; i++ )
	{
		if( xval[nd][i] < x && x < xval[nd][i+1] )
		{
			*ind1 = i;
			*ind2 = i+1;
			return;
		}
	}

	/* this should never be reached ! */
	TotalInsanity();
}

STATIC bool lgFileReadable(const char *chFnam, process_counter& pc, access_scheme scheme)
{
	DEBUG_ENTRY( "lgFileReadable()" );

	string fname = chFnam;
	bool lgASCIIfile = ( fname.substr(max(fname.length(),6)-6) == ".ascii" );
	FILE *ioIN = open_data( chFnam, "r", scheme );
	if( ioIN != NULL )
	{
		fclose( ioIN );
		if( lgASCIIfile )
			++pc.nFound;
		return true;
	}
	else
	{
		return false;
	}
}

/* check that the stored frequency mesh matches what is used in Cloudy */
STATIC void ValidateMesh(const stellar_grid *grid,
			 const vector<Energy>& anu)
{
	DEBUG_ENTRY( "ValidateMesh()" );

	if( !lgValidMesh( anu ) ) 
	{
		fprintf( ioQQQ, " ValidateMesh: the compiled atmospheres file is produced"
			 " with an incompatible version of Cloudy.\n" );
		fprintf( ioQQQ, " ValidateMesh: Please recompile the stellar"
			 " atmospheres file with the command: %s.\n", grid->command.c_str() );
		cdEXIT(EXIT_FAILURE);
	}
}

STATIC bool lgValidMesh(const vector<Energy>& anu)
{
	DEBUG_ENTRY( "lgValidMesh()" );

	for( long i=0; i < rfield.nflux_with_check; ++i )
	{
		// check with 32-bit FP precision since the binary file stores realnums
		// also use 32-bit precision with -DFLT_IS_DBL since the fundamental constants
		// may have changed since the file was written...
		if( !fp_equal_tol( anu[i].Ryd(), rfield.anu(i), 3.*double(FLT_EPSILON)*anu[i].Ryd() ) )
			return false;
	}
	return true;
}

/*ValidateGrid: check each model in the grid to see if it has the correct Teff */
STATIC void ValidateGrid(const stellar_grid *grid,
			 double toler)
{
	DEBUG_ENTRY( "ValidateGrid()" );

	if( strcmp( grid->names[0], "Teff" ) != 0 )
		return;

	vector<Energy> anu(rfield.nflux_with_check);
	vector<realnum> flux(rfield.nflux_with_check);

	GetBins( grid, anu );

	for( long i=0; i < grid->nmods; i++ ) 
	{
		fprintf( ioQQQ, "testing model %ld ", i+1 );
		for( long nd=0; nd < grid->npar; nd++ )
			fprintf( ioQQQ, " %s %g", grid->names[nd], grid->telg[i].par[nd] );

		GetModel( grid, i, get_ptr(flux), lgSILENT, lgLINEAR );

		if( lgValidModel( anu, flux, grid->telg[i].par[0], toler ) )
			fprintf( ioQQQ, "   OK\n" );
	}
}

STATIC bool lgValidModel(const vector<Energy>& anu,
			 const vector<realnum>& flux,
			 double Teff,
			 double toler)
{
	DEBUG_ENTRY( "lgValidModel()" );

	ASSERT( Teff > 0. );

	double lumi = 0.;
	/* rebinned models are in cgs F_nu units */
	for( long k=1; k < rfield.nflux_with_check; k++ )
		lumi += (anu[k].Ryd() - anu[k-1].Ryd())*(flux[k] + flux[k-1])/2.;

	/* now convert luminosity to effective temperature */
	double chk = powpq(lumi*FR1RYD/STEFAN_BOLTZ,1,4);
	/* the allowed tolerance is set by the caller in toler */
	bool lgPassed = true;
	if( fabs(Teff - chk) > toler*Teff ) {
		fprintf( ioQQQ, "\n*** WARNING, Teff discrepancy for this model, expected Teff %.2f, ", Teff);
		fprintf( ioQQQ, "integration yielded Teff %.2f, delta %.2f%%\n", chk, (chk/Teff-1.)*100. );
		lgPassed = false;
	}
	return lgPassed;
}

/*RebinAtmosphere: generic routine for rebinning atmospheres onto Cloudy grid */
STATIC void RebinAtmosphere(const vector<realnum>& StarEner, /* StarEner[nCont], the freq grid for the model, in Ryd*/
			    const vector<realnum>& StarFlux, /* StarFlux[nCont], the original model flux */
			    long nCont,
			    const realnum Edges[],   /* Edges[nEdges], energies of the edges */
			    long nEdges,             /* the number of bound-free continuum edges in AbsorbEdge */
			    realnum CloudyFlux[])    /* CloudyFlux[], the model flux on the cloudy grid */
{
	DEBUG_ENTRY( "RebinAtmosphere()" );

	/* cut off that part of the Wien tail that evaluated to zero */
	for( long j=nCont-1; j >= 0; j-- )
	{
		if( StarFlux[j] > 0.f )
		{
			nCont = j+1;
			break;
		}
	}
	ASSERT( nCont > 0 );

	vector<long> EdgeInd;
	vector<realnum> EdgeLow, EdgeHigh;
	for( long j=0; j < nEdges; j++ )
	{
		long ind = RebinFind(get_ptr(StarEner), nCont, Edges[j]);

		if( ind >= 1 && ind+2 < nCont )
		{
			EdgeInd.push_back( ind );
			EdgeLow.push_back( StarEner[ind] );
			EdgeHigh.push_back( StarEner[ind+1] );
		}
		else
		{
			EdgeInd.push_back( -1 );
			EdgeLow.push_back( 0. );
			EdgeHigh.push_back( 0. );
		}
	}

	vector<realnum> StarPower(nCont-1);

	for( long j=0; j < nCont-1; j++ )
	{
		/* >>chng 05 nov 22, add sanity check to prevent invalid fp operations */
		ASSERT( StarEner[j+1] > StarEner[j] );

		/* >>chng 06 aug 11, on some systems (e.g., macbook pro) y/x can get evaluated as y*(1/x);
		 * this causes overflows if x is a denormalized number, hence we force a cast to double, PvH */
		double ratio_x = (double)StarEner[j+1]/(double)StarEner[j];
		if( StarFlux[j] == 0.f )
			StarPower[j] = FLT_MAX;
		else if( StarFlux[j+1] == 0.f )
			StarPower[j] = -FLT_MAX;
		else
		{
			double ratio_y = (double)StarFlux[j+1]/(double)StarFlux[j];
			StarPower[j] = (realnum)(log(ratio_y)/log(ratio_x));
		}
	}

	for( long j=0; j < rfield.nflux_with_check; j++ )
	{
		/* >>chng 00 aug 14, take special care not to interpolate over major edges,
		 * the region in between EdgeLow and EdgeHigh should be avoided,
		 * the spectrum is extremely steep there, leading to significant roundoff error, PvH */
		bool lgDone = false;
		for( long k=0; k < nEdges; k++ )
		{
			if( EdgeInd[k] > 0 && rfield.anumax(j) > EdgeLow[k] && rfield.anumin(j) < EdgeHigh[k] )
			{
				long ipLo;
				if( rfield.anu(j) < Edges[k] )
					ipLo = EdgeInd[k]-1; // extrapolate from lower cell
				else
					ipLo = EdgeInd[k]+1; // extrapolate from higher cell
				// the cell with the edge should have a steeper gradient than the adjacent cell
				if( fabs(StarPower[ipLo]) < fabs(StarPower[EdgeInd[k]]) &&
				    fabs(StarPower[ipLo]) != FLT_MAX )
				{
					CloudyFlux[j] = StarFlux[ipLo]*
						pow(rfield.anu(j)/StarEner[ipLo],(double)StarPower[ipLo]);
					lgDone = true;
				}
				break;
			}
		}

		/* default case when we are not close to an edge */
		if( !lgDone )
			CloudyFlux[j] = RebinSingleCell(j, get_ptr(StarEner), get_ptr(StarFlux), StarPower, nCont);
	}
}

STATIC realnum RebinSingleCell(long j,
			       const realnum StarEner[],         /* StarEner[nCont] */
			       const realnum StarFlux[],         /* StarFlux[nCont] */
			       const vector<realnum>& StarPower, /* StarPower[nCont-1] */
			       long nCont)
{
	DEBUG_ENTRY( "RebinSingleCell()" );

	double BinLow = rfield.anumin(j);
	double BinHigh = rfield.anumax(j);
	double anu = rfield.anu(j);
	/* >>chng 05 nov 22, reduce widflx if cell sticks out above highest frequency in model, PvH */
	double widflx = MIN2(BinHigh,StarEner[nCont-1])-BinLow;
	double retval;

	if( BinLow < StarEner[0] )
	{
		/* this is case where Cloudy's continuum is below stellar continuum,
		 * (at least for part of the cell), so we do Rayleigh Jeans extrapolation */
		retval = (realnum)(StarFlux[0]*pow2(anu/StarEner[0]));
	}
	else if( BinLow > StarEner[nCont-1] )
	{
		/* case where cloudy continuum is entirely above highest stellar point */
		retval = 0.0e00;
	}
	else
	{
		/* now go through stellar continuum to find bins corresponding to
		 * this cloudy cell, stellar continuum defined through nCont cells */
		long ipLo = RebinFind(StarEner,nCont,BinLow);
		long ipHi = RebinFind(StarEner,nCont,BinHigh);
		/* sanity check */
		ASSERT( ipLo >= 0 && ipLo < nCont-1 && ipHi >= ipLo );

		/* when either Fnu(i) or Fnu(i+1) is zero we cannot use logarithmic integration.
		 * in this case we use linear integration using the following assumptions:
		 *
		 * if Fnu(i+1) is zero we assume that the flux is constant at Fnu(i) up to the
		 * midpoint between nu(i) and nu(i+1) and zero between the midpoint and nu(i+1).
		 *
		 * if Fnu(i) is zero we assume that the flux is zero up to the midpoint between
		 * nu(i) and nu(i+1) and Fnu(i+1) between the midpoint and nu(i+1). */

		if( ipHi == ipLo )
		{
			/* Do the case where the cloudy cell and its edges are between
			 * two adjacent stellar model points */
			if( StarPower[ipLo] == -FLT_MAX )
			{
				// case where bin at ipLo+1 has zero flux
				double StarHigh = min((StarEner[ipLo]+StarEner[ipLo+1])/2.,BinHigh);
				retval = StarFlux[ipLo]*max(StarHigh-BinLow,0.)/widflx;
			}
			else if( StarPower[ipLo] == FLT_MAX )
			{
				// case where bin at ipLo has zero flux (ipLo+1 possibly also)
				double StarLow = max((StarEner[ipLo]+StarEner[ipLo+1])/2.,BinLow);
				retval = StarFlux[ipLo+1]*max(BinHigh-StarLow,0.)/widflx;
			}
			else
			{
				/* do power-law interpolation  */
				retval = (realnum)(StarFlux[ipLo]*pow(anu/StarEner[ipLo],(double)StarPower[ipLo]));
			}
		}
		else
		{
			/* Do the case where the cloudy cell and its edges span two or more
			 * stellar model cells:  add segments with power-law interpolation up to
			 * do the averaging.*/

			double sum = 0.;

			/* ipHi points to stellar point at high end of cloudy continuum cell,
			 * if the Cloudy cell extends beyond the stellar grid, ipHi == nCont-1
			 * and the MIN2(ipHi,nCont-2) prevents access beyond allocated memory
			 * ipLo points to low end, above we asserted that 0 <= ipLo < nCont-1 */
			for( long i=ipLo; i <= MIN2(ipHi,nCont-2); i++ )
			{
				if( StarPower[i] == -FLT_MAX )
				{
					// case where bin at i+1 has zero flux
					double upper = min((StarEner[i]+StarEner[i+1])/2.,BinHigh);
					double wid = max(upper-max(BinLow,StarEner[i]),0.);
					sum += StarFlux[i]*wid;
					continue;
				}
				if( StarPower[i] == FLT_MAX )
				{
					// case where bin at i has zero flux (i+1 possibly also)
					double lower = max((StarEner[i]+StarEner[i+1])/2.,BinLow);
					double wid = max(min(BinHigh,StarEner[i+1])-lower,0.);
					sum += StarFlux[i+1]*wid;
					continue;
				}

				double v1, val, x1, x2;
				double pp1 = StarPower[i] + 1.;

				if( i == ipLo )
				{
					x1 = BinLow;
					x2 = StarEner[i+1];
					v1 = StarFlux[i]*pow(x1/StarEner[i],(double)StarPower[i]);
					/*v2 = StarFlux[i+1];*/
				}
				else if( i == ipHi )
				{
					x2 = BinHigh;
					x1 = StarEner[i];
					/*v2 = StarFlux[i]*pow(x2/StarEner[i],StarPower[i]);*/
					v1 = StarFlux[i];
				}
				else /*if( i > ipLo && i < ipHi )*/
				{
					x1 = StarEner[i];
					x2 = StarEner[i+1];
					v1 = StarFlux[i];
					/*v2 = StarFlux[i+1];*/
				}

				if( fabs(pp1) < 0.001 )
				{
					double z = log(x2/x1);
					double zp = z*pp1;
					val = x1*v1*z*(((zp/24.+1./6.)*zp+1./2.)*zp+1.);
				}
				else
				{
					val = pow(x2/x1,pp1) - 1.;
					val *= x1*v1/pp1;
				}
				sum += val;
			}

			retval = sum/widflx;
		}
	}
	return (realnum)retval;
}

inline long RebinFind(const realnum array[], /* array[nArr] */
		      long nArr,
		      realnum val)
{
	/* return ind(val) such that array[ind] <= val < array[ind+1],
	 *
	 * NB NB: this routine assumes that array[] increases monotonically !
	 *
	 * the first two clauses indicate out-of-bounds conditions and
	 * guarantee that when val1 <= val2, also ind(val1) <= ind(val2) */

	ASSERT( nArr > 1 );
	if( val < array[0] )
		return -1;
	else if( val > array[nArr-1] )
		return nArr-1;
	else
		return hunt_bisect(array, nArr, val);
}

template<class T>
void DumpAtmosphere(const char *fnam,
		    long imod,
		    long npar,
		    char names[MDIM][MNAM+1],
		    const vector<mpp>& telg,
		    long nflux,
		    const T anu[],
		    const realnum flux[])
{
	DEBUG_ENTRY( "DumpAtmosphere()" );

	// this will create the file if it doesn't exist yet...
	FILE *ioBUG = open_data( fnam, "a" );

	fprintf( ioBUG, "######## MODEL %ld", imod+1 );
	for( long nd=0; nd < npar; nd++ )
		fprintf( ioBUG, " %s %g", names[nd], telg[imod].par[nd] );
	fprintf( ioBUG, " ####################\n" );

	for( long k=0; k < nflux; ++k )
		fprintf( ioBUG, "%e %e\n", anu[k], flux[k] );

	fclose( ioBUG );
}
