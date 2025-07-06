/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef STARS_H_
#define STARS_H_


/** This is the largest number of dimensions that the grid can have */
static const int MDIM = 4;

/** This is the maximum length a dimension label may have [e.g., Teff, log(g)] */
static const int MNAM = 6;

/** interpolation mode; usually IM_RECT_GRID, but for CoStar models 4 other options
 * exist; each name gives the first and second parameter used in the interpolation */
typedef enum {
	IM_ILLEGAL_MODE=-1, IM_RECT_GRID, IM_COSTAR_TEFF_MODID,
	IM_COSTAR_TEFF_LOGG, IM_COSTAR_MZAMS_AGE, IM_COSTAR_AGE_MZAMS
} IntMode;

typedef enum {
	TL_OBSTAR, TL_BSTAR, TL_OSTAR
} tl_grid;

typedef enum {
	SB_TOTAL, SB_STELLAR, SB_NEBULAR
} sb_mode;

struct process_counter
{
	int nFound;
	int notProcessed;
	int nOK;
	int nFail;
	process_counter() : nFound(0), notProcessed(0), nOK(0), nFail(0) {}
};

/** List all the available TABLE STAR <grid> commands by checking installed *.mod files */
void AtmospheresAvail();

/** AtlasInterpolate interpolate on atlas model atmospheres, by K Volk */
long AtlasInterpolate(double val[], /* val[nval] */
		      long *nval,
		      long *ndim,
		      const string& chMetalicity,
		      const string& chODFNew,
		      bool lgList,
		      double *Tlow,
		      double *Thigh);

/** CoStarCompile rebin costar stellar models to match energy grid of code*/
bool CoStarCompile();
/** CoStarInterpolate read in and interpolate on Werner grid of PN atmospheres, by K Volk */
long CoStarInterpolate(double val[], /* requested model parameters */
		       long *nval,
		       long *ndim,
		       IntMode imode, /* which interpolation mode is requested */
		       bool lgHalo,  /* flag indicating whether solar (==0) or halo (==1) abundances */
		       bool lgList,
		       double *val0_lo,
		       double *val0_hi);

/** GridCompile rebin user supplied stellar models to match energy grid of code */
bool GridCompile(const string& InName);
bool GridCompile(const string& fnam, process_counter& pc);
/** GridInterpolate read in and interpolate on user supplied grid of atmospheres */
long GridInterpolate(double val[], /* val[nval] */
		     long *nval,
		     long *ndim,
		     const string& InName,
		     bool lgList,
		     double *Tlow,
		     double *Thigh);

/** HaardtMadauInterpolate read in and interpolate on Haardt & Madau SEDs */
long HaardtMadauInterpolate(double val,
			    int version,
			    bool lgQuasar,
			    double *zlow,
			    double *zhigh);

/** KhaireSrianandInterpolate read in and interpolate on Khaire & Srianand SEDs */
long KhaireSrianandInterpolate(double val,
			       int Q,
			       double *zlow,
			       double *zhigh);

/** Kurucz79Interpolate read in and interpolate on Kurucz 1979 grid of atmospheres */
long Kurucz79Interpolate(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgList,
			 double *Tlow,
			 double *Thigh);

/** MihalasInterpolate read in and interpolate on Mihalas grid of atmospheres */
long MihalasInterpolate(double val[], /* val[nval] */
			long *nval,
			long *ndim,
			bool lgList,
			double *Tlow,
			double *Thigh);

/** RauchCompile create ascii and mod files for Rauch atmospheres
 * return 0 if success, 1 if failure */
bool RauchCompile();
/** RauchInterpolateHydr get one of the Rauch pure hydrogen model atmospheres */
long RauchInterpolateHydr(double val[], /* val[nval] */
			  long *nval,
			  long *ndim,
			  bool lgList,
			  double *Tlow,
			  double *Thigh);
/** RauchInterpolateHelium get one of the Rauch pure helium model atmospheres */
long RauchInterpolateHelium(double val[], /* val[nval] */
			    long *nval,
			    long *ndim,
			    bool lgList,
			    double *Tlow,
			    double *Thigh);
/** RauchInterpolateHpHe get one of the Rauch hydrogen plus helium model atmospheres */
long RauchInterpolateHpHe(double val[], /* val[nval] */
			  long *nval,
			  long *ndim,
			  bool lgList,
			  double *Tlow,
			  double *Thigh);
/** RauchInterpolatePG1159 get one of the Rauch PG1159 model atmospheres */
long RauchInterpolatePG1159(double val[], /* val[nval] */
			    long *nval,
			    long *ndim,
			    bool lgList,
			    double *Tlow,
			    double *Thigh);
/** RauchInterpolateCOWD get one of the Rauch C/O white dwarf model atmospheres */
long RauchInterpolateCOWD(double val[], /* val[nval] */
			  long *nval,
			  long *ndim,
			  bool lgList,
			  double *Tlow,
			  double *Thigh);
/** RauchInterpolateHCa get one of the Rauch H-Ca model atmospheres, originally by K. Volk */
long RauchInterpolateHCa(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgHalo,
			 bool lgList,
			 double *Tlow,
			 double *Thigh);
/** RauchInterpolateHNi get one of the Rauch H-Ni model atmospheres */
long RauchInterpolateHNi(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgHalo,
			 bool lgList,
			 double *Tlow,
			 double *Thigh);

/** Create .ascii file out of Starburst99 output */
bool StarburstInitialize(const string& chInName,
			 const string& chOutName,
			 sb_mode mode);
/** StarburstCompile, rebin Starburst99 model output to match energy grid of code */
bool StarburstCompile();

/** TlustyInterpolate get one of the Tlusty OSTAR2002 model atmospheres */
long TlustyInterpolate(double val[], /* val[nval] */
		       long *nval,
		       long *ndim,
		       tl_grid tlg,
		       const string& chMetalicity,
		       bool lgList,
		       double *Tlow,
		       double *Thigh);

/** WernerInterpolate read in and interpolate on Werner grid of PN atmospheres, by K Volk */
long WernerInterpolate(double val[], /* val[nval] */
		       long *nval,
		       long *ndim,
		       bool lgList,
		       double *Tlow,
		       double *Thigh);

/** WMBASICInterpolate read in and interpolate on WMBASIC grid of hot star atmospheres */
long WMBASICInterpolate(double val[], /* val[nval] */
			long *nval,
			long *ndim,
			bool lgList,
			double *Tlow,
			double *Thigh);

#endif /* STARS_H_ */
