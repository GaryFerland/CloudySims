/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GRAINS_H_
#define GRAINS_H_

#include "rfield.h"
#include "grainvar.h"
#include "monointerp.h"

/** GrainDrive main routine to converge grains thermal solution */
void GrainDrive();

/** GrainDrift computes grains drift velocity */
void GrainDrift();

/** this routine is called by IterStart() */
void GrainStartIter();

/** this routine is called by IterRestart() */
void GrainRestartIter();

/** this routine is called by ParseSet() */
void SetNChrgStates(long);

/** startup routine for grains, called before first calculations, but after parsecommands */
void GrainsInit();

/** main routine for generating the grain diffuse emission */
void GrainMakeDiffuse();

/** main routine for quantum heating */
void qheat(/*@out@*/vector<double>&,/*@out@*/vector<double>&,/*@out@*/long*,size_t);

/** initialize interpolation arrays for grain enthalpy */
void InitEnthalpy();

struct GrainPar;

/** check validity of a refractive index file by checking the magic number */
bool lgValidRfiFile(const string& fnam);

/** check validity of a mixed medium file by checking the magic number */
bool lgValidMixFile(const string& fnam);

/** check validity of a size distribution file by checking the magic number */
bool lgValidSzdFile(const string& fnam);

/** check validity of an opacity file by checking the magic number */
bool lgValidOpcFile(const string& fnam);

/** mie_write_opc
 \param [in] *rfi_file 
 \param [in] *szd_file 
*/
void mie_write_opc(/*@in@*/const char*,/*@in@*/const char*,long int);
/**read in the *.opc file with opacities and other relevant information
\param *chFile
\param gp
*/
void mie_read_opc(/*@in@*/const char*,const GrainPar&);
/**set up Gaussian quadrature for arbitrary interval
\param nn
\param xbot
\param xtop
\param x[]
\param a[]
\param rr[]
\param ww[]
*/
void gauss_init(long int,double,double,const vector<double>&,const vector<double>&,vector<double>&,vector<double>&);
/** set up abscissas and weights for Gauss-Legendre intergration of arbitrary even order 
\param nn
\param x[]
\param a[]
*/
void gauss_legendre(long int,vector<double>&,vector<double>&);
/** find index ind such that min(xa[ind],xa[ind+1]) <= x <= max(xa[ind],xa[ind+1]).
 * xa is assumed to be strictly monotically increasing or decreasing.
 * if x is outside the range spanned by xa, lgOutOfBounds is raised and ind is set to -1
 * n is the number of elements in xa. 
 \param x
 \param xa[]
 \param n
 \param [out] *ind
 \param [out] *lgOutOfBounds
*/
void find_arr(double,const vector<double>&,long int,/*@out@*/long int*,/*@out@*/bool*);

/* grain_interpolate: interpolate on an array on the grain frequency mesh to create an array on the standard mesh */
template<typename T>
inline long grain_interpolate(const T arr1[], T arr2[], long n1) // arr1[n1], n1 <= gv.nflux
{
	DEBUG_ENTRY( "grain_interpolate()" );

	vector<double> d(n1-1), h(n1), g(n1);
	for( long k=0; k < n1-1; ++k ) 
	{
		h[k] = (gv.anu(k+1)-gv.anu(k));
		d[k] = (arr1[k+1]-arr1[k])/h[k];
	}
	MI_SetupData(d, h, g, n1);

	long i1=0, i2;
	double hh = h[0];
	for( i2=0; i2 < rfield.nflux; ++i2 )
	{
		double x = rfield.anu(i2);
		while( i1 < n1-1 && x >= gv.anu(i1+1) )
			hh = h[++i1];
		if( i1 == n1-1 )
			break;
		// apply Fritsch & Butland monotonic interpolation algorithm
		double t = (x - gv.anu(i1))/hh;
		arr2[i2] = arr1[i1]*h00(t) + hh*g[i1]*h10(t) +
			arr1[i1+1]*h00(1.0-t) - hh*g[i1+1]*h10(1.0-t);
	}
	// return the number of elements of arr2 that have been filled in
	return i2;
}
		
#endif /* GRAINS_H_ */
