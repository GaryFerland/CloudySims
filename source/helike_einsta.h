/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HELIKE_EINSTA_H_
#define HELIKE_EINSTA_H_

#define N_HE1_TRANS_PROB	651

#define MAX_TP_INDEX	110

void HelikeTransProbSetup( void );

/** compute energy diffference in wn and Aul for given line
 * return is 0 for success, 1 for failure 
\param nelem charge on the C scale, 1 is helium
\param Enerwn energy difference in wavenumber
\param Eff_nupper upper quantum numbers
\param Eff_nlower lower quantum numbers
\param lHi
\param sHi
\param jHi
\param lLo
\param sLo
\param jLo
\param *error1
\param *error2
*/
double he_1trans(  
		long nelem, 
		double Enerwn, 
		double Eff_nupper, long lHi, long sHi, long jHi,
		double Eff_nlower, long lLo, long sLo, long jLo,
		realnum *error1, realnum *error2 );

/** Every bit of this routine is based upon the singlet-triplet mixing formalism given in
* >>refer He   FSM     Drake, G. W. F. 1996, in Atomic, Molecular, \& Optical Physics Handbook,
* >>refercon   ed. G. W. F. Drake (New York: AIP Press).       
* That formalism mixes the levels themselves, but since this code is not J-resolved, we simulate
* that by mixing only the transition probabilities.  We find results comparable to those calculated
* in the fully J-resolved model spearheaded by Rob Bauman, and described in
* >>refer      He      FSM     Bauman, R. P., Porter, R. L., Ferland, G. J., \& MacAdam, K. B. 2005, ApJ, accepted 
\param nelem
\param ipLoSing
\param ipHiSing
*/
void DoFSMixing( long nelem, long ipLoSing, long ipHiSing );

#endif /* HELIKE_EINSTA_H_ */
