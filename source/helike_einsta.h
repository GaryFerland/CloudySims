/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HELIKE_EINSTA_H_
#define HELIKE_EINSTA_H_

const int N_HE1_TRANS_PROB = 651;

void HelikeTransProbSetup();

/** compute energy diffference in wn and Aul for given line
 * return is 0 for success, 1 for failure 
\param nelem charge on the C scale, 1 is helium
\param Eff_nupper upper quantum numbers
\param lHi
\param sHi
\param jHi
\param Eff_nlower lower quantum numbers
\param lLo
\param sLo
\param jLo
*/
double he_1trans(  
		long nelem, 
		double Eff_nupper, long lHi, long sHi, long jHi,
		double Eff_nlower, long lLo, long sLo, long jLo);

/** helike_transprob_collapsed_to_collapsed calculates Einstein A coefficients for collapsed levels of he-like atoms
\param nelem
\param nHi
\param nLo
\param Enerwn
*/
double helike_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo, double Enerwn );

/**helike_transprob_collapsed_to_resolved calculates Einstein A coefficients for collapsed to resolved levels of he-like atoms
\param nelem
\param nHi
\param nLo
\param lLo
\param sLo
\param jLo
\param Enerwn
*/
double helike_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo, long sLo, long jLo, double Enerwn );



#endif /* HELIKE_EINSTA_H_ */
