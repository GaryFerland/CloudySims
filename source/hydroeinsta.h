/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HYDROEINSTA_H_
#define HYDROEINSTA_H_

/**HydroEinstA calculates Einstein A's from  osillator strengths
\param lower
\param iupper
*/
double HydroEinstA(double lower, double iupper);

realnum hydro_transprob( long nelem, long ipHi, long ipLo );

#endif /* HYDROEINSTA_H_ */

/**hydro_transprob_collapsed_to_resolved calculates Einstein A's using Drake routine for collapsed to resolved levels
\param nelem
\param nHi
\param nLo
\param lLo
 */
realnum hydro_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo );

/**hydro_transprob_collapsed_to_collapsed calculates Einstein A's using Drake routine for collapsed levels
\param nelem
\param nHi
\param nLo
 */
realnum hydro_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo );
