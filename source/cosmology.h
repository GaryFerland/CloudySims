/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COSMOLOGY_H_
#define COSMOLOGY_H_

/** this is the current temperature of the cosmic background radiation<BR>
 >>refer	CMB	temp	Mather, J.C., Fixsen, D.J., Shafer, R.A., Mosier, C., & <BR>
 * >>refercon	Wilkinson, D.T. 1999, ApJ, 512, 511 */
#define	CMB_TEMP	2.725

/**GetDensity find the baryonic density at the given redshift 
\param z
*/
realnum GetDensity(realnum z);

/**GetHubbleFactor find the Hubble factor at the given redshift 
\param z
*/
realnum GetHubbleFactor(realnum z);

/** cosmology.h saves options and parameters relating to cosmology */
struct t_cosmology {

	realnum 
		redshift_current,
		redshift_start,
		redshift_step;

	realnum 
		omega_baryon,
		omega_rad,
		omega_lambda,
		omega_matter,
		omega_k;

	realnum
		h,
		H_0;

	realnum 
		f_He;

	bool lgDo;

	};
extern t_cosmology cosmology;


#endif /* COSMOLOGY_H_ */
