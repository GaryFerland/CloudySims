/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HYDRO_VS_RATES_H_
#define HYDRO_VS_RATES_H_

/** VS80 stands for Vriens and Smeets 1980<BR>
    This routine calculates thermally-averaged collision strengths.<BR>
    \param nHi 
    \param gHi
    \param IP_Ryd_Hi
    \param nLo
    \param gLo
    \param IP_Ryd_Lo
    \param Aul
    \param nelem
    \param Collider
    \param temp
*/
double CS_VS80( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, 
	double Aul, long nelem, long Collider, double temp );

/**hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients 
 \param ionization_energy_Ryd
 \param Te
 \param stat_level
 \param stat_ion
 */
double hydro_vs_coll_recomb( double ionization_energy_Ryd, double Te, double stat_level, double stat_ion );

/**hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients 
 \param ionization_energy_Ryd
 \param Te
 */
double hydro_vs_ioniz( double ionization_energy_Ryd, double Te );


/**Hion_coll_ioniz_ratecoef calculate hydrogenic ionization rates for all n, and Z
 \param ipISO the isoelectronic sequence 
 \param nelem element, >=1 since only used for ions<BR>
              nelem = 1 is helium the least possible charge
 \param n 	 principal quantum number, > 1<BR>
		 since only used for excited states<BR> 
 \param ionization_energy_Ryd
 \param temperature
 */
double Hion_coll_ioniz_ratecoef(
		long int ipISO ,
		long int nelem,
		long int n,
		double ionization_energy_Ryd,
		double temperature );

/**hydro_vs_deexcit generate hydrogenic collisional ionization rate coefficients 
 * for quantum number n 
 \param nHi
 \param gHi
 \param IP_Ryd_Hi 
 \param nLo
 \param gLo
 \param IP_Ryd_Lo 
 \param Aul
 */
double hydro_vs_deexcit( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul );

#endif /* HYDRO_VS_RATES_H_ */
