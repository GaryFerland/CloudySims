/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

 /*  Created on: Sep 28, 2019
 *      Author: fran
 */

#ifndef SOURCE_ISO_CONDENSED_H_
#define SOURCE_ISO_CONDENSED_H_

#include "cddefines.h"
#include "container_classes.h"

void add_equilibrium_levels( long nelem, long ipISO, long mmax);
void get_mj(long mmax, vector<long>& mj);
void SourceTerms(vector<double> enerj,vector<double>& srcterm,multi_arr<double,2,C_TYPE>& einji);
double bfn(long j, double En);
double Kbn( long nelem, long ipISO, double Nhi, double gHi, long iplevel, long ntop, double IP_Ryd_Hi);
class integrand_matrix_element;



#endif /* SOURCE_ISO_CONDENSED_H_ */
