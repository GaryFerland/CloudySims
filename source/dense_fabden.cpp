/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*dense_fabden called by dlaw command, returns density for any density law */
#include "cddefines.h"
#include "rfield.h"
#include "dense.h"

/*dense_fabden implements the dlaw command, returns density using
 * current position and up to ten parameters on dlaw command line */
double dense_fabden(double radius, 
  double depth)
{

	fprintf(ioQQQ,"PROBLEM The DLAW command requires a user-defined version of dense_fabden.cpp\n");
	fprintf(ioQQQ,"dense_fabden was called with radius=%e.2 and depth=%.2e\n", radius, depth);
	fprintf(ioQQQ,"Edit this file to create your own function and remove this warning.\n");
	cdEXIT(EXIT_FAILURE);
	//return fabden_v;
}
