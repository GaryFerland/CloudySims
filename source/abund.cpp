/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "abund.h"
t_abund abund;

void t_abund::zero()
{
	DEBUG_ENTRY( "t_abund::zero()" );
	/* some vars dealing with depletions due to grains */
	for( long nelem=0; nelem < LIMELM; nelem++ )
	{
		/* depletion scale factors */
		depset[nelem] = 1.;
	}
	lgDepln = false;
	ScaleMetals = 1.;
	
}
