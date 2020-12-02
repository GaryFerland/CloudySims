/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HeLikeError fills uncertainty arrays */
#include "cddefines.h" 
#include "iso.h"

/* This routine handles errors when that option is turned on (via the command
 * "atom he-like error generation" */
void iso_put_error(long int ipISO,
			  long int nelem,
			  long int ipHi,
			  long int ipLo,
			  long int whichData,
			  realnum errorOpt,
			  realnum errorPess)
{

	DEBUG_ENTRY( "iso_put_error()" );

	if( iso_ctrl.lgRandErrGen[ipISO] )
	{
		/* whichData is either IPRAD, IPCOLLIS, or IPENERGY */
		ASSERT( whichData <= 2 );
		ASSERT( ipISO < NISO );
		ASSERT( nelem < LIMELM );
		ASSERT( ipHi <= iso_sp[ipISO][nelem].numLevels_max );
		ASSERT( ipLo <= iso_sp[ipISO][nelem].numLevels_max );
		ASSERT( errorOpt >= 0. );
		ASSERT( errorPess >= 0. );
		
		if( !iso_ctrl.lgPessimisticErrors )
			iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[whichData] = errorOpt;
		else
			iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[whichData] = errorPess;
	}
	return;
}

void iso_error_generation( long ipISO, long nelem )
{
	long ipHi, ipLo, typeOfRate;

	DEBUG_ENTRY( "iso_error_generation()" );

	//iso_sp[ipISO][nelem].ex[iso_sp[ipISO][nelem].numLevels_max][iso_sp[ipISO][nelem].numLevels_max].ErrorFactor[IPRAD] =
		//(realnum)MyGaussRand( iso_sp[ipISO][nelem].ex[iso_sp[ipISO][nelem].numLevels_max][iso_sp[ipISO][nelem].numLevels_max].Error[IPRAD] );

	for( ipHi=1; ipHi<= iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
	{
		/* >>chng 06 mar 15, the upper limit incorrectly went to numLevels_max */
		for( ipLo=0; ipLo< ipHi; ipLo++ )
		{
			for( typeOfRate=0; typeOfRate<=1; typeOfRate++ )
			{
				if( iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[typeOfRate] >= 0. )
				{
					iso_sp[ipISO][nelem].ex[ipHi][ipLo].ErrorFactor[typeOfRate] =
						(realnum)MyGaussRand( iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[typeOfRate] );
					ASSERT( iso_sp[ipISO][nelem].ex[ipHi][ipLo].ErrorFactor[typeOfRate] > 0. );
				}
				else
				{
					iso_sp[ipISO][nelem].ex[ipHi][ipLo].ErrorFactor[typeOfRate] = 1.0f;
				}
			}
		}
	}

	/* set flag saying that error generation has been done.  */
	iso_sp[ipISO][nelem].lgErrGenDone = true;
	return;
}
