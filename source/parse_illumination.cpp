/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseStop parse the stop command */
#include "cddefines.h"
#include "optimize.h"
#include "rfield.h"
#include "geometry.h"
#include "input.h"
#include "parser.h"

void ParseIllumination(Parser &p)
{
	DEBUG_ENTRY( "ParseIllumination()" );

	if( rfield.nShape < 1 )
	{
		fprintf(ioQQQ,"PROBLEM the illumination command has come before any "
			"radiation field shape commands.\nThis must come after the field"
			" is specified.\nSorry.\n");
		cdEXIT(EXIT_FAILURE);
	}

	rfield.Illumination[rfield.nShape-1] = Illumination::FORWARD;// default
	if( p.nMatch( "FORW" ) )
		rfield.Illumination[rfield.nShape-1] = Illumination::FORWARD;
	else if( p.nMatch( "REVE" ) )
	{
		rfield.Illumination[rfield.nShape-1] = Illumination::REVERSE;
		fprintf(ioQQQ, " Caution: the ILLUMINATION REVERSE command has not been implemented and has no effect.\n");
	}
	else if( p.nMatch( "SYMM" ) )
	{
		rfield.Illumination[rfield.nShape-1] = Illumination::SYMMETRIC;
		fprintf(ioQQQ, " Caution: the ILLUMINATION SYMMETRIC command has not been implemented and has no effect.\n");
	}

	// isotropic case
	if( p.nMatch( "ISOT" ) )
	{
		// isotropic illumination - following not meaningful
		rfield.OpticalDepthScaleFactor[rfield.nShape-1] = 1.;
		rfield.lgBeamed[rfield.nShape-1] = false;
	}
	else
	{
		double AngleIllumRadian = 0.;

		/* default is to specify an illumination angle */
		rfield.lgBeamed[rfield.nShape-1] = true;
		double a = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("illumination angle");

		/* default is degrees, but radian key accepted */
		if( p.nMatch( "RADI" ) )
			AngleIllumRadian = a;
		else
			AngleIllumRadian = a/RADIAN;
		if( AngleIllumRadian < 0. || AngleIllumRadian >= PI/2. )
		{
			fprintf( ioQQQ, " Angle of illumination must be between 0 and 90 degrees "
				"or 0 and pi/2 radians.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* this is effective path 1. / COS( theta ) -  so that
		 * dTau_eff = dTau_normal * geometry.DirectionalCosin */
		geometry.DirectionalCosin = (realnum)(1./cos(AngleIllumRadian));

		rfield.OpticalDepthScaleFactor[rfield.nShape-1] = 
			(realnum)(1./cos(AngleIllumRadian));

		/* vary option */
		if( optimize.lgVarOn )
		{
			/* no luminosity options on vary */
			optimize.nvarxt[optimize.nparm] = 1;
			if( p.nMatch( "RADI" ) )
			{
				strcpy( optimize.chVarFmt[optimize.nparm], "ILLUMINATION %f RADIAN" );
				optimize.vparm[0][optimize.nparm] = (realnum)AngleIllumRadian;
				optimize.varang[optimize.nparm][0] = 0.f;
				optimize.varang[optimize.nparm][1] = realnum(PI/2.);
			}
			else
			{
				strcpy( optimize.chVarFmt[optimize.nparm], "ILLUMINATION %f" );
				optimize.vparm[0][optimize.nparm] = (realnum)(AngleIllumRadian*RADIAN);
				optimize.varang[optimize.nparm][0] = 0.f;
				optimize.varang[optimize.nparm][1] = 90.f;
			}

			if( rfield.Illumination[rfield.nShape-1] == Illumination::FORWARD )
				strcat( optimize.chVarFmt[optimize.nparm], " FORWARD" );
			else if( rfield.Illumination[rfield.nShape-1] == Illumination::REVERSE )
				strcat( optimize.chVarFmt[optimize.nparm], " REVERSE" );
			else if( rfield.Illumination[rfield.nShape-1] == Illumination::SYMMETRIC )
				strcat( optimize.chVarFmt[optimize.nparm], " SYMMETRIC" );

			optimize.lgOptimizeAsLinear[optimize.nparm] = true;
			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			/* the increment in the first steps away from the original value */
			optimize.vincr[optimize.nparm] = 0.1f;
			optimize.nparm += 1;
		}
	}
}

