/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseTest parse the test command */
#include "cddefines.h"
#include "iso.h"
#include "monitor_results.h"
#include "parser.h"

void ParseTest(Parser &p)
{
	char chStuff[INPUT_LINE_LENGTH];	

	DEBUG_ENTRY( "ParseTest()" );

	/* do a smoke test - 
	* generate test input stream for rapid testing of code */

	/* optional keyword PRINT will cause all the commands to be printed */
	int nPrintTest = p.nMatch("PRIN"  );

	bool lgH2 = p.nMatch(" H2 ");
	bool lgLARG = p.nMatch("LARG");
	bool lgMOLE = p.nMatch("MOLE");

	/* option to turn on the H2 molecule */
	if( lgH2 )
	{
		/* this will both enable the molecule and also set the threshold for actually
		* computing it to a very low value, so that it is actually used in ionized gas */
		sprintf( chStuff , "DATABASE H2 LIMIT -20 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseDatabaseH2(p);
	}

	/* use largest possible hydrogen atom */
	if( lgLARG )
	{
		sprintf( chStuff , "DATABASE H-LIKE ELEMENT HYDROGEN LEVELS LARGER  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseDatabaseISO(ipH_LIKE,p);
	}

	if( lgMOLE )
	{
#		if 0
		{
			/* hydrogen density */
			sprintf( chStuff , "TRACE TEMPERATURE CONVERGENCE " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseTrace(chStuff);
		}
#		endif

		/* hydrogen density */
		sprintf( chStuff , "HDEN 5 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseHDEN(p);

		/* make a constant temperature model */
		sprintf( chStuff , "CONSTANT TEMPER 50K  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseConstant(p);

		/* continuum to include full energy range */
		sprintf( chStuff , "TABLE ISM  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseTable( p);

		/* extinguish this continuum */
		sprintf( chStuff , "EXTINGUISH 23  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseExtinguish( p );
		
		/* stop in second zone, so we do use the zone increment logic */
		sprintf( chStuff , "STOP ZONE 2  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseStop(p);

		/* set thickness */
		sprintf( chStuff , "SET DR 0  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseSet(p);
		
		/* do Case B so lyman line pumping of H is not important */
		sprintf( chStuff , "CASE B  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseCaseB( p );

		/* set cosmic rays */
		//>>chng 12 apr 09, this uses old Williams CR backbround for continuity
		sprintf( chStuff , "COSMIC RAY BACKGROUND linear 0.1266 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseCosmicRays( p );

		/* add grains and ISM abundances */
		sprintf( chStuff , "ABUNDANCES ISM  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseAbundances( p );

		/* add grains and ism abundances */
		sprintf( chStuff , "CONSTANT GRAIN TEMPERATURE 20K  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseConstant( p );

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR EDEN 0.625 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* create series of monitor commands */
		if( lgH2 )
		{
			/*>>chng 13 apr 13, -1.528 to -1.450 we are not monitoring this in test suite */
			sprintf( chStuff , "MONITOR MOLECULAR FRACTION H2 -1.450 error 0.1 " );
		}
		else
		{
			/*>>chng 13 apr 13, -1.528 to -1.542 we are not monitoring this in test suite */
			sprintf( chStuff , "MONITOR MOLECULAR FRACTION H2 -1.542 error 0.1 " );
		}
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* create series of monitor commands */
		/*>>chng 13 apr 13, 0.987 to 0.954 we are not monitoring this in test suite */
		sprintf( chStuff , "MONITOR COLUMN CO 0.954 error 0.1 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);			

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR EDEN 0.625 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR HYDROGEN 1 TEMPERATURE 50K " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);
	}
	else
	{
		/* hydrogen density */
		sprintf( chStuff , "HDEN 4 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseHDEN( p );

		/* make a constant temperature model */
		sprintf( chStuff , "CONSTANT TEMPER 4  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseConstant( p );

		/* continuum to include full energy range */
		sprintf( chStuff , "TABLE AGN  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseTable( p);

		/* set ionization parameter */
		sprintf( chStuff , "IONIZATION PARAMETER -2  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseIonParI( p);

		/* use old abundances */
		sprintf( chStuff , "ABUNDANCES OLD SOLAR 84 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseAbundances(p);

		/* >>chng 02 apr 19, add this */
		/* set total Lyman continuum depth - this is to prevent caution
		 * that Lyman continuum was thin but expected to be thick */
		sprintf( chStuff , "STOP LYMAN OPTICAL -4  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseStop(p);
		
		/* stop in second zone, so we do use the zone increment logic */
		sprintf( chStuff , "STOP ZONE 2  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseStop(p);

		/* set thickness */
		sprintf( chStuff , "SET DR 0  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseSet(p);

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR HYDROGEN 1 IONIZATION -3.052 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		sprintf( chStuff , "MONITOR HELIUM 2 IONIZATION -1.076 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 13 feb 01, from -2.377 to -2.319. undo r6703 DR suppression */
		/*>>chng 13 apr 13, from -2.319 to -2.359, redo DR suppression */
		sprintf( chStuff , "MONITOR CARBON 2 IONIZATION -2.359 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.653 to -0.560. Badnell DR by default */
		/*>>chng 13 feb 01, from -0.596 to -0.565. undo r6703 DR suppression */
		/*>>chng 13 apr 13, from -0.565 to -0.586, redo DR suppression */
		sprintf( chStuff , "MONITOR CARBON 3 IONIZATION -0.586 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.348 to -0.373. Badnell DR by default */
		sprintf( chStuff , "MONITOR CARBON 4 IONIZATION -0.361 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.490 to -0.530. Badnell DR by default */
		sprintf( chStuff , "MONITOR CARBON 5 IONIZATION -0.514 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.800 to -0.861. Badnell DR by default */
		/*>>chng 11 jul 12, from -0.861 to -0.935, DR coll suppression */
		/*>>chng 13 feb 01, from -0.935 to -0.863. undo r6703 DR suppression */
		/*>>chng 13 apr 13, from -0.863 to -0.894, redo DR suppression */
		/*>>chng 15 may 14, from -0.894 to -0.864, CollisSuppres, final version from Dragan Nikolic */
		/*>>chng 15 oct 14, from -0.864 to -0.894, CollisSuppres, corrected version from Nikolic */
		/*>>chng 16 dec 20, from -0.894 to -0.865, add secondary autoionization to DR suppression */
		sprintf( chStuff , "MONITOR OXYGEN 3 IONIZATION  -0.865 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.180 to -0.157. Badnell DR by default */
		sprintf( chStuff , "MONITOR OXYGEN 4 IONIZATION -0.148 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.770 to -0.808. Badnell DR by default */
		/*>>chng 13 feb 01, from -0.790 to -0.807. undo r6703 DR suppression */
		/*>>chng 13 mar 13, from -0.807 to -0.798. undo DR suppression */
		sprintf( chStuff , "MONITOR OXYGEN 5 IONIZATION -0.798 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* >>chng 02 apr 19, from 0.7258 to 0.946, due to adding Lyman cont depth */
		/* >>chng 07 oct 22, from 0.946  to 1.108, resolve l-levels of h-like sequence */
		sprintf( chStuff , "MONITOR LINE \"CA B\" 4861.33 1.108 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		/* must have copy of this in chOrgCard, which is used by the routine to get lab */
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* >>chng 02 apr 19, from 2.4603 to 3.25, due to adding Lyman cont depth 
		 * >>chng 06 nov 17, asserted value from 3.25 to 3.11, drift over last year */
		/* >>chng 06 dec 01, from 3.11 to 2.72. Badnell DR by default */
		/* >>chng 07 oct 22, from 2.72 to 3.18, resolve l-levels of h-like sequence */
		/* >>chng 11 jul 12, from 3.18 to 2.695, DR coll suppression */
		/* >>chng 13 feb 01, from 2.695 to 3.189. undo r6703 DR suppression */
		/* >>chng 13 apr 13, from 3.189 to 2.961, redo DR suppression */
		/* >>chng 15 may 14, from 2.961 to 3.275, CollisSuppres, final version from Dragan Nikolic */
		/* >>chng 15 oct 14, from 3.275 to 3.063, CollisSuppres, corrected version from Nikolic */
		/* >>chng 16 dec 20, from 3.063 to 3.207, add secondary autoionization to DR suppression */
		sprintf( chStuff , "MONITOR LINE \"O  3\" 5006.84 3.207 " );

		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		/* must have copy of this in chOrgCard, which is used by the routine to get lab */
		p.setline(chStuff);
		ParseMonitorResults(p);

		sprintf( chStuff , "MONITOR HTOT -15.019" );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);
	}
	
	return;
}
