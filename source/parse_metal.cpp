/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseMetal parse parameters on metal command */
#include "cddefines.h"
#include "input.h"
#include "optimize.h"
//#include "service.h"
#include "elementnames.h"
#include "grainvar.h"
#include "called.h"
#include "abund.h"
#include "parser.h"

STATIC int lgPrintMetalsDeplete;

STATIC void GetMetalsDeplete( Parser &p )
{
	/* save depletions read in from external file */
	STATIC realnum DepleteClassicSave[LIMELM];
	static int lgFirst=true;
	if( lgFirst )
	{
		lgFirst = false;
		string chFile;	/*file name for table read */
		if( p.nMatch( "\"" ) )
		{
			/*
			 * if a quote occurs on the line then get the ini file name
			 * this will also set the name in chCard and OrgCard to spaces
			 * so later keywords do not key off it
			 */
			if( p.GetQuote( chFile ) )
				p.StringError();
		}
		else
		{
			/* no quote appeared, so this is the default name, cloudy.ini */
			chFile = "ISM_CloudyClassic.dep";
		}
		string chPath = "abundances" + cpu.i().chDirSeparator() + chFile;
		STATIC FILE *ioDATA = open_data( chPath, "r" );	// will abort if not found

		if( lgPrintMetalsDeplete )
			fprintf(ioQQQ," First call, GetMetalsDeplete opened file %s \n", chPath.c_str() );

		// init with no depletion set, equal to 1
		for(int nelem=0; nelem<LIMELM; ++nelem)
			DepleteClassicSave[nelem] = 1.;

		string chLine;
		while( read_whole_line( chLine, ioDATA ) )
		{
			if( lgPrintMetalsDeplete )
				fprintf(ioQQQ, "line: %s", chLine.c_str() );

			/* field of stars or empty line end data */
			if( chLine.size() == 0 || chLine[0]=='\n' || chLine[0]=='\r' || chLine[0]=='*' )
				break;

			/* skip comment */
			if( chLine[0]=='#' )
				continue;

			size_t pp;
			/* erase EOL character */
			if( (pp = chLine.find_first_of("\n\r")) != string::npos )
				chLine.erase(pp);

			string chCAPS = chLine;
			caps(chCAPS);

			bool lgFound = false;
			for(int nelem=0; nelem<LIMELM; nelem++)
			{
				if( chCAPS.find(elementnames.chElementNameShort[nelem]) != string::npos )
				{
					lgFound = true;
					long int i = 1;
					bool lgEOL;

					DepleteClassicSave[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					if( lgPrintMetalsDeplete )
						fprintf(ioQQQ, " Parse:\t%s\t%.2f\n",
								elementnames.chElementNameShort[nelem] , DepleteClassicSave[nelem] );
					if( DepleteClassicSave[nelem]<0. || DepleteClassicSave[nelem]>1. )
					{
						fprintf(ioQQQ," The grain depletion for element %s was %.2e, it must be between 0 and 1\n",
								elementnames.chElementNameShort[nelem] , DepleteClassicSave[nelem] );
						cdEXIT(EXIT_FAILURE);
					}

					//we shouldn't need to continue once an element name is found on the line...
					break;
				}
			}
			if( !lgFound )
			{
				fprintf(ioQQQ, "PROBLEM in METALS DEPLETE: did not identify element name on this line: %s\n",
					chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* use derived depletions */
	for( long int i=0; i < LIMELM; i++ )
		abund.Depletion[i] = DepleteClassicSave[i];

}

STATIC double AX[LIMELM] , BX[LIMELM] , ZX[LIMELM];
STATIC int lgSetJenkins09[LIMELM];

STATIC void EvalJenkins(double DepJenkins09[LIMELM] , double Fstar )
{
	for( int nelem=0; nelem<LIMELM; ++nelem )
	{
		if( lgSetJenkins09[nelem] )
			DepJenkins09[nelem] = pow(10., BX[nelem] + AX[nelem]*(Fstar-ZX[nelem]) );
		else
			DepJenkins09[nelem] = 1.;
	}
}

STATIC void GetJenkins09(int lgPrtJenkins09 , double DepJenkins09[LIMELM] , double Fstar , Parser &p )
{
	static int lgFirst=true;
	if( lgFirst )
	{
		lgFirst = false;
		string chFile;	/*file name for table read */
		if( p.nMatch( "\"" ) )
		{
			/*
			 * if a quote occurs on the line then get the ini file name
			 * this will also set the name in chCard and OrgCard to spaces
			 * so later keywords do not key off it
			 */
			if( p.GetQuote( chFile ) )
				p.StringError();
		}
		else
		{
			/* no quote appeared, so this is the default name, cloudy.ini */
			chFile = "ISM_Jenkins09_Tab4.dep";
		}
		string chPath = "abundances" + cpu.i().chDirSeparator() + chFile;
		STATIC FILE *ioDATA = open_data( chPath, "r" );	// will abort if not found
		if( lgPrtJenkins09 )
			fprintf(ioQQQ," First call, opened file %s \n", chPath.c_str() );
		// init with no depletion set
		for(int nelem=0; nelem<LIMELM; ++nelem)
			lgSetJenkins09[nelem] = false;

		string chLine;
		while( read_whole_line( chLine, ioDATA ) )
		{
			if( 0 && lgPrtJenkins09 )
				fprintf(ioQQQ, "line: %s", chLine.c_str() );

			/* field of stars or empty line end data */
			if( chLine.size() == 0 || chLine[0]=='\n' || chLine[0]=='\r' || chLine[0]=='*' )
				break;

			/* skip comment */
			if( chLine[0]=='#' )
				continue;

			size_t pp;
			/* erase EOL character */
			if( (pp = chLine.find_first_of("\n\r")) != string::npos )
				chLine.erase(pp);

			string chCAPS = chLine;
			caps(chCAPS);

			bool lgFound = false;
			for(int nelem=0; nelem<LIMELM; nelem++)
			{
				if( chCAPS.find(elementnames.chElementNameShort[nelem]) != string::npos )
				{
					lgFound = true;
					long int i = 1;
					bool lgEOL;

					//STATIC double AX[LIMELM] , BX[LIMELM] , ZX[LIMELM];
					lgSetJenkins09[nelem] = true;
					AX[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					BX[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					ZX[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					if( lgPrtJenkins09 )
						fprintf(ioQQQ, " Parse:\t%s\t%.2f\t%.2f\t%.2f\n",
								elementnames.chElementNameShort[nelem] , AX[nelem] , BX[nelem] , ZX[nelem] );

					//we shouldn't need to continue once an element name is found on the line...
					break;
				}
			}
			if( !lgFound )
			{
				fprintf(ioQQQ, "PROBLEM in ABUNDANCES: did not identify element name on this line: %s\n",
					chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* option to print table of depletions */
	if( lgPrtJenkins09 )
	{
		fprintf(ioQQQ," Parse: end of data in file, report of range of depletions follows:\n\n Fstar");
		for( int nelem=0; nelem<LIMELM; ++nelem )
			fprintf(ioQQQ,"\t%s" , elementnames.chElementNameShort[nelem] );
		fprintf(ioQQQ,"\n");

		/* print range of Fstar, do not change original value */
		for( double FstarLoc=0; FstarLoc<1.01; FstarLoc+=0.1 )
		{
			fprintf(ioQQQ,"%.2f",FstarLoc);
			EvalJenkins( DepJenkins09 , FstarLoc );
			for( int nelem=0; nelem<LIMELM; ++nelem )
				fprintf(ioQQQ,"\t%.2e" , DepJenkins09[nelem] );
			fprintf(ioQQQ,"\n");
		}
	}

	/* return correct depletions for specified Fstar */
	EvalJenkins( DepJenkins09 , Fstar );
	if( lgPrtJenkins09 )
		fprintf(ioQQQ," Jenkins09 later call \n");

}
void ParseMetal(Parser &p)
{
	DEBUG_ENTRY( "ParseMetal()" );

	/* parse the metals command */
	abund.lgAbnSolar = false;	
	if( p.nMatch("DEPL") )
	{
		/* keyword depletion is present
		 * deplete by set of scale factors */
		abund.lgDepln = true;

		/* option to use Jenkins 2009ApJ...700.1299J fits to ISM depletion
		 * syntax is METALS DEPLETE JENKINS 2009 FSTAR=0.5 */
		if( p.nMatch("JENKINS") && p.nMatch("2009") )
		{
			/* read Fstar, Jenkins' basic parameter, and confirm it is in range */
			realnum Fstar = 0;
			Fstar = (realnum)p.FFmtRead();
			if( fabs( Fstar-2009) > 0.1 )
			{
				fprintf(ioQQQ," The first number of METALS DEPLETE JENKINS 2009 "
						"must be 2009, it was %.2f\n Sorry\n\n ", Fstar);
				cdEXIT( EXIT_FAILURE );
			}
			Fstar = (realnum)p.FFmtRead();

			/* sort out whether log */
			bool lgLogOn = false;
			if( p.nMatch(" LOG") )
				lgLogOn = true;
			else if( p.nMatch("LINE") )
				lgLogOn = false;

			/* save log of Fstar since needed for vary option */
			double FsLog;
			if( Fstar <0. || lgLogOn )
			{
				FsLog = Fstar;
				Fstar = exp10(Fstar);
			}
			else
			{
				FsLog = log10(Fstar);
			}

			/* check range of final Fstar */
			if(Fstar<0 || Fstar>1. )
			{
				fprintf(ioQQQ, "Fstar on METALS DEPLETE JENKINS command must be "
						"between 0 and 1 and was %.2f\n Sorry\n\n", Fstar );
				cdEXIT( EXIT_FAILURE );
			}


			/* print option, to report details and print table of values */
			STATIC int lgPrtJenkins09=false;
			if( p.nMatch("PRINT")  )
				lgPrtJenkins09 = true;

			if( lgPrtJenkins09 )
				fprintf(ioQQQ,"\n Jenkins 2009, print set, found Fstar = %.2e\n" , Fstar);

			STATIC double DepJenkins09[LIMELM];
			GetJenkins09(lgPrtJenkins09 , DepJenkins09 , Fstar , p);

			/* use derived depletions */
			for( long int i=0; i < LIMELM; i++ )
				abund.Depletion[i] = DepJenkins09[i];

			/* vary option */
			if( optimize.lgVarOn )
			{
				strcpy( optimize.chVarFmt[optimize.nparm], "METALS DEPLETE JENKINS 2009=%f LOG " );
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				optimize.vparm[0][optimize.nparm] = (realnum)FsLog;
				optimize.nvarxt[optimize.nparm] = 1;
				/* initial increment is 0.1 and allowed range i 0 to 1 */
				optimize.vincr[optimize.nparm] = 0.1;
				optimize.varang[optimize.nparm][0] = 0.f;
				optimize.varang[optimize.nparm][1] = 1.f;

				++optimize.nparm;
			}
		}
		else
		{
			/* metals deplete command - not Jenkins */
			if( p.nMatch("PRINT")  )
				lgPrintMetalsDeplete = true;
			else
				lgPrintMetalsDeplete = false;

			/* no keyword - use stored abundances */
			GetMetalsDeplete( p );
			for( long int i=0; i < LIMELM; i++ )
			{
				if( lgPrintMetalsDeplete )
					fprintf(ioQQQ,"DEBUGGG depnew %s\t%.3e\n", elementnames.chElementName[i], abund.Depletion[i] );
			}
		}

	}
	else
	{
		/* metal depletion factor, if negative then it is the log */
		abund.ScaleMetals = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("the depletion factor");

		/* sort out whether log */
		bool lgLogOn = false;
		if( p.nMatch(" LOG") )
		{
			lgLogOn = true;
		}
		else if( p.nMatch("LINE") )
		{
			lgLogOn = false;
		}

		double dmlog;
		if( abund.ScaleMetals <= 0. || lgLogOn )
		{
			dmlog = abund.ScaleMetals;
			abund.ScaleMetals = exp10(abund.ScaleMetals);
		}
		else
		{
			dmlog = log10(abund.ScaleMetals);
		}

		/* option to vary grain abundance as well */
		bool lgGrains;
		if( p.nMatch("GRAI") )
		{
			lgGrains = true;
			gv.GrainMetal = abund.ScaleMetals;
		}
		else
		{
			lgGrains = false;
			gv.GrainMetal = 1.;
		}

		/* vary option */
		if( optimize.lgVarOn )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "METALS= %f LOG" );
			if( lgGrains )
				strcat( optimize.chVarFmt[optimize.nparm], " GRAINS" );
			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			optimize.vparm[0][optimize.nparm] = (realnum)dmlog;
			optimize.vincr[optimize.nparm] = 0.5;
			optimize.nvarxt[optimize.nparm] = 1;
			++optimize.nparm;
		}
	}
	return;
}
