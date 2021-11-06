/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program that calls cloudy when used as a stand-alone program */
#include "cddefines.h"
#include "cddrive.h"

/*int main( int argc, char *argv[] )*/
int main( void )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	try {
		double hdenLimit , hdenInit , hden, TeInc , temp, TeLimit, TeInit , hdenInc;

		FILE *ioGRD ;
		char chLine[100];

		/* initialize the code for this run */
		cdInit();

		cdOutput( "CaseB.out" );

		// one of these must be selected
		string chName;
		int lgH1do = true;
		if( lgH1do )
			chName = "H1";
		int lgHe1do = false;
		if( lgHe1do )
			chName = "He1";
		int lgHe2do = false;
		if( lgHe2do )
			chName = "He2";

		string chBaseFilename = chName+"CaseB";

		/* calculation's grid points (density, temperature) will go to this file*/
		fprintf(ioQQQ, "Filename = %s\n",chBaseFilename.c_str());
		ioGRD = open_data(chBaseFilename+".grd","w");  
		fprintf(ioGRD,	"#Density\tTemperature\n");

		string chLineList = "LineList_"+chName+"_CaseB.dat";
		fprintf(stderr, "line list input is %s\n",chLineList.c_str());
		string chSaveData = chName+"CaseB.dat";
		fprintf(stderr, "saving predictions in %s\n",chSaveData.c_str());

		string chSaveLineListCommand = "save linelist absolute \""+chName+"CaseB.dat\" \"LineList_"+chName+"_CaseB.dat\" last no hash no clobber ";
		fprintf(stderr, "save line list command string is %s\n",chSaveLineListCommand.c_str());
		//exit(1);


		/* the range of density, and the increment in density, for this grid */
		hdenInit = 0.;
		hdenLimit = 6.;
		hdenInc = 2.;

		/* the range of log temperatures for the grid */
		TeInit = 4.;
		TeLimit = 5.;
		/* multiplicative inc - log */
		TeInc = 1.;

		/* set the density and temperature to the initial values */
		hden = hdenInit;
		temp = TeInit;


		while( temp < 1.01*TeLimit )
		{
			while( hden < 1.01 * hdenLimit )
			{
				/* initialize the code for this run */
				cdInit();

				/* option to not execute the code, uncomment when debugging setup */
				/*cdNoExec( );*/

				/* gas temperature for this calculation */
				sprintf(chLine,"constant temperature %f", temp);
				cdRead( chLine );

				/* gas density for this calculation */
				sprintf(chLine,"hden %f", hden);
				cdRead( chLine );

				/* only want the first zone */
				cdRead( "stop zone 1" );

				/* only H and He */
				cdRead( "init \"hheonly.ini\" " );

				/* set upline optical depth  */
				cdRead( "case B" );

				/* save linelist command */
				cdRead( chSaveLineListCommand.c_str() );	

				/* adjust thickness to get 4 pi J in line */
				sprintf(chLine,"set dr linear %e", 1./(pow(10.,2.*hden)*1.1 ) ) ;
				cdRead( chLine );


				/* an incident continuum must be specified to get the code
				 * to run at all - not very important since we will set
				 * the gas temperature */
				if( lgH1do )
				{
					cdRead( "laser 2 Ryd" );
					cdRead( "database levels H-like element Hydrogen resolved 10" );
					cdRead( "database levels H-like element Hydrogen collapsed 200" );
				}
				else if( lgHe1do )
				{
					cdRead( "laser 2 Ryd" );
					cdRead( "database levels He-like element Helium resolved 10" );
					cdRead( "database levels He-like element Helium collapsed 200" );
				}
				else if( lgHe2do )
				{
					cdRead( "laser 5 Ryd" );
					cdRead( "database levels H-like element Helium resolved 10" );
					cdRead( "database levels H-like element Helium collapsed 200" );
				}
				else
				{
					fprintf(ioQQQ,"did not find a Case B flag for H1, He1, or He2\n");
					exit(1);
				}

				cdRead( "ionization parameter 0" );
				/* actually call the code */
				if( cdDrive() )
					exit_status = ES_FAILURE;
				/* flush the output so we see it on the screen */
				fflush(ioQQQ);

				fprintf(ioGRD,"%.3e\t%.5f\n",	hden, temp );
				fprintf(stderr,"%.3e\t%.5f\n",	hden, temp );

				/************************* end lines with lf and flush it ***************/
				fflush(ioGRD );
				hden += hdenInc;
			}
			temp += TeInc;
			hden = hdenInit;
		}

		cdEXIT(exit_status);
	}
	catch( bad_alloc )
	{
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Most likely your computer "
			 "ran out of memory.\n Try monitoring the memory use of your run. Bailing out...\n" );
		exit_status = ES_BAD_ALLOC;
	}
	catch( out_of_range& e )
	{
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_OUT_OF_RANGE;
	}
	catch( domain_error& e )
	{
		fprintf( ioQQQ, " DISASTER - A vectorized math routine threw a domain_error. Bailing out...\n" );
		fprintf( ioQQQ, " What() = %s", e.what() );
		exit_status = ES_DOMAIN_ERROR;
	}
	catch( bad_assert& e )
	{
		MyAssert( e.file(), e.line() , e.comment() );
		exit_status = ES_BAD_ASSERT;
	}
#ifdef CATCH_SIGNAL
	catch( bad_signal& e )
	{
		if( ioQQQ != NULL )
		{
			if( e.sig() == SIGINT || e.sig() == SIGQUIT )
			{
				fprintf( ioQQQ, " User interrupt request. Bailing out...\n" );
				exit_status = ES_USER_INTERRUPT;
			}
			else if( e.sig() == SIGTERM )
			{
				fprintf( ioQQQ, " Termination request. Bailing out...\n" );
				exit_status = ES_TERMINATION_REQUEST;
			}
			else if( e.sig() == SIGILL )
			{
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
				exit_status = ES_ILLEGAL_INSTRUCTION;
			}
			else if( e.sig() == SIGFPE )
			{
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
				exit_status = ES_FP_EXCEPTION;
			}
			else if( e.sig() == SIGSEGV )
			{
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
				exit_status = ES_SEGFAULT;
			}
#			ifdef SIGBUS
			else if( e.sig() == SIGBUS )
			{
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
				exit_status = ES_BUS_ERROR;
			}
#			endif
			else
			{
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );
				exit_status = ES_UNKNOWN_SIGNAL;
			}

		}
	}
#endif
	catch( cloudy_exit& e )
	{
		if( ioQQQ != NULL )
		{
			ostringstream oss;
			oss << " [Stop in " << e.routine();
			oss << " at " << e.file() << ":" << e.line();
			if( e.exit_status() == 0 )
				oss << ", Cloudy exited OK]";
			else
				oss << ", something went wrong]";
			fprintf( ioQQQ, "%s\n", oss.str().c_str() );
		}
		exit_status = e.exit_status();
	}
	catch( std::exception& e )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}
	// generic catch-all in case we forget any specific exception above... so this MUST be the last one.
	catch( ... )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}

	cdPrepareExit(exit_status);

	return exit_status;
}

