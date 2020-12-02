/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseDrive parse the drive command - drive calls to various subs */
/*DrvCaseBHS allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
/*DrvHyas allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
/*dgaunt drive gaunt factor routines by letting user query values */
#include "cddefines.h"
#include "trace.h"
#include "hydro_bauman.h"
#include "atmdat.h"  
#include "abund.h"
#include "rt_escprob.h"
#include "rt.h"
#include "mc_escape.h"
#include "parser.h"
#include "thirdparty.h"
#include "atmdat_gaunt.h"

/*dgaunt drive gaunt factor routines by letting user query values */
STATIC void dgaunt(void);

/*DrvHyas allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
STATIC void DrvHyas(void);

/* drive escape probability routines */
STATIC void DrvEscP( void );

void ParseDrive(Parser &p )
{
	bool lgEOL;
	long int n, 
	  i;
	double fac, 
	  zed;

	DEBUG_ENTRY( "ParseDrive()" );

	/* NB evolve all following names to style DrvSomething */

	/* option to drive cloudy, which one? */
	if( p.nMatch("FFMT") )
	{
		/* free format parser */
		char chInput[INPUT_LINE_LENGTH];
		fprintf( ioQQQ, " FFmtRead ParseDrive entered.  Enter number.\n" );
		lgEOL = false;
		while( !lgEOL )
		{
			if( read_whole_line( chInput , (int)sizeof(chInput) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " ParseDrive.dat error getting magic number\n" );
				cdEXIT(EXIT_FAILURE);
			}
			i = 1;
			fac = FFmtRead(chInput,&i,sizeof(chInput),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " FFmtRead hit the EOL with no value, return=%10.2e\n", 
				  fac );
				break;
			}
			else if( fac == 0. )
			{
				break;
			}
			else
			{
				fprintf( ioQQQ, " FFmtRead returned with value%11.4e\n", 
				  fac );
			}
			fprintf( ioQQQ, " Enter 0 to stop, or another value.\n" );
		}
		fprintf( ioQQQ, " FFmtRead ParseDrive exits.\n" );
	}

	else if( p.nMatch("CDLI") )
	{
		/* drive cdLine to check that it finds all the right lines, routine is in lines.c */
		trace.lgDrv_cdLine = true;
	}

	else if( p.nMatch(" E1 ") )
	{
		// option to drive exponential integral routines
		// first, special case given in Abramowitz & Stegun
		double tau = 1.275;
		for( i=0; i<50; ++i )
		{
			fprintf(ioQQQ,"tau\t%.3e\t exp-tau\t%.5e\t e1 tau\t%.5e  \t e2 "
				"\t%.5e \te2n %.5e \t e3\t%.5e \t e4\t%.5e \n",
				tau, sexp(tau), e1(tau), e2(tau), expn(2, tau),
				expn(3, tau), expn(4, tau) );
			tau = exp10(  ((double)i/4. - 9.) );
		}
		cdEXIT(EXIT_SUCCESS);
	}

	else if( p.nMatch("ESCA") )
	{
		if ( p.nMatch("VALI") )
		{
			double tau = p.FFmtRead();
			double a = p.FFmtRead();
			if (p.lgEOL())
			{
				a = 0.0;
			}
			double beta = p.FFmtRead();
			if (p.lgEOL())
			{
				beta = 0.0;
			}
			mc_escape(tau, a, beta);
			cdEXIT(EXIT_SUCCESS);
		}
		
		/* option to drive escape probability routines */
		DrvEscP( );
	}

	else if( p.nMatch("HYAS") )
	{
		/* option to drive Jason's hydrogen transition probabilities */
		DrvHyas();
	}

	else if( p.nMatch("GAUN") )
	{
		/* drive gaunt factor routine */
		if( p.nMatch("CHEC") )
		{
			double Z = p.FFmtRead();
			if( p.lgEOL() )
				Z = 1.;
			else
			{
				Z = floor(Z);
				if( Z <= 0. || Z > LIMELM )
				{
					fprintf( ioQQQ, " invalid value for charge %ld\n", long(Z) );
					cdEXIT(ES_FAILURE);
				}
			}
			dgaunt_check(Z);
		}
		else
			dgaunt();
	}

	else if( p.nMatch("PUMP") )
	{
		if( p.nMatch("VALI") )
		{
			double damp = p.FFmtRead();
			if (p.lgEOL())
			{
				damp = 0.0;
			}
			DrivePump(damp);
			cdEXIT(EXIT_SUCCESS);
		}
		char chInput[INPUT_LINE_LENGTH];
		lgEOL = false;
		fprintf( ioQQQ, " Continuum pump ParseDrive entered - Enter log tau\n" );
		while( !lgEOL )
		{
			if( read_whole_line( chInput , (int)sizeof(chInput) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " Parse Drive error getting optical depth\n" );
				cdEXIT(EXIT_FAILURE);
			}
			/* get tau */
			i = 1;
			double tau = FFmtRead(chInput,&i,sizeof(chInput),&lgEOL);
			if( lgEOL )
				break;
			tau = exp10(tau);
			fprintf( ioQQQ, " Tau =%11.4e\n", tau );
			double damp = 0.;
			fac = DrvContPump(tau, damp); 
			fprintf( ioQQQ, " ContPump =%11.4e\n", fac );
			fprintf( ioQQQ, " Enter null to stop, or another value.\n" );
		}
		fprintf( ioQQQ, " ContPump ParseDrive exits.\n" );
	}

	else if( p.nMatch("STAR") )
	{
		char chInput[INPUT_LINE_LENGTH];
		/* get starburst abundances */
		for( i=0; i < 40; i++ )
		{
			zed = ((double)i+1.)/4. + 0.01;
			sprintf( chInput, "starburst, zed=%10.4f", zed );
			p.setline(chInput);
			abund_starburst(p);
			fprintf( ioQQQ, "%8.1e", zed );
			for(n=0; n < LIMELM; n++)
				fprintf( ioQQQ, "%8.1e", abund.solar[n] );
			fprintf( ioQQQ, "\n" );
		}
	}

	else if( p.nMatch("VOIGT") )
	{
		string file;
		bool hasstr = ( p.GetQuote(file) == 0 );
		FILE *ioVOIGT = ioQQQ;
		if (hasstr)
			ioVOIGT = open_data(file.c_str(), "w");

		if( p.nMatch("TABLE") )
		{
			/* create tab-delimited table giving Voigt function */
			fprintf(ioVOIGT,"x \\ a");
			const realnum DampLogMin = -4., DampLogMax = 4.01;
			for( realnum damplog=DampLogMin; damplog<DampLogMax; ++damplog)
				fprintf(ioVOIGT,"\ta=%.3e",exp10(damplog));
			fprintf(ioVOIGT , "\n");
			
			for( realnum x=-2.; x<5.;x+=0.05)
			{
				realnum xlin = exp10(x);
				fprintf(ioVOIGT,"%.3e",xlin);
				for( realnum damplog=DampLogMin; damplog<DampLogMax; ++damplog)
				{
					realnum xval[1];
					xval[0] = xlin;
					realnum damp = exp10(damplog);
					realnum yval[1];
					VoigtH(damp,xval,yval,1);
					fprintf(ioVOIGT , "\t%.3e",yval[0]);
				}
				fprintf(ioVOIGT , "\n");
			}
		}
		else
		{
			/* Voigt function debugging print - parameter is damping constant a */
			realnum damp = (realnum)p.FFmtRead();
			if( p.lgEOL() )
			{
				fprintf( ioVOIGT, " The damping constant must appear on the print voigt command.  Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			
			const long NVOIGT=100;
			realnum xprofile[NVOIGT], profileVoigtH[NVOIGT];
			for( long i=0; i<NVOIGT; ++i )
				xprofile[i] = (realnum)i * 10.f / (realnum)NVOIGT;
			
			VoigtH( damp, xprofile, profileVoigtH, NVOIGT );
			
			fprintf(ioVOIGT,"\n    x        VoigtH\n");
			for( long int i=0; i<NVOIGT; ++i )
			{
				fprintf(ioVOIGT,"%.4e %.4e\n", xprofile[i], profileVoigtH[i] );
			}
		}

		if (ioVOIGT != ioQQQ)
			fclose(ioVOIGT);
		cdEXIT(EXIT_SUCCESS);
	} 

	else
	{
		fprintf( ioQQQ, 
			" Unrecognized key; keys are CDLIne, E1, ESCApe, FFMTread, GAUNt, "
			"HYAS, PUMP, STAR, and VOIGt.  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

/*DrvEscP user queries escape probability routines, which return values */
STATIC void DrvEscP( void )
{
	char chCard[INPUT_LINE_LENGTH];
	bool lgEOL;
	long i;
	double tau;

	DEBUG_ENTRY( "DrvEscP()" );

	/* this routine is enterd with the command escape probability, and
	 * drives the escape probability routine to compare answers */
	fprintf( ioQQQ, " Enter the log of the one-sided optical depth; line with no number to stop.\n" );

	lgEOL = false;
	while( !lgEOL )
	{
		if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
		{
			break;
		}

		i = 1;
		tau = FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
		if( lgEOL )
		{
			break;
		}

		tau = exp10(tau);
		fprintf( ioQQQ, "tau was %e\n", tau );
		fprintf( ioQQQ, " ESCINC=%11.3e\n", esc_PRD_1side(tau,1e-4) );
		fprintf( ioQQQ, " ESCCOM=%11.3e\n", esc_CRDwing_1side(tau,1e-4 ) );
		fprintf( ioQQQ, " ESCA0K2=%11.3e\n", esca0k2(tau) );

	}
	return;
}

/*DrvHyas allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
STATIC void DrvHyas(void)
{
	char chCard[INPUT_LINE_LENGTH];
	bool lgEOL;
	long int i, nHi, lHi, nLo, lLo;

	DEBUG_ENTRY( "DrvHyas()" );

	/* this routine is entered with the command DRIVE HYAS, and
	 * drives Jason's hydrogen einstein A routines */

	nHi = 1;
	/* nHi never lt 1 */
	while( nHi != 0 )
	{
		fprintf( ioQQQ, " Enter four quantum numbers (n, l, n', l'), null line to stop.\n" );
		if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
		{
			fprintf( ioQQQ, " error getting drvhyas \n" );
			cdEXIT(EXIT_FAILURE);
		}

		i = 1;
		nHi = (long int)FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
		if( lgEOL )
			break;

		lHi = (long int)FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " must be four numbers!\n" );
			break;
		}

		if( lHi >= nHi )
		{
			fprintf( ioQQQ, " l must be less than n!\n" );
			break;
		}

		nLo = (long int)FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " must be four numbers!\n" );
			break;
		}

		lLo = (long int)FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " must be four numbers!\n" );
			break;
		}

		if( lLo >= nLo )
		{
			fprintf( ioQQQ, " l must be less than n!\n" );
			break;
		}

		if( nLo > nHi )
		{
			long nTemp, lTemp;

			/* swap hi and lo */
			nTemp = nLo;
			lTemp = lLo;
			nLo = nHi;
			lLo = lHi;
			nHi = nTemp;
			lHi = lTemp;
		}

		fprintf( ioQQQ, " A(%3ld,%3ld->%3ld,%3ld)=%11.3e\n", 
			nHi, lHi, nLo, lLo,
			H_Einstein_A( nHi, lHi, nLo, lLo, 1 ) );

	}
	fprintf( ioQQQ, " Driver exits, enter next line.\n" );

	return;
}

/*dgaunt drive Gaunt factor routines by letting user query values */
STATIC void dgaunt()
{
	DEBUG_ENTRY( "dgaunt()" );

	/* this routine is entered with the command DRIVE GAUNT, and
	 * drives the Gaunt factor routine to check range
	 * */
	char chCard[INPUT_LINE_LENGTH];
	bool lgEOL;

	fprintf( ioQQQ, " Enter 0 to input temp, energy, and net charge, or 1 for u, gamma**2, and net charge.\n" );
	if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
	{
		fprintf( ioQQQ, " dgaunt error getting input line\n" );
		cdEXIT(EXIT_FAILURE);
	}
	long i = 1;
	int inputflag = (int)FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);

	if( inputflag == 0 )
	{
		fprintf( ioQQQ, " Enter the temperature (log if <=10), energy (Ryd), and net charge. Null line to stop.\n" );
		/* >>chng 96 july 07, got rid of statement labels replacing with do while
		 * */
		long ierror = 0;
		while( ierror == 0 )
		{
			if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " dgaunt error getting input line\n" );
				cdEXIT(EXIT_FAILURE);
			}
			i = 1;
			double alogte = FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
			/* the line may be trash but ierror will pick it up  */
			if( lgEOL  )
			{
				fprintf( ioQQQ, " Gaunt driver exits, enter next line.\n" );
				break;
			}
			/* numbers less than or equal to 10 are the log of the temperature */
			double TeNew;
			if( alogte > 10. )
				TeNew = alogte;
			else
				TeNew = exp10(alogte);

			double enerlin = FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
			if( lgEOL || enerlin == 0. )
				fprintf( ioQQQ, " Sorry, but there should be two more numbers, energy and charge.\n" );

			double z = FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
			if( lgEOL || z == 0. )
				fprintf( ioQQQ, " Sorry, but there should be a third number, charge.\n" );

			/* This is thermally averaged Gaunt factor */
			double mygaunt = t_gaunt::Inst().gauntff( long(z), TeNew, enerlin );
			fprintf( ioQQQ, " Using my routine, Gff= %.4e\n", mygaunt );
		}
	}
	else
	{
		/* this routine is entered with the command DRIVE GAUNT, and
		 * drives the Gaunt factor routine to check range
		 * */
		fprintf( ioQQQ, " Enter log u, log gamma2, and net charge. Null line to stop.\n" );
		long ierror = 0;
		while( ierror == 0 )
		{
			if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " dgaunt error getting input line\n" );
				cdEXIT(EXIT_FAILURE);
			}
			i = 1;
			double logu = FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
			/* the line may be trash but ierror will pick it up  */
			if( lgEOL )
			{
				fprintf( ioQQQ, " Gaunt driver exits, enter next line.\n" );
				break;
			}

			double loggamma2 = FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
			if( lgEOL )
				fprintf( ioQQQ, " Sorry, but there should be two more numbers, log gamma2 and charge.\n" );

			double z = FFmtRead(chCard,&i,sizeof(chCard),&lgEOL);
			if( lgEOL )
				fprintf( ioQQQ, " Sorry, but there should be another number, charge.\n" );

			/* This is my attempt to calculate thermally averaged Gaunt factors. */
			double mygaunt = t_gaunt::Inst().gauntff( long(z), TE1RYD*pow2(z)/exp10(loggamma2),
								  pow2(z)*exp10(logu-loggamma2) );
			fprintf( ioQQQ, " Using my routine, Gff= %.4e\n", mygaunt );
		}
	}
}
