/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*plot master routine to generate some sort of plot */
#include "cddefines.h"
#include "iterations.h"
#include "called.h"
#include "input.h"
#include "rfield.h"
#include "trace.h"
#include "radius.h"
#include "geometry.h"
#include "opacity.h"
#include "dense.h"
#include "hcmap.h"
#include "plot.h"

static const int IHI = 59;
static const int IWID = 121;

t_plotCom plotCom;

/*pltcon generate plot of continuum array */
STATIC void pltcon(long int np, 
  const char *chCall);

/*pltmap generate plot of heating and cooling map */
STATIC void pltmap(long int np, 
  const char *chCall);

/*pltopc generate plot of local gas opacity */
STATIC void pltopc(long int np, 
  const char *chCall);

/* this is the base routine that actually makes the plots, called by above */
STATIC void pltr(const double[],const double[],long,double,double,double,double,
	char,const char*,long,bool);


void plot(const char *chCall)
{
	long int np;

	DEBUG_ENTRY( "plot()" );

	/* return if this is not the last iteration, or a plot not required,
	 * or we are not speaking */
	if( !plotCom.lgPlotON || !called.lgTalk )
	{ 
		return;
	}

	if( !iterations.lgLastIt && (strcmp(chCall,"FIRST") != 0) )
	{ 
		return;
	}

	/* loop over all the requested plots */
	for( np=0; np < plotCom.nplot; np++ )
	{
		/* series of tests to determine which type of plot we will do */
		if( strcmp(plotCom.chPType[np]," MAP") == 0 )
		{
			/* thermal map */
			pltmap(np,chCall);
		}
		else if( strcmp(plotCom.chPType[np] ,"CONT") == 0 || 
			      strcmp(plotCom.chPType[np] ,"CRAW") == 0 || 
					strcmp(plotCom.chPType[np] ,"DIFF") == 0 || 
					strcmp(plotCom.chPType[np] ,"REFL") == 0 || 
					strcmp(plotCom.chPType[np] ,"EMIT") == 0 || 
					strcmp(plotCom.chPType[np] ,"CPHT") == 0 || 
					strcmp(plotCom.chPType[np] ,"OUTW") == 0 )
		{
			/* this is a contiuum plot of some kind */
			pltcon(np,chCall);
		}

		else if( 
			strcmp(plotCom.chPType[np] ,"OPAA") == 0 || 
			strcmp(plotCom.chPType[np] ,"OPAS") == 0 || 
			strcmp(plotCom.chPType[np] ,"OPAT") == 0 )
		{
			/*  absorption, scattering, or total opacity */
			pltopc(np,chCall);
		}
		else
		{
			fprintf( ioQQQ, " PLOT type=%4.4s not known.  STOP\n", 
			  plotCom.chPType[np] );
			cdEXIT(EXIT_FAILURE);
		}
	}

	return;
}

/*pltcon generate plot of continuum array */

STATIC void pltcon(
  long int np, 
  const char *chCall)
{
	char chSymPlt2[3], 
	  chXtitle[23];
	char chSym, 
	  chSymPlt1;
	long int i;
	double contin, 
	  ymin2;
	static double xmax, 
	  xmin, 
	  ymax, 
	  ymin;
	static double *y/*[rfield.nflux_with_check]*/, 
	  *y2/*[rfield.nflux_with_check]*/; 

	DEBUG_ENTRY( "pltcon()" );

	if( strcmp(chCall,"FIRST") == 0 )
	{ 
		return;
	}

	xmin = rfield.anulog10(0);
	xmin = MAX2((double)plotCom.pltxmn[np],xmin);
	xmax = rfield.anulog10(rfield.nflux-1);
	xmax = MIN2(xmax,(double)plotCom.pltxmx[np]);

	if( plotCom.lgPltTrace[np] )
	{
		fprintf( ioQQQ, " XMIN, XMAX=%12.4e%12.4e NFLUX=%4ld\n", 
		  xmin, xmax, rfield.nflux );
	}

	if( strcmp(plotCom.chPType[np],"REFL") == 0 && geometry.lgSphere )
	{
		fprintf( ioQQQ, " Reflected continuum not computed when SPHERE set.\n" );
		return;
	}

	y = (double*)MALLOC((size_t)rfield.nflux_with_check*sizeof(double) );
	y2 = (double*)MALLOC((size_t)rfield.nflux_with_check*sizeof(double) );

	/* these will be the default symbols for first and second plot */
	chSymPlt1 = '.';
	strcpy( chSymPlt2, "o " );
	ymin = FLT_MAX;
	ymin2 = FLT_MAX;
	ymax = -FLT_MAX;
	for( i=0; i < rfield.nflux; i++ )
	{
		if( rfield.anulog10(i) > xmin && rfield.anulog10(i) < xmax )
		{
			if( strcmp(plotCom.chPType[np],"CONT") == 0 )
			{
				y[i] = log10(MAX2(rfield.flux_total_incident[0][i]/rfield.widflx(i)*
				  rfield.anu2(i),1e-37));
				/* >>chng 01 jul 13, add rfield.ConEmitReflec[0][i] */
				contin = rfield.flux[0][i] + rfield.ConEmitOut[0][i]*geometry.covgeo + rfield.ConEmitReflec[0][i];
				y2[i] = MAX2((contin/rfield.widflx(i)+(rfield.outlin[0][i]+rfield.outlin_noplot[i])/
					rfield.anu(i)*geometry.covgeo)*
					rfield.anu2(i)*radius.r1r0sq,1e-37);
				y2[i] = log10(y2[i]);
			}
			else if( strcmp(plotCom.chPType[np],"CPHT") == 0 )
			{
				/*    plot continuum as photons */
				y[i] = log10(MAX2(rfield.flux_total_incident[0][i]/rfield.widflx(i),
				  1e-37));
				contin = rfield.flux[0][i] + rfield.ConEmitOut[0][i]*geometry.covgeo/rfield.widflx(i);
				y2[i] = MAX2((contin+(rfield.outlin[0][i]+rfield.outlin[0][i])/rfield.anu(i)*
					geometry.covgeo)* radius.r1r0sq,1e-37);
				y2[i] = log10(y2[i]);
			}
			else if( strcmp(plotCom.chPType[np],"REFL") == 0 )
			{
				/*    plot "reflected" continuum from last zone only */
				y[i] = log10(MAX2((rfield.ConEmitReflec[0][i]/rfield.widflx(i)+
				  rfield.reflin[0][i])*rfield.anu2(i),1e-37));
				y2[i] = y[i];
			}
			else if( strcmp(plotCom.chPType[np],"EMIT") == 0 )
			{
				/*    plot "emitted" continuum from both sides of cloud */
				y[i] = log10(MAX2(
					((rfield.ConEmitReflec[0][i]+rfield.ConEmitOut[0][i])/
					rfield.widflx(i)+
					(rfield.outlin[0][i]+rfield.outlin_noplot[i]+rfield.reflin[0][i])/rfield.anu(i) )*
					rfield.anu2(i),1e-37));
				y2[i] = y[i];
			}
			else if( strcmp(plotCom.chPType[np],"OUTW") == 0 )
			{
				/*    plot outward and attenuated incident continuum */
				chSymPlt1 = 'i';
				y[i] = log10(MAX2(rfield.flux[0][i]*opac.opacity_abs[i],
				  1e-37));
				strcpy( chSymPlt2, "o " );
				y2[i] = log10(MAX2((rfield.outlin[0][i]+rfield.outlin_noplot[i]+rfield.ConEmitOut[0][i])*
					opac.opacity_abs[i],1e-37));
			}
			else if( strcmp(plotCom.chPType[np],"DIFF") == 0 )
			{
				/*    plot "diffuse" continuum from last zone only */
				y[i] = log10(MAX2(rfield.ConEmitLocal[nzone][i]*rfield.anu2(i)/
				  rfield.widflx(i),1e-37));
				y2[i] = y[i];
			}
			else if( strcmp(plotCom.chPType[np],"CRAW") == 0 )
			{
				y[i] = log10(MAX2(rfield.flux_total_incident[0][i],1e-37));
				y2[i] = MAX2((rfield.flux[0][i]+
				  rfield.otscon[i]+rfield.otslin[i]+rfield.outlin[0][i]+rfield.outlin_noplot[i]+
				  rfield.ConEmitOut[0][i])*radius.r1r0sq,1e-37);
				y2[i] = log10(y2[i]);
			}

			if( y[i] > -36.9 )
			{
				ymin = MIN2(ymin,(double)y[i]);
			}

			if( y2[i] > -36.9 )
			{
				ymin2 = MIN2(ymin2,(double)y2[i]);
			}

			ymax = MAX2(ymax,(double)y[i]);
			ymax = MAX2(ymax,(double)y2[i]);
		}
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " PLOT called for the first time, YMAX, MIN=%10.2e%10.2e\n", 
		  ymax, ymin );
	}

	/* lower min by at most 5 dex below peak */
	ymin2 = MAX3(ymax-5.,-35.,ymin2);

	/* make sure there is room at the bottom */
	ymin = MIN3(ymin2,ymin,ymax-1.);

	/* emitted continuum is thermal, so goes to zero */
	if( strcmp(plotCom.chPType[np],"EMIT") == 0 )
	{
		ymin = MAX2(ymin,ymax-4.);
	}

	if( plotCom.lgPltTrace[np] )
	{
		fprintf( ioQQQ, " YMAX, MIN=%14.4e%14.4e Npnts=%4ld\n", 
		  ymax
		  , ymin, rfield.nflux );
	}
	strcpy( chXtitle, "Log(nu fnu) vs LOG(nu)" );

	chSym = chSymPlt1;

	pltr(rfield.anulog10ptr(),y,rfield.nflux,xmin,xmax,ymin,ymax,chSym,chXtitle
	  ,1,plotCom.lgPltTrace[np]);

	chSym = chSymPlt2[0];

	pltr(rfield.anulog10ptr(),y2,rfield.nflux,xmin,xmax,ymin,ymax,chSym,chXtitle
	  ,3,plotCom.lgPltTrace[np]);

	free( y );
	free( y2 );
	return;
}

/*pltmap generate plot of heating and cooling map */

STATIC void pltmap(
  long int np, 
  const char *chCall)
{
	char chXtitle[23];
	static bool lgTlkSav;
	char chSym;

	long int i;

	static double xmax, 
	  xmin, 
	  ymax, 
	  ymin;

	DEBUG_ENTRY( "pltmap()" );

	if( strcmp(chCall,"FIRST") == 0 )
	{ 
		return;
	}

	lgTlkSav = called.lgTalk;
	called.lgTalk = false;
	hcmap.lgMapBeingDone = true;
	hcmap.RangeMap[0] = exp10(plotCom.pltxmn[np]);
	hcmap.RangeMap[1] = exp10(plotCom.pltxmx[np]);
	map_do(ioQQQ, " map");
	called.lgTalk = lgTlkSav;

	for( i=0; i < hcmap.nmap; i++ )
	{
		hcmap.temap[i] = log10(hcmap.temap[i]);
	}

	xmin = MIN2(hcmap.temap[0],hcmap.temap[hcmap.nmap-1]);
	xmin = MAX2((double)plotCom.pltxmn[np],xmin);
	xmax = MAX2(hcmap.temap[0],hcmap.temap[hcmap.nmap-1]);
	xmax = MIN2(xmax,(double)plotCom.pltxmx[np]);

	if( plotCom.lgPltTrace[np] )
	{
		fprintf( ioQQQ, " xmin, xmax=%12.4e%12.4e nmap=%4ld\n", 
		  xmin, xmax, hcmap.nmap );
	}

	ymin = FLT_MAX;
	ymax = -FLT_MAX;

	for( i=0; i < hcmap.nmap; i++ )
	{
		if( hcmap.temap[i] > xmin && hcmap.temap[i] < xmax )
		{
			hcmap.hmap[i] = log10(MAX2(hcmap.hmap[i],1e-35));
			hcmap.cmap[i] = log10(MAX2(hcmap.cmap[i],1e-35));
			if( hcmap.cmap[i] > -34. )
			{
				ymin = MIN3(ymin,hcmap.hmap[i],hcmap.cmap[i]);
			}
			else
			{
				ymin = MIN2(ymin,hcmap.hmap[i]);
			}
			ymax = MAX3(ymax,hcmap.hmap[i],hcmap.cmap[i]);
		}
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " PLOT called for the first time, YMAX, MIN=%10.2e%10.2e\n", 
		  ymax, ymin );
	}

	if( plotCom.lgPltTrace[np] )
	{
		fprintf( ioQQQ, " YMAX, MIN=%14.4e%14.4e Npnts=%4ld\n", 
		  ymax, ymin, hcmap.nmap );
	}

	chSym = 'H';
	strcpy( chXtitle, "heating - cooling v te" );

	pltr(hcmap.temap,hcmap.hmap,hcmap.nmap,xmin,xmax,ymin,ymax,chSym,
	  chXtitle,1,plotCom.lgPltTrace[np]);

	chSym = 'C';

	pltr(hcmap.temap,hcmap.cmap,hcmap.nmap,xmin,xmax,ymin,ymax,chSym,
	  chXtitle,3,plotCom.lgPltTrace[np]);

	return;
}

/*pltopc generate plot of local gas opacity */
STATIC void pltopc(
  long int np, 
  const char *chCall)
{
	char chXtitle[23];
	char chSym;
	long int i;
	double arg1, 
	  arg2;
	static double xmax, 
	  xmin, 
	  ymax, 
	  ymin;
	static double *y/*[rfield.nflux_with_check]*/, 
	  *y2/*[rfield.nflux_with_check]*/;

	DEBUG_ENTRY( "pltopc()" );

	if( strcmp(chCall,"FIRST") == 0 )
	{ 
		return;
	}

	y = (double*)MALLOC((size_t)rfield.nflux_with_check*sizeof(double) );
	y2 = (double*)MALLOC((size_t)rfield.nflux_with_check*sizeof(double) );

	xmin = rfield.anulog10(0);
	xmin = MAX2((double)plotCom.pltxmn[np],xmin);
	xmax = rfield.anulog10(rfield.nflux-1);
	xmax = MIN2(xmax,(double)plotCom.pltxmx[np]);

	if( plotCom.lgPltTrace[np] )
	{
		fprintf( ioQQQ, " XMIN, XMAX=%12.4e%12.4e NFLUX=%4ld\n", 
		  xmin, xmax, rfield.nflux );
	}

	ymin = FLT_MAX;
	ymax = -FLT_MAX;

	for( i=0; i < rfield.nflux; i++ )
	{
		if( strcmp(plotCom.chPType[np],"OPAA") == 0 )
		{
			/*  absorption opacity */
			arg1 = opac.opacity_abs_savzon1[i];
			arg2 = opac.opacity_abs[i];
		}

		else if( strcmp(plotCom.chPType[np],"OPAS") == 0 )
		{
			/*  scattering opacity */
			arg1 = opac.opacity_sct_savzon1[i];
			arg2 = opac.opacity_sct[i];
		}

		else if( strcmp(plotCom.chPType[np],"OPAT") == 0 )
		{
			/*  total opacity */
			arg1 = opac.opacity_abs_savzon1[i] + opac.opacity_sct_savzon1[i];
			arg2 = opac.opacity_abs[i] + opac.opacity_sct[i];
		}

		else
		{
			/* this cannot happen since type was set to one of above */
			fprintf( ioQQQ, " pltopc type=%4.4s not known.  STOP\n", 
			  plotCom.chPType[np] );
			cdEXIT(EXIT_FAILURE);
		}

		y[i] = log10(MAX2(arg1/dense.gas_phase[ipHYDROGEN],1e-35));
		y2[i] = log10(MAX2(arg2/dense.gas_phase[ipHYDROGEN],1e-35));

		if( rfield.anulog10(i) > xmin && rfield.anulog10(i) < xmax )
		{
			ymin = MIN3(ymin,y[i],y2[i]);
			ymax = MAX3(ymax,y[i],y2[i]);
		}
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " PLOT called for the first time, YMAX, MIN=%10.2e%10.2e\n", 
		  ymax, ymin );
	}

	/* lower min by factor of 10 to show absorption in next plot */
	ymin = MAX2(ymin-1.,-35.);
	ymax += 1.;
	if( plotCom.lgPltTrace[np] )
	{
		fprintf( ioQQQ, " YMAX, MIN=%14.4e%14.4e Npnts=%4ld\n", 
		  ymax, ymin, rfield.nflux );
	}

	strcpy( chXtitle, "Log(opacity) vs log(n)" );

	chSym = '.';
	pltr(rfield.anulog10ptr(),y,rfield.nflux,xmin,xmax,ymin,ymax,chSym,chXtitle
	  ,1,plotCom.lgPltTrace[np]);

	chSym = 'o';
	pltr(rfield.anulog10ptr(),y2,rfield.nflux,xmin,xmax,ymin,ymax,chSym,chXtitle
	  ,3,plotCom.lgPltTrace[np]);

	free(y);
	free(y2);
	return;
}

/*pltr core plotting routine for generating line printer plots */

STATIC void pltr(
	/* the x-axis */
	const double x[], 
	/* the y-axi */
	const double y[], 
	/* number of points */
	long int npnts, 
	/* mins and maxs, log of min and max of x-axis */
	double xmin, 
	double xmax, 
	double ymin, 
	double ymax, 
	/* plot symbol */
	char chSym, 
	const char *chXtitle, 
	long int itim, 
	bool lgTrace)
{
	static char chPage[59][122];

	long int i, 
	  ix, 
	  iy, 
	  j, 
	  nc;

	/* the max number of decades we can plot */
	static const int NDECAD = 19;

	static long int jpnt[NDECAD], 
	  lowx, 
	  lx;

	static double xdec, 
	  xinc, 
	  ydown, 
	  yinc;

	/* this is log of smallestnumer in following set */
	const realnum xAxisMin = -9.f;

	static char chLab[NDECAD][5]={"1E-9","1E-8","1E-7","1E-6","1E-5",
		"1E-4",".001","0.01"," 0.1","  1 ",
	  " 10 "," 100","1000","1E4 ","1E5 ","1E6 ","1E7 ","1E8 ","1E9 "};

	DEBUG_ENTRY( "pltr()" );

	/* ITIM=1, first call, =2 intermediate calls, =3 for last call*/
	if( itim == 1 )
	{
		/* first call, set left border of plot and clear out array */
		for( i=1; i < IHI; i++ )
		{
			chPage[i][0] = 'l';
			for( j=1; j < IWID; j++ )
			{
				chPage[i][j] = ' ';
			}
		}

		/* centered label for plot */
		strcpy( chPage[1], "                        " );
		strcat( chPage[1],  chXtitle );
		strcat( chPage[1], input.chTitle );

		/* one dex increments in x and y marked special */
		i = 1;
		ydown = 0.;
		yinc = (realnum)(IHI-2)/(ymax - ymin);
		nc = 0;

		while( i <= IHI && nc < 200 )
		{
			chPage[i-1][1] = '-';
			ydown += yinc;
			i = (long)(ydown + 1);
			nc += 1;
		}

		/* bottom increments of plot */
		for( i=0; i < IWID; i++ )
		{
			chPage[IHI-1][i] = '-';
		}

		if( xmin < xAxisMin )
		{
			fprintf(ioQQQ," plts: xmin is less than min value in array\n");
			cdEXIT(EXIT_FAILURE);
		}
		/* LX is pointer to label for number in x-axis in chLab */
		if( xmin < 0. )
		{
			lx = (long)(4.999-fabs(xmin));
			/* lx is the offset within the array of x-axis values */
			/* >>chng 99 jun 11 change to allow any min value of x-axis */
			lx = (long)(fabs(xAxisMin)-0.001-fabs(xmin));
			lx = MAX2(0,lx);
			/* this is lowest decade on plot */
			xdec = -floor(fabs(xmin)+1e-5);
		}
		else
		{
			double aa;
			lx = (long)MAX2(0.,4.+xmin);
			/* lx is the offset within the array of x-axis values */
			lx = (long)MAX2(0.,4.+xmin);
			/* >>chng 99 jun 11 change to allow any min value of x-axis */
			aa = fabs(xAxisMin);
			lx = (long)MAX2(0., aa-1. + xmin );
			xdec = floor(xmin+1e-5);
		}

		lowx = lx + 1;
		xinc = (realnum)(IWID-1)/(xmax - xmin);
		i = (long)MAX2(1.,(xdec-xmin)*xinc+1.);
		nc = 0;

		while( i < IWID && nc < 100 )
		{
			chPage[IHI-2][i - 1] = 'l';

			/*  fix position of x labels */
			lx = MIN2(lx+1,NDECAD);

			/*  slight offset to center label */
			jpnt[lx-1] = MAX2(0,i-3);
			jpnt[lx-1] = MIN2((long)IWID-4,jpnt[lx-1]);
			xdec += 1.;
			i = (long)MAX2(1.,(xdec-xmin)*xinc+1.);
			nc += 1;
		}
	}

	/* everything falls down through here */
	/* now fill in data, symbol is chSym */
	for( i=0; i < npnts; i++ )
	{
		if( (double)x[i] > xmin && (double)x[i] < xmax )
		{
			iy = (long)(IHI - MAX2(y[i]-ymin,0.)*yinc);
			iy = MAX2(1,iy);
			ix = (long)((x[i] - xmin)*xinc + 1);

			if( lgTrace )
			{
				fprintf( ioQQQ, " x, y, ix, iy=%7.3f%7.3f%4ld%4ld\n", 
				  x[i], y[i], ix, iy );
			}
			chPage[iy-1][ix - 1] = chSym;
		}
	}

	if( itim == 3 )
	{
		/* make the plot */
		fprintf( ioQQQ, "1\n" );
		for( i=1; i < IHI; i++ )
		{
			fprintf( ioQQQ, "     %121.121s\n", chPage[i] );
		}

		/* now put on label for X-axis */
		for( i=0; i < IWID; i++ )
		{
			chPage[0][i] = ' ';
		}

		for( i=lowx-1; i < lx; i++ )
		{
			/* copy the four char of the numeric string */
			strncpy(chPage[0]+jpnt[i] , chLab[i+1] , 4);
		}
		fprintf( ioQQQ, "     %121.121s\n", chPage[0] );
	}
	return;
}
