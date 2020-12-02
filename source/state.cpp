/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*state_get_put get or save state - called by cloudy - job is either "get" or "put" */
/*state_do - worker to actually get or put the structure -
 * called by state_get_put below */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "struc.h"
#include "rfield.h"
#include "iterations.h"
#include "h2.h"
#include "dense.h"
#include "opacity.h"
#include "state.h"
#include "mole.h"

t_state state;
static bool lgGet;
static FILE *ioSTATE;

void t_state::zero()
{
	DEBUG_ENTRY( "t_state::zero()" );
	lgGet_state = false;
	lgPut_state = false;
	lgState_print = false;
}

/*state_do - worker to actually get or put the structure -
 * called by state_get_put below */
STATIC void state_do( void *pnt , size_t sizeof_pnt )
{
	size_t n;
	double sanity = 1.,
		chk_sanity;
	size_t sizeof_sanity =sizeof( sanity );

	DEBUG_ENTRY( "state_do()" );

	/* nothing to do if this is not positive */
	if( sizeof_pnt == 0 )
		return;

	if( lgGet )
	{
		/* get option - read in data */
		if( (n=fread( pnt , 1 , sizeof_pnt , ioSTATE )) - sizeof_pnt )
		{
			fprintf( ioQQQ, " state_do failed reading state file, wanted %lu got %lu\n",
				(unsigned long)sizeof_pnt ,
				(unsigned long)n);
			cdEXIT(EXIT_FAILURE);
		}
		/* perform sanity check */
		if( (n=fread( &chk_sanity , 1 , sizeof_sanity , ioSTATE )) - sizeof_sanity )
		{
			fprintf( ioQQQ, " state_do failed reading sanity par of state file, wanted %lu got %lu\n",
				(unsigned long)sizeof_sanity ,
				(unsigned long)n);
			cdEXIT(EXIT_FAILURE);
		}
		if( ! fp_equal( sanity, chk_sanity ) )
		{
			fprintf( ioQQQ, " state_do sanity fails in state file, wanted %g got %g\n",
				sanity ,
				chk_sanity);
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		/* put option - write out data */
		fwrite( pnt , 1 , sizeof_pnt , ioSTATE );
		/* write sanity check */
		fwrite( &sanity , 1 , sizeof_sanity , ioSTATE );
	}

	return;
}

/*state_get_put get or save state - called by cloudy - job is either "get" or "put" */
void state_get_put( const char chJob[] )
{
	long int ipISO , nelem , ipHi, i;
	
	DEBUG_ENTRY( "state_get_put()" );

	if( (strcmp( chJob , "get" ) == 0) )
	{
		lgGet = true;
		ioSTATE = open_data( state.chGetFilename, "rb", AS_LOCAL_ONLY );
	}
	else if( (strcmp( chJob , "put" ) == 0) )
	{
		lgGet = false;
		char chFilename[INPUT_LINE_LENGTH];
		if( !state.lgPutAll && iteration <= iterations.itermx )
		{
			/* not last iteration and do not want to save state for all iters, so simply quit */
			return;
		}

		/* get base of file name */
		strcpy( chFilename , state.chPutFilename );
		/* append iteration number if ALL keyword appears */
		if( state.lgPutAll )
		{
			char chIteration[INPUT_LINE_LENGTH];
			sprintf( chIteration , "_%li", iteration );
			strcat( chFilename , chIteration );
		}
		ioSTATE = open_data( chFilename, "wb" );
	}
	else
		TotalInsanity();

	if( state.lgState_print )
		fprintf(ioQQQ," Print state quantities, start iso seq \n");

	/* actually do the read / write */
	fixit("Wouldn't actually work, as these classes contain pointers");
#if 0
	iso_sp.state_do( ioSTATE, lgGet );
	ExtraLymanLines.state_do( ioSTATE, lgGet );
#endif

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* loop over all iso-electronic sequences */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				/* arrays are dim'd iso_sp[ipH_LIKE][nelem].numLevels_max+1 */
				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ++ipHi )
				{
					if( state.lgState_print )
					{
						fprintf(ioQQQ," start ISO ipISO= %li, nelem= %li, ipHi %li \n", 
							ipISO , 
							nelem ,
							ipHi);
						for( long n=0; n< ipHi; ++n )
						{
							fprintf(ioQQQ," ISO %li %li %li %li %.4e %.4e \n",
								ipISO , nelem , ipHi , n , 
								iso_sp[ipISO][nelem].trans(ipHi,n).Emis().TauIn() ,
								iso_sp[ipISO][nelem].trans(ipHi,n).Emis().TauTot() );
						}
						fprintf(ioQQQ," end ISO ipISO\n");
					}
				}

				if( state.lgState_print )
				{
					fprintf(ioQQQ," start Ext ipISO= %li, nelem= %li, got %li \n", 
						ipISO , 
						nelem ,
						iso_ctrl.nLyman_malloc[ipISO] );
				}
				if( state.lgState_print )
				{
					for( long n=2; n< iso_ctrl.nLyman_malloc[ipISO]; ++n )
					{
						TransitionList::iterator tr = ExtraLymanLines[ipISO][nelem].begin()+ipExtraLymanLines[ipISO][nelem][n];
						fprintf(ioQQQ," Ext %li %li %li %.4e %.4e \n",
							ipISO , nelem , n , 
							(*tr).Emis().TauIn() ,
							(*tr).Emis().TauTot() );
					}
					fprintf(ioQQQ," end Ext ipISO\n");
				}
			}
		}
	}

	fixit("Broken -- contains pointers...");
	//state_do( TauLine2.begin()+0, (size_t)nWindLine  *sizeof(transition) );
	if( state.lgState_print )
	{
		for( long n=0; n< nWindLine; ++n )
		{
			fprintf(ioQQQ," TauLine2 %li %.4e %.4e \n",
				n , 
				TauLine2[n].Emis().TauIn() ,
				TauLine2[n].Emis().TauTot() );
		}
	}

	// this implicitly assumes that a vector stores its data contiguously
	//state_do( UTALines.begin()+0, (size_t)nUTA   *sizeof(transition) );
	if( state.lgState_print )
	{
		for( size_t n=0; n< UTALines.size(); ++n )
		{
			fprintf(ioQQQ," UTALines %lu %.4e %.4e \n",
			   (unsigned long) n , 
				UTALines[n].Emis().TauIn() ,
				UTALines[n].Emis().TauTot() );
		}
	}

	//state_do( HFLines.begin()+0, (size_t)nHFLines   *sizeof(transition) );
	if( state.lgState_print )
	{
		for( size_t n=0; n< HFLines.size(); ++n )
		{
			fprintf(ioQQQ," HFLines %lu %.4e %.4e \n",
					  (unsigned long) n , 
				HFLines[n].Emis().TauIn() ,
				HFLines[n].Emis().TauTot() );
		}
	}

	for( i=0; i<2; ++i )
	{
		state_do( opac.TauAbsGeo[i]  , (size_t)rfield.nflux_with_check*sizeof(realnum) );
		if( state.lgState_print )
		{
			for( long n=0; n< rfield.nflux_with_check; ++n )
			{
				fprintf(ioQQQ," TauAbsGeo %li %li %.4e \n",
					i , n , 
					opac.TauAbsGeo[i][n] );
			}
		}

		state_do( opac.TauScatGeo[i] , (size_t)rfield.nflux_with_check*sizeof(realnum) );
		if( state.lgState_print )
		{
			for( long n=0; n< rfield.nflux_with_check; ++n )
			{
				fprintf(ioQQQ," TauScatGeo %li %li %.4e \n",
					i , n , 
					opac.TauAbsGeo[i][n] );
			}
		}

		state_do( opac.TauTotalGeo[i], (size_t)rfield.nflux_with_check*sizeof(realnum) );
		if( state.lgState_print )
		{
			for( long n=0; n< rfield.nflux_with_check; ++n )
			{
				fprintf(ioQQQ," TauTotalGeo %li %li %.4e \n",
					i , n , 
					opac.TauAbsGeo[i][n] );
			}
		}

	}

	/* the large H2 molecule, only if on */
	if( h2.lgEnabled )
	{
		//
		fixit("Wouldn't have worked -- transition contained pointers...");
		//state_do( &h2.sys.trans, h2.sys.trans.size() * sizeof(transition) );
		//state_do( &h2.states, h2.states.size() * sizeof( quantumState ) );
	}

	/* the struc variables */
	state_do( &struc.nzlim, sizeof(struc.nzlim ) );
	state_do( &struc.dr_ionfrac_limit, sizeof(struc.dr_ionfrac_limit ) );

	state_do( struc.testr, (size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.volstr, (size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.drad_x_fillfac, (size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.drad, (size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.histr, (size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.hiistr, (size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.ednstr, (size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.o3str, (size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.pressure,(size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.GasPressure ,(size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.pres_radiation_lines_curr ,(size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.hden ,(size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.DenParticles ,(size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.DenMass,(size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.depth,(size_t)(struc.nzlim)*sizeof(realnum ) );
	state_do( struc.xLyman_depth , (size_t)(struc.nzlim)*sizeof(realnum ) );

	state_do( struc.coolstr,(size_t)(struc.nzlim)*sizeof(double ) );
	state_do( struc.heatstr , (size_t)(struc.nzlim)*sizeof(double ) );

	for( long nz=0; nz<struc.nzlim; ++nz)
	{
		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			state_do( struc.xIonDense[nz][nelem] , (size_t)(LIMELM+1)*sizeof(realnum ) );
		}
	}
	for( long nz=0; nz<struc.nzlim; ++nz)
	{
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					state_do( struc.StatesElem[nz][nelem][nelem-ipISO] , (size_t)(iso_sp[ipISO][nelem].numLevels_local)*sizeof(realnum ) );
				}
			}
		}
	}

	for( long nz=0; nz<struc.nzlim; ++nz)
	{
		state_do( struc.molecules[nz] , (size_t)(mole_global.num_calc)*sizeof(realnum ) );
	}
	state_do( struc.H2_abund , (size_t)(struc.nzlim)*sizeof(realnum ) );

	for( long nz=0; nz<struc.nzlim; ++nz)
	{
		state_do( struc.gas_phase[nz] , (size_t)(LIMELM)*sizeof(realnum ) );
	}

	/*fprintf(ioQQQ,"DEBUG done\n");*/

	/* close the file */
	fclose( ioSTATE );
	return;
}
