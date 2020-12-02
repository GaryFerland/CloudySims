/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "lines.h"
#include "prt.h"
#include "lines_service.h"
#include "version.h"
#include "service.h"

t_LineSave LineSave;
/* these are the definitions of the line save arrays in lines.h */

void t_LineSave::zero()
{
	DEBUG_ENTRY( "t_LineSave::zero()" );

	/* index within the line in the line stack 
	 * default is Hbeta total - the third line in the stack
	 * 0th is a zero for sanity, 1st is unit, 2nd is a comment */
	/* >>chng 02 apr 22 from 2 to 3 since added unit at 1 */
	/* >>chng 06 mar 11, from 3 to -1 will now set to "H  1" 4861 */
	ipNormWavL = -1;
	WavLNorm = 4861.33f;
	lgNormSet = false;
	sig_figs = sig_figs_max;

	/* the label for the normalization line */
	strcpy( chNormLab, "    " );

	/* The scale factor for the normalization line */
	ScaleNormLine = 1.;

}

void LinSv::prt(FILE* ioPUN) const
{
	DEBUG_ENTRY( "LinSv::prt()" );
	fprintf(ioPUN,"%s",label().c_str());
}

string LinSv::label() const
{
	DEBUG_ENTRY( "LinSv::label()" );
	string val = chALab();
	val.resize( NCHLAB-1, ' ' );
	val += " ";
	char buf[100];
	sprt_wl(buf, wavelength());
	val += buf;
	return val;
}

string LinSv::biglabel() const
{
	DEBUG_ENTRY( "LinSv::biglabel()" );
	string val = label();
	if (m_tr.associated())
	{
		val += " ";
		val += "[";
		val += (*m_tr.Hi()).chConfig();
		val += "->";
		val += (*m_tr.Lo()).chConfig();
		val += "]";
	}
	return val;
}

bool LinSv::isCat(const char *s) const
{
	DEBUG_ENTRY( "LinSv::isCat()" );
	return strncmp(chALab(),s,strlen(s)) == 0;
}

void LinSv::chALabSet(const char *that)
{
	/* check that null is correct, string overruns have 
	 * been a problem in the past */
	ASSERT( (int) strlen( that )< NCHLAB );
	strncpy(m_chALab,that,NCHLAB-1);
	m_chALab[NCHLAB-1] = '\0';
	trimTrailingWhiteSpace( m_chALab );

	strncpy(m_chCLab,that,NCHLAB-1);
	m_chCLab[NCHLAB-1] = '\0';
	trimTrailingWhiteSpace( m_chCLab );
	caps(m_chCLab);
}

void LinSv::init(long index, char chSumTyp, const char *chComment, const char *label,
					  const TransitionProxy& tr)
{
	DEBUG_ENTRY( "LinSv::init()" );
	/* first call to stuff lines in array, confirm that label is one of
	 * the four correct ones */
	ASSERT( (chSumTyp == 'c') || (chSumTyp == 'h') || (chSumTyp == 'i') || (chSumTyp == 'r')  || (chSumTyp == 't') );
	ASSERT(!m_tr.associated());
	ASSERT(index >= 0);
	m_index = index;
	/* then save it into array */
	m_chSumTyp = chSumTyp;
	emslinZero();
	m_chComment = chComment;
	chALabSet( label );
	if (isCat("####"))
	{
		m_type = SEPARATOR;
	}
	else if (isCat("Unit"))
	{
		m_type = UNIT;
	}
	else if (isCat("UntD"))
	{
		m_type = UNITD;
	}
	else if (isCat("Inwd"))
	{
		m_type = INWARD;
	}
	else if (isCat("InwC"))
	{
		m_type = INWARDCONTINUUM;
	}
	else if (isCat("InwT"))
	{
		m_type = INWARDTOTAL;
	}
	else if (isCat("Coll"))
	{
		m_type = COLLISIONAL;
	}
	else if (isCat("Pump"))
	{
		m_type = PUMP;
	}
	else if (isCat("Heat"))
	{
		m_type = HEAT;
	}
	else if (isCat("Ca A"))
	{
		m_type = CASEA;
	}
	else if (isCat("Ca B"))
	{
		m_type = CASEB;
	}
	else if (isCat("nInu"))
	{
		m_type = NINU;
	}
	else if (isCat("nFnu"))
	{
		m_type = NFNU;
	}
	else if (isCat("Pho+"))
	{
		m_type = PHOPLUS;
	}
	else if (isCat("Pcon"))
	{
		m_type = PCON;
	}
	else if (isCat("Q(H)"))
	{
		m_type = QH;
	}
	else
	{
		m_type = DEFAULT;
	}

	SumLineZero();
	m_tr = tr;
}

void LinSv::addComponentID(long id)
{
	m_component.push_back(id);
	if (m_component.size() == 1)
		m_chComment += ": ";
	else
		m_chComment += "+";
	if (id > 0)
	{
		m_chComment += "\""+LineSave.lines[id].label()+"\"";
	}
	else
	{
		m_chComment += "N/A";
	}
}
void LinSv::addComponent(const char* species,const double wavelength1)
{
	if ( LineSave.ipass == 0 )
	{
		long id = LineSave.findline(species,wlAirVac(wavelength1));
		if (id <= 0)
		{
			fprintf( ioQQQ, "ERROR: A component to line blend \"%s\" %.3f was not identified.\n",
				species, wavelength1 );
			cdEXIT( EXIT_FAILURE );
		}
		addComponentID(id);
	}
}

// Automatically generate blend for specified species, at wavelength +/- width
void LinSv::makeBlend(const char* chLabel,const double wavelength1, const double width)
{
	DEBUG_ENTRY("LinSv::makeBlend()");
	if ( LineSave.ipass == 0 )
	{
		realnum wlo = wlAirVac(wavelength1-width),
			whi = wlAirVac(wavelength1+width);
			
		/* check that chLabel[4] is null - supposed to be 4 char + end */
		if( strlen(chLabel) > NCHLAB-1 )
		{
			fprintf( ioQQQ, " makeBlend called with insane species \"%s\", must be %d or less characters long.\n",
						chLabel, NCHLAB-1 );
			cdEXIT(EXIT_FAILURE);
		}
		
		char chCARD[INPUT_LINE_LENGTH];
		strcpy( chCARD, chLabel );
		
		/* make sure chLabel is all caps */
		caps(chCARD);/* convert to caps */
		
		/* this replaces tabs with spaces. */
		/* \todo	2 keep this in, do it elsewhere, just warn and bail? */
		for( char *s=chCARD; *s != '\0'; ++s )
		{
			if( *s == '\t' )
			{
				*s = ' ';
			}
		}

		for( long j=1; j < LineSave.nsum; j++ )
		{
			/* use pre-capitalized version of label to be like input chLineLabel */
			const char *chCaps = LineSave.lines[j].chCLab();

			if (LineSave.wavelength(j) >= wlo && 
				 LineSave.wavelength(j) <= whi &&
				 strcmp(chCaps,chCARD) == 0)
			{
				addComponentID(j);
				fprintf( ioQQQ,"Adding \"%s\" to blend\n", 
							LineSave.lines[j].label().c_str() );
			}
		}
	}
}
string LinSv::chComment() const
{
	return m_chComment;
}

static bool wavelength_compare (long a, long b)
{
	LinSv* a1 = &LineSave.lines[a];
	LinSv* b1 = &LineSave.lines[b];
	/* comparison is b-a so we get inverse wavelength order (increasing energy order) */
	if( b1->wavelength() < a1->wavelength() )
		return true;
	else 
		return false;
}

static bool wavelength_compare_realnum (size_t a, realnum wavelength)
{
	/* comparison is b-a so we get inverse wavelength order (increasing energy order) */
	if( wavelength < LineSave.wavelength(a) )
		return true;
	else 
		return false;
}

void t_LineSave::setSortWL()
{
	DEBUG_ENTRY( "t_LineSave::setSortWL()" );
	SortWL.resize((unsigned)nsum);
	for (long nlin=0; nlin < nsum; ++nlin)
	{
		SortWL[nlin] = nlin;
	}
	stable_sort(SortWL.begin(), SortWL.end(), wavelength_compare);
}

long t_LineSave::findline(const char *chLabel, realnum wavelength1)
{
	DEBUG_ENTRY( "t_LineSave::findline()" );
	
	ASSERT( nsum >= 0);
	if (nsum == 0)
		return -1;

	bool lgDEBUG = false;
	
	if( strlen(chLabel) > NCHLAB-1 )
	{
		fprintf( ioQQQ, " findline called with insane chLabel (between quotes) \"%s\", must be no more than %d characters long.\n",
					chLabel, NCHLAB-1 );
		return -2;
	}
 
	char chCARD[INPUT_LINE_LENGTH];
	strcpy( chCARD, chLabel );

	/* make sure chLabel is all caps */
	caps(chCARD);/* convert to caps */
	trimTrailingWhiteSpace( chCARD );

	/* this replaces tabs with spaces. */
	/* \todo	2 keep this in, do it elsewhere, just warn and bail? */
	for( char *s=chCARD; *s != '\0'; ++s )
	{
		if( *s == '\t' )
		{
			*s = ' ';
		}
	}

	/* get the error associated with specified significant figures; add
	 * FLT_EPSILON*wavelength to broaden bounds enough to allow for
	 * cancellation error
	 */
	realnum errorwave = WavlenErrorGet( wavelength1, LineSave.sig_figs ) + FLT_EPSILON*wavelength1, 
		smallest_error=BIGFLOAT,
		smallest_error_w_correct_label=BIGFLOAT;

	// Possibly more 'user friendly' line search mode -- not active
	// NB falls through to previous implementation if ipass < 1, as
	// lines will not yet be sorted
	if (ipass == 1)
	{
		char wlbuf[INPUT_LINE_LENGTH];
		// Need to search for lines within allowed band, and plausible confusions

		// Find position in list of lines
		vector<size_t>::iterator first =
			lower_bound(SortWL.begin(),SortWL.end(),
							wavelength1+errorwave, wavelength_compare_realnum);

		// first is now first line below upper limit
		if (first == SortWL.end())
		{
			sprt_wl(wlbuf,wavelength1);
			fprintf(ioQQQ,"Didn't find anything at %s\n",wlbuf);
			cdEXIT(EXIT_FAILURE);
		}

		// look for first line below lower limit -- may be the same as first if there are no matches
		vector<size_t>::iterator second;
		for(second=first; second != SortWL.end(); ++second)
		{
			if (wavelength(*second) < wavelength1-errorwave)
				break;
		}

		vector<size_t>::iterator found = SortWL.end();
		int nmatch=0;
		realnum dbest = 0.;
		for (vector<size_t>::iterator pos = first; pos != second; ++pos)
		{
			if ( strcmp(lines[*pos].chCLab(),chCARD) == 0 )
			{
				++nmatch;
				realnum dwl = wavelength(*pos)-wavelength1;
				if(nmatch >= 2 && !t_version::Inst().lgReleaseBranch
					&& !t_version::Inst().lgRelease )
				{
					if (nmatch == 2)
					{
						sprt_wl(wlbuf,wavelength1);
						fprintf(ioQQQ,"WARNING: multiple matching lines found in search for \"%s\" %s\n",
								  chLabel,wlbuf);
						fprintf(ioQQQ,"WARNING: match 1 is \"%s\" (dwl=%gA)\n",
								  lines[*found].biglabel().c_str(),wavelength(*found)-wavelength1);
					}
					fprintf(ioQQQ,"WARNING: match %d is \"%s\" (dwl=%gA)\n",
							  nmatch, lines[*pos].biglabel().c_str(),dwl);
				}
				if ( found == SortWL.end() )
				{
					found = pos;
					dbest = fabs(wavelength(*pos)-wavelength1);
				}
				else if ( fabs(dwl) < dbest )
				{
					found = pos;
					dbest = fabs(dwl);
				}
			}
		}
		if ( found != SortWL.end())
		{
			if (0)
				fprintf(ioQQQ,"Found %s\n", lines[*found].label().c_str());
			return *found;
		}

		sprt_wl(wlbuf,wavelength1);
		fprintf(ioQQQ,"WARNING: no exact matching lines found for \"%s\" %s\n",chLabel,wlbuf);
		for (vector<size_t>::iterator pos = first; pos != second; ++pos)
		{
			fprintf(ioQQQ,"WARNING: Line with incorrect label found close \"%s\"\n",
					  lines[*pos].label().c_str());
		}
		// Haven't found a match with correct label
		vector<size_t>::iterator best = SortWL.end();
		realnum besterror = 0.;
		for (;;)
		{
			realnum errordown = wavelength(*(first-1))-wavelength1;
			realnum errorup = wavelength1- (second == SortWL.end() ? 0.0 : wavelength(*second)) ;
			realnum error = 0.;
			vector<size_t>::iterator next;
			if ( errordown < errorup || second == SortWL.end())
			{
				error = errordown;
				--first;
				next = first;
			}
			else
			{
				error = errorup;
				next = second;
				++second;
			}
			if ( strcmp(lines[*next].chCLab(),chCARD) == 0 )
			{
				if (best == SortWL.end())
				{
					best = next;
					besterror = error;
					fprintf(ioQQQ,"Taking best match as \"%s\"\n",lines[*next].label().c_str());
				}
				else
				{
					fprintf(ioQQQ,"Possible ambiguity with \"%s\"\n",lines[*next].label().c_str());
				}
			}
			// Assume this is clearly unambiguous
			if (best != SortWL.end() && error > 100.*besterror)
				break;
			// Assume this is clearly unmatched
			if (error > 0.01*wavelength1)
				break;
		}
		if (best != SortWL.end())
			return *best;

		sprt_wl(wlbuf,wavelength1);
		fprintf(ioQQQ,"PROBLEM: no matching line found in search for \"%s\" %s\n",chLabel,wlbuf);
		cdEXIT(EXIT_FAILURE);
	}

	// Falls through to previous implementation if ipass < 1.  Can
	// be more restrictive with match handling, as this will be from
	// a request internal to the code infrastructure, rather than
	// user data

	long int j, index_of_closest=LONG_MIN,
		index_of_closest_w_correct_label=-1;
	
	for( j=1; j < nsum; j++ )
	{
		realnum current_error = (realnum)fabs(wavelength(j)-wavelength1);
		/* use pre-capitalized version of label to be like input chLineLabel */
		const char *chCaps = lines[j].chCLab();

		if( current_error < smallest_error )
		{
			index_of_closest = j;
			smallest_error = current_error;
		}
			
		if( current_error < smallest_error_w_correct_label &&
			 (strcmp(chCaps,chCARD) == 0) )
		{
			index_of_closest_w_correct_label = j;
			smallest_error_w_correct_label = current_error;
		}

		/* check wavelength and chLabel for a match */
		/*if( fabs(lines[j].wavelength- wavelength)/MAX2(DELTA,wavelength)<errorwave && 
			strcmp(chCaps,chCARD) == 0 ) */
		if( lgDEBUG && (current_error <= errorwave || 
			  fp_equal( wavelength1 + errorwave, wavelength(j) ) ||
			  fp_equal( wavelength1 - errorwave, wavelength(j) ))  
			 && strcmp(chCaps,chCARD) == 0 )
		{
			/* match, so set emiss to emissivity in line */
			/* and announce success by returning line index within stack */
			printf("Matched %s %15.8g %ld %18.11g %s\n",
					 chCaps,wavelength1,j,wavelength(j),
					 lines[j].biglabel().c_str());
		}
	}


	// Didn't find line, handle error
	if( index_of_closest_w_correct_label == -1 ||
		smallest_error_w_correct_label > errorwave )
	{
		/* >>chng 05 dec 21, report closest line if we did not find exact match, note that
		 * exact match returns above, where we will return negative number of lines in stack */
		fprintf( ioQQQ," PROBLEM findline did not find line " );
		prt_line_err( ioQQQ, chCARD, wavelength1 );
		if( index_of_closest >= 0 )
		{
			fprintf( ioQQQ,"  The closest line (any label) was   \"%s\"\n", 
						lines[index_of_closest].label().c_str() );
			if( index_of_closest_w_correct_label >= 0 )
			{
				fprintf( ioQQQ,"  The closest with correct label was \"%s\"\n", 
							lines[index_of_closest_w_correct_label].label().c_str() );
				fprintf( ioQQQ,"  Error was %15.8g vs. tolerance %15.8g\n", 
							smallest_error_w_correct_label, errorwave );
			}
			else
				fprintf( ioQQQ,"\n  No line found with label \"%s\".\n", chCARD );
			fprintf( ioQQQ,"\n" );
		}
		else
		{
			fprintf( ioQQQ,".\n PROBLEM No close line was found\n" );
			TotalInsanity();
		}
		return -3;
	}

	if (lgDEBUG)
		fprintf(ioQQQ,"Identified %ld\n",index_of_closest_w_correct_label);


	return index_of_closest_w_correct_label;
}

/*************************************************************************
 *
 * cdEmis obtain the local emissivity for a line with known index
 *
 ************************************************************************/
void cdEmis(
	const LinSv* line,
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// intrinsic or emergent
	bool lgEmergent )
{
	DEBUG_ENTRY( "cdEmis()" );

	if (line)
		*emiss = line->emslin(lgEmergent);
	else
		*emiss = 0.;
}

