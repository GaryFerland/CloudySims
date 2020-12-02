/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "parse_species.h"

#include "parser.h"
#include "mole.h"
#include "species.h"

// This file groups the handler for the species command with the
// routines which interpret the options.  It would be better if this
// were done in just one routine, to allow the code to warn on invalid
// options rather than just ignore them

void ParseSpecies(Parser &p)
{
	DEBUG_ENTRY( "ParseSpecies()" );
	string chString_quotes_original;
	bool lgQuotesFound = true;
	if (p.GetQuote(chString_quotes_original))
		lgQuotesFound = false;

	if( ! lgQuotesFound )
	{
		fprintf(ioQQQ,"Need to provide species name to SPECIES command.\n"
				  "The previous SPECIES command has been renamed DATABASE -- did you mean to use that?\n");
		cdEXIT(EXIT_FAILURE);
	}

	// Break up the remainder of the input line into options and (optionally) values
	// String or numeric values can be given without an '=', named options require '='
	Symbol s=p.getSymbol();
	while(1)
	{
		if (s.toktype == Symbol::EOSTAT)
		{
			break;
		}
		if ( s.toktype != Symbol::NAME )
		{
			fprintf(ioQQQ,"Value %s not understood\n",s.value.c_str());
			cdEXIT(EXIT_FAILURE);
		}
		Symbol t=p.getSymbol();
		if (t.toktype == Symbol::EOSTAT || t.toktype == Symbol::NAME)
		{
			mole_global.speciesProperties[chString_quotes_original].p.push_back(
				pair<string,count_ptr<Option> >(s.value,count_ptr<Option>(new Option(true))));
			s = t;
		}
		else if (t.toktype == Symbol::NUMBER)
		{
			long num = strtol(t.value.c_str(),NULL,10);
			mole_global.speciesProperties[chString_quotes_original].p.push_back(
				pair<string,count_ptr<Option> >(s.value,count_ptr<Option>(new Option(num))));
			s = p.getSymbol();
		}
		else if (t.toktype == Symbol::STRING)
		{
			mole_global.speciesProperties[chString_quotes_original].p.push_back(
				pair<string,count_ptr<Option> >(s.value,count_ptr<Option>(new Option(t.value,Option::QUOTED))));
			s = p.getSymbol();
		}
		else if (t.toktype == Symbol::OPERATOR && t.value == "=")
		{
			t = p.getSymbol();
			if (t.toktype == Symbol::NAME)
			{
				mole_global.speciesProperties[chString_quotes_original].p.push_back(
				pair<string,count_ptr<Option> >(s.value,count_ptr<Option>(new Option(t.value,Option::NOTQUOTED))));
			}
			else if (t.toktype == Symbol::STRING)
			{
				mole_global.speciesProperties[chString_quotes_original].p.push_back(
				pair<string,count_ptr<Option> >(s.value,count_ptr<Option>(new Option(t.value,Option::QUOTED))));
			}
			else if (t.toktype == Symbol::NUMBER)
			{
				long num = strtol(t.value.c_str(),NULL,10);
				mole_global.speciesProperties[chString_quotes_original].p.push_back(
				pair<string,count_ptr<Option> >(s.value,count_ptr<Option>(new Option(num))));
			}
			else
			{
				fprintf(ioQQQ,"Option '%s' value '%s' not understood\n",s.value.c_str(),t.value.c_str());
				fprintf(ioQQQ,"WARNING: This usage of the SPECIES command was not understood.\n"
						  "The previous SPECIES command has been renamed DATABASE -- did you mean to use that?\n");
				cdEXIT(EXIT_FAILURE);
			}
			s = p.getSymbol();
		}
		else
		{
			fprintf(ioQQQ,"Option '%s' not understood\n",s.value.c_str());
			fprintf(ioQQQ,"WARNING: This usage of the SPECIES command was not understood.\n"
					  "The previous SPECIES command has been renamed DATABASE -- did you mean to use that?\n");
			cdEXIT(EXIT_FAILURE);
		}
	}
}

void setProperties(species& sp)
{
	DEBUG_ENTRY("setProperties()");
	char chLabelChemical[CHARS_SPECIES];
	if( sp.lgMolecular )
	{
		strcpy( chLabelChemical, sp.chLabel );
	}
	else
	{
		long nelem = 0, IonStg;
		parsespect(sp.chLabel,nelem,IonStg);
		makeChemical(chLabelChemical, nelem, IonStg-1);
	}

	map<string,Properties >::iterator props = 
		mole_global.speciesProperties.find(chLabelChemical);
	if (props == mole_global.speciesProperties.end())
		return;

	props->second.setDone();

	for(vector<pair<string,count_ptr<Option> > >::iterator prop=props->second.p.begin();
		 prop!=props->second.p.end();++prop)
	{
		if (prop->first == "LEVELS")
		{
			if (prop->second->opttype == Option::LONG)
			{
				if (prop->second->i > sp.numLevels_max)
				{
					fprintf(ioQQQ,"Species '%s', %ld levels requested, only %ld available.  Using all.\n",
							  chLabelChemical,prop->second->i,
							  sp.numLevels_max);
					sp.setLevels = sp.numLevels_max;
				}
				else
				{
					sp.setLevels = prop->second->i;
				}
			}
			else if (prop->second->opttype == Option::OPTION && prop->second->s == "ALL")
			{
				sp.setLevels = sp.numLevels_max;
			}
			else
			{
				fprintf(ioQQQ,"Incorrect type for 'LEVELS' option to species\n"
						  "Expecting 'LEVELS <number>', 'LEVELS=<number>' or 'LEVELS=ALL'\n");
			}
		}
		else if (prop->first == "LTE")
		{
			if (prop->second->opttype == Option::BOOL)
			{
				sp.lgLTE = prop->second->l;
			}
			else
			{
				fprintf(ioQQQ,"Incorrect type for 'LTE' option to species\n");
			}
		}
		else if (prop->first == "OFF")
		{
			; // Do nothing, already handled
		}
		else
		{
			fprintf(ioQQQ,"Option '%s' not understood for species '%s'\n",
					  prop->first.c_str(),sp.chLabel);
			cdEXIT(EXIT_FAILURE);
		}
	}
}

void speciesCheck()
{
	for (map<string,Properties>::iterator props = 
			  mole_global.speciesProperties.begin(); 
		  props != mole_global.speciesProperties.end();
		  ++props)
	{
		if (! props->second.isDone())
		{
			fprintf(ioQQQ,"\n\nWarning: SPECIES \"%s\" command has no effect since species not found.\n"
				"Is species inactive or misspelt?  Names of species are case-sensitive.  Consult species list in docs/SpeciesLabels.txt.\n\n",
				props->first.c_str());
		}
	}
}

bool speciesOff(const string& label)
{
	map<string,Properties>::iterator props = 
		mole_global.speciesProperties.find(label);
	if (props != mole_global.speciesProperties.end())
	{
		props->second.setDone();

		for(vector<pair<string,count_ptr<Option> > >::iterator prop=props->second.p.begin();
			 prop!=props->second.p.end();++prop)
		{
			if (prop->first == "OFF")
			{
				if (prop->second->opttype != Option::BOOL)
				{
					fprintf(ioQQQ,"Incorrect type for 'OFF' option to species\n");
				}
				
				if (prop->second->l || prop->second->opttype != Option::BOOL)
				{
					return true;
				}
			}
		}
	}
	return false;
}
