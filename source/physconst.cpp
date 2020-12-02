/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"

static int register_physconst(const char* name, double value);

// First define all the constants
#include "physconst.h"

// Load template file with alternate action
#define NEW_CONSTANT(NAME, VALUE)						\
	static int UNUSED physconst_##NAME = register_physconst(#NAME,VALUE)
#include "physconst_template.h"
#undef NEW_CONSTANT

static int register_physconst(const char* name, double value)
{ 
	t_physconst::Inst().add_constant( name, value );
	return 1;
}

t_physconst::t_physconst()
{
#ifndef HAVE_CONSTEXPR
//	dprintf( ioQQQ, "MILNE_CONST %.16e SAHA %.16e\n", SPEEDLIGHT*sqrt(pow3(FINE_STRUCTURE2)*pow3(TE1RYD)/PI),
//		 sqrt(pow3(HION_LTE_POP)) );
	if (!fp_equal(MILNE_CONST, 
					  SPEEDLIGHT*sqrt(pow3(FINE_STRUCTURE2)*pow3(TE1RYD)/PI) ))
	{
		fprintf(stderr," PROBLEM - Failed consistency check for MILNE_CONSTANT\n");
		exit(1);
	}
	if (!fp_equal(SAHA, 
					  sqrt(pow3(HION_LTE_POP)))) 
	{
		fprintf(stderr," PROBLEM - Failed consistency check for SAHA\n");
		exit(1);
	}
#endif
}

//
// the macro NEW_CONSTANT calls a function to set a static global
// integer for every constant.  the constructor of that class will
// populate the list of constants in the t_physconst singleton by
// calling the add_constant() method.
//
// these integers are not used anywhere else: need to check if the
// linker would optimize this module away if the routine below
// wouldn't exist and be used in parse_print.cpp.
//
// this routine can be removed if this file is merged with other files that
// instantiate global data.
//
void prt_phys_constants(FILE* io)
{
	t_physconst::Inst().prt_constants(io);
}
