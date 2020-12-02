/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GENERIC_STATE_H_
#define GENERIC_STATE_H_

#include "quantumstate.h"

class molezone;
class genericState
{
public:
   genericState() : sp(NULL), val(NULL) {}
	explicit genericState(const molezone* sp) : sp(sp), val(NULL) {}
	explicit genericState(const qStateConstProxy& q) 
		: sp(NULL), q(q), val(NULL) {}
	explicit genericState(const char *valLabel, const double *val) 
		: sp(NULL), val(val), valLabel(valLabel) {}
	const molezone* sp;
	qStateConstProxy q;
	const double* val;
	string valLabel;
	bool associated() const;
	string label() const;
	string database() const;
};

double column(const genericState&);
double density(const genericState&);
double depart(const genericState&);
double energy(const genericState&);
double levels(const genericState&);

/** getLevelsGeneric -- get the requested species levels
 * \param chLabel	(in) species label
 * \param lgValidate	(in) flag to report errors in the level request
 * \param LevelList	(out) list of species levels
 */
const molezone *getLevelsGeneric(const char *chLabel, bool lgValidate, vector<long> &LevelList);

/** matchGeneric -- get a list of all species that match request
 * \param chLabel	(in) species label
 * \param lgValidate	(in) flag to report errors in the level request
 * \return 		(out) list of genericState species
 */
vector<genericState> matchGeneric(const char *chLabel, bool lgValidate);

#endif // GENERIC_STATE_H_
