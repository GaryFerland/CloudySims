/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MONOINTERP_H_
#define MONOINTERP_H_

// Internal utility functions
namespace {
	// Hermite polynomial basis functions
	inline double h00(double t)
	{
		return (1.0+2.0*t)*(1.0-t)*(1.0-t);
	}
	inline double h10(double t)
	{
		return t*(1.0-t)*(1.0-t);
	}
}

// Usage:
// Monointerp m(xvals,yvals, npt); -- Constructor
// y_interp = m(x_interp);         -- Interpolate value

void MI_SetupData(const vector<double>& d, const vector<double>& h, vector<double>& m_g, long n);

class Monointerp {
	const std::vector<double> m_x, m_y;
	std::vector<double> m_g;
public:
	// CONSTRUCTORS
	Monointerp(const double x[], const double y[], long n);
private:
	Monointerp(const Monointerp&);              // Not implemented
	Monointerp& operator= (const Monointerp&);  // Not implemented

	// MANIPULATORS -- None

	// ACCESSORS
public:
	double operator() (double xval) const;
};

#endif
