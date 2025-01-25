/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DEPTH_TABLE_H_
#define DEPTH_TABLE_H_

class DepthTable
{
public:
	/* lg is true if depth, false if radius to be used*/
	bool lgDepth;
	/**dist is log radius in cm, val is log value*/
	vector<double> dist;
	vector<double> val;

	/**number of values in above table */
	long int nvals() const
	{
		ASSERT( dist.size() == val.size() );
		return long(dist.size());
	}

	/**tabval, adapted from dense_tabden interpolate on table of points for density with dlaw table command, by K Volk 
		\param r0    current radius in cm
		\param depth current depth in cm
	*/
	double tabval( double r0, double depth) const;
	void clear()
	{
		lgDepth = false;
		dist.resize(0);
		val.resize(0);
	}
};

#endif // DEPTH_TABLE_H_
