/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HCMAP_H_
#define HCMAP_H_

/**punt produce map of heating-cooling space for specified zone 
\param *io
\param *chType
*/
void map_do(
			FILE *io, 
			const char *chType);

/** hcmap.h */
struct t_hcmap {

	/**parameters doing with map command, MapZone is zone number,
	 * -1 for very start, init to 9999 to prevent map*/
	long int MapZone;

	/** RangeMap is range, parsed from same command*/
	realnum RangeMap[2];

	/** nMapStep is number of steps, set with set nmap command*/
	long int nMapStep;

	/** logical flag indicating whether results of map were ok -
	 * set false if inflection points occurred */
	bool lgMapOK;

	/** this is set true when we are in process of doing a map */
	bool lgMapBeingDone;

	/** logical flag indicating whether a map has been done */
	bool lgMapDone;

	/** number of points actually done */
	long int nmap;

	/** saved temperatures for the map */
	double *temap;

	/** saved heating for the map */
	double *hmap; 

	/** saved cooling for the map */
	double *cmap;

	/** flag saying whether previous three vectors have had space allocated,
	 * this starts out zero and is set to first number needed when alloc happens, 
	 * and then never reset */
	long int nMapAlloc;

	};
extern t_hcmap hcmap;

#endif /* HCMAP_H_ */
