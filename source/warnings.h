/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef WARNINGS_H_
#define WARNINGS_H_

/* warnings.h */

#include "module.h"

/** this is the limit to the number of warnings, cautions, and
 * notes that we can save */
static const int LIMWCN = 2000;

class t_warnings : public module {
public:
	const char *chName() const
	{
		return "warnings";
	}

	/**wcnint initialize stack or warnings, cautions, notes */
	void zero(void);
	void comment(t_warnings&) {}

	/** this are counters for the number of warnings,
	 * cautions, notes and surprises in the calculation*/
	long int nwarn, 
	  ncaun, 
	  nnote, 
	  nbang;

	/** a comment about the geometry after model stops */
	char chRgcln[2][INPUT_LINE_LENGTH];

	/** these are the strings that contain the warnings, cautions,
	 * and notes about the calculation */
	char chWarnln[LIMWCN][INPUT_LINE_LENGTH], 
	  chCaunln[LIMWCN][INPUT_LINE_LENGTH], 
	  chNoteln[LIMWCN][INPUT_LINE_LENGTH], 
	  chBangln[LIMWCN][INPUT_LINE_LENGTH];

	/** flags set if warnings or cautions present */
	bool lgWarngs, 
	  lgCautns;

	/**warnin enter warnings at the end of the calculations into large stack 
	\param *chLine
	*/
	void warnin(const char *chLine);
	
	/**notein enter a note about calculation into comment array 
	\param *chLine
	*/
	void notein(const char *chLine);
	
	/**bangin called by routine comment to enter surprise into comment stack 
		\param *chLine
	*/
	void bangin(const char *chLine);
	
	/**caunin called by comment to enter caution into comment stack 
		\param *chLine
	*/
	void caunin(const char *chLine);

	};
extern t_warnings warnings;


/**warnin enter warnings at the end of the calculations into large stack 
\param *chLine
*/
inline void warnin(const char *chLine)
{
	warnings.warnin(chLine);
}

/**notein enter a note about calculation into comment array 
\param *chLine
*/
inline void notein(const char *chLine)
{
	warnings.notein(chLine);
}

/**bangin called by routine comment to enter surprise into comment stack 
\param *chLine
*/
inline void bangin(const char *chLine)
{
	warnings.bangin(chLine);
}

/**caunin called by comment to enter caution into comment stack 
\param *chLine
*/
inline void caunin(const char *chLine)
{
	warnings.caunin(chLine);
}


#endif /* WARNINGS_H_ */
