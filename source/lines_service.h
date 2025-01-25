/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef LINES_SERVICE_H_
#define LINES_SERVICE_H_

class LinSv;

/** linadd - Add line to line stack, but don't transfer it
 *
 * This version of the function must be used with the 'Inwd' parts of
 * transferred lines.  See Issue #491 for problems that led to its inception.
 *
 * \param xInten [in] - local line emissivity, no filling factor
 * \param xIntenIsoBkg [in] - as above, but corrected for isotropic continua
 * \param wavelength [in] - line wavelength
 * \param chLab [in] - line label, e.g., 'Fe26'
 * \param chInfo [in] - line type: 'i' info, 'c' cooling, 'h' heating, 'r' recomb
 * \param chComment [in] - line comment, shown on output of 'save line labels'
 *
 * \returns A pointer to the newly added line on the line stack.
 */
LinSv *linadd(
  double xInten,
  double xIntenIsoBkg,
  t_wavl wavelength,
  const char *chLab,
  char chInfo,
  const char *chComment );

/** linadd - Add line to line stack, but don't transfer it
 *
 * The original version of the function, uses the local emissivity as the value
 * corrected for isotropic continua.  For limitations to its use, see Issue #491.
 *
 * \param xInten [in] - local line emissivity, no filling factor
 * \param wavelength [in] - line wavelength
 * \param chLab [in] - line label, e.g., 'Fe26'
 * \param chInfo [in] - line type: 'i' info, 'c' cooling, 'h' heating, 'r' recomb
 * \param chComment [in] - line comment, shown on output of 'save line labels'
 *
 * \returns A pointer to the newly added line on the line stack.
 */
LinSv *linadd(
  double xInten,
  t_wavl wavelength,
  const char *chLab,
  char chInfo ,	
  const char *chComment );

/*outline_base - adds line photons to reflin and outlin */
void outline_base(double dampXvel, double damp, bool lgTransStackLine, long int ip, double phots, realnum inwd,
						double nonScatteredFraction);

/*outline_base_bin - adds line photons to bins of reflin and outlin */
void outline_base_bin(bool lgTransStackLine, long int ip, double phots, realnum inwd,
						double nonScatteredFraction);

/** put forbidden line into stack, using index derived below 
\param xInten - local emissivity per unit vol
\param wavelength wavelength Angstroms
\param *chLab string label for ion
\param ipnt offset of line in continuum mesh
\param chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
\param lgOutToo should line be included in outward beam?
\param *chComment string explaining line 
*/
void lindst(double xInten,
  t_wavl wavelength, 
  const char *chLab, 
  long int ipnt, 
  char chInfo, 
  bool lgOutToo,
  const char *chComment);

/** put forbidden line into stack, using index derived below
\param dampXvel - damping constant times Doppler velocity
\param damp - damping constant
\param xInten - local emissivity per unit vol
\param wavelength wavelength Angstroms
\param *chLab string label for ion
\param ipnt offset of line in continuum mesh
\param chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
\param lgOutToo should line be included in outward beam?
\param *chComment string explaining line
*/
void lindst(double dampXvel,
  double damp,
  double xInten,
  t_wavl wavelength,
  const char *chLab,
  long int ipnt,
  char chInfo,
  bool lgOutToo,
  const char *chComment);

/** put forbidden line into stack, using index derived below 
\param t -- pointer to transition data
\param wavelength wavelength Angstroms
\param *chLab string label for ion
\param ipnt offset of line in continuum mesh
\param chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
\param lgOutToo should line be included in outward beam?
\param *chComment string explaining line 
*/
class TransitionProxy;
class ExtraInten;
void lindst(
  const TransitionProxy &t,
  const ExtraInten &extra,
  const char *chLab, 
  char chInfo, 
  bool lgOutToo,
  const char *chComment);

/** absorption due to continuous opacity 
\param emissivity [erg cm-3 s-1] in inward direction 
\param emissivity [erg cm-3 s-1] in outward direction 
\param array index for continuum frequency 
*/
double emergent_line( 
	 /* emissivity [erg cm-3 s-1] in inward direction */
	 double emissivity_in , 
	 /* emissivity [erg cm-3 s-1] in outward direction */
	 double emissivity_out , 
	 /* array index for continuum frequency */
	 long int ipCont );

/**PntForLine generate pointer for forbidden line 
\param wavelength wavelength of line in Angstroms 
\param *chLabel label for the line
\param *ipnt this is array index on the f, not c scale,
   for the continuum cell holding the line
*/
void PntForLine(t_wavl wavelength, 
  const char *chLabel, 
  long int *ipnt);

/**GetGF convert Einstein A into oscillator strength 
\param eina
\param enercm
\param gup
*/
double GetGF(double eina, 
	  double enercm, 
	  double gup);

/** S2Aul convert line strength S into transition probability Aul
\param S line strength
\param waveAng vacuum wavelength in Angstrom
\param gup statistical weight of the upper level
\param transType transition type, "E1", "M1", "E2", etc.
*/
double S2Aul(double S,
	     double waveAng,
	     double gup,
	     const string& transType);

/**eina convert a gf into an Einstein A 
\param gf
\param enercm
\param gup
*/
double eina(double gf, 
	  double enercm, 
	  double gup);

/**abscf convert gf into absorption coefficient 
\param gf
\param enercm
\param gl
*/
double abscf(double gf, 
	  double enercm, 
	  double gl);

/** setting true will use low-density Lyman branching ratios */
#define LOWDEN_LYMAN 0


/**WavlenErrorGet - given the real wavelength in A for a line
 * routine will find the error expected between the real 
 * wavelength and the wavelength printed in the output, with 4 sig figs,
  \param wavelength
  \return function returns difference between exact and 4 sig fig wl, so 
  we have found correct line is fabs(d wl) < return 
  */
realnum WavlenErrorGet( realnum wavelength, long sig_figs );

/** convert down coll rate back into electron cs in case other parts of code need this for reference 
\param gHi - stat weight of upper level 
\param rate - deexcitation rate, units s-1 
*/
double ConvRate2CS( realnum gHi , realnum rate );

/** convert collisional deexcitation cross section for into collision strength 
\param CrsSectCM2 - the cross section
\param gLo - statistical weight of lower level of transition
\param E_ProjectileRyd - initial projectile energy in Rydbergs
\param reduced_mass_grams - reduced mass MpMt/(Mp+Mt) of projectile-target system
*/
double ConvCrossSect2CollStr( double CrsSectCM2, double gLo, double E_ProjectileRyd, double reduced_mass_grams );

/**totlin sum total intensity of cooling, recombination, or intensity lines 
\param chInfo chInfor is 1 char,<BR>
	'i' information, <BR>
	'r' recombination or <BR>
	'c' collision
*/
double totlin(
	int chInfo);


/**FndLineHt search through line heat arrays to find the strongest heat source 
\param *level
*/
const TransitionProxy FndLineHt(long int *level);

/**set_xIntensity: compute gross and net number of emitted line photons */
void set_xIntensity( const TransitionProxy &t );

/**wn2angVac convert energy in wavenumbers to vacuum wavelength in angstrom
 \param fenergyWN energy in wavenumbers, cm^-1
 \return vacuum wavelength in angstrom
*/
inline realnum wn2angVac( double fenergyWN )
{
	return safe_div( 1e+8_r, realnum(fenergyWN) );
}

#endif /* LINES_SERVICE_H_ */
