/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PARSER_H_
#define PARSER_H_

 /**nWord determine whether match to a keyword occurs on command line,
   return value is 0 if no match, and position of match within string if hit 
	  \param *chKey
	  \param *chCard
 */ 

#include <cstdio>
#include <map>
#include "count_ptr.h"
#include "input.h"

const char * nWord(const char *chKey, 
	    const char *chCard);

class Parser;
class DepthTable;
class LineID;

typedef void (*OptionParser)(Parser &);

struct CloudyCommand {
	const char *name;
	OptionParser action;
};

bool isBoundaryChar(char c);

class Symbol {
public:
	enum tokens { ERROR, EOSTAT, NUMBER, NAME, STRING, OPERATOR, NTOKS };
	enum tokens toktype;
	string value;
   Symbol(enum tokens t, const string& v) : toktype(t), value(v) {}
};

/** Parser class holds pointer to string currently being analysed */
class Parser
{
	char m_card[INPUT_LINE_LENGTH];         // all-caps version of m_card_raw
	char m_card_raw[INPUT_LINE_LENGTH];     // raw input line with all comments stripped
	char m_card_comment[INPUT_LINE_LENGTH]; // input line with visible comments included, used for output
	long int m_len;
	long int m_off;
	bool m_lgEOL;
	const CloudyCommand * const m_Commands;
	std::map<string,double> m_symtab;
public:
	long int m_nqh;
	bool m_lgDSet, m_lgEOF;

	explicit Parser(void) :	m_Commands(NULL)
	{
		init();
	}
	explicit Parser(const CloudyCommand *commands) : m_Commands(commands) 
	{
		init();
	}
private:
	void init()
	{
		m_nqh = 0;
		m_lgDSet = m_lgEOF = false;
		setline("");
	}
	void newlineProcess();
	bool at_end() const;
	char current( void ) const
	{
		return m_card[m_off];
	}
	char current_raw() const;
	void skip_whitespace();
public:
	bool getline();
	void setline(const char * const card)
		{
			ASSERT(INPUT_LINE_LENGTH > 0);
			ASSERT(strlen(card) < (unsigned) INPUT_LINE_LENGTH);
			strncpy(m_card_comment,card,INPUT_LINE_LENGTH);
			string m_card_local( card );
			StripComment( m_card_local, true );
			// the ASSERTs above guarantee that this will not overflow buffer
			// as StripComment() can only make the string shorter
			strcpy(m_card_raw,m_card_local.c_str());
			newlineProcess();
		}

	void set_point(long int ipnt)
	{
		m_off = ipnt;
	}
	const char * nWord(const char *chKey) const;
	bool lgReachedEnd();
	void showLocation(FILE *io = ioQQQ) const;
public:
	long int GetElem( void ) const;
	double FFmtRead( void );
	double getNumberPlain( const char *chDesc );
	double getNumberCheck( const char *chDesc );
	double getNumberDefault( const char *chDesc, double fdef );
	double getNumberCheckLogLinNegImplLog( const char *chDesc );
	double getNumberCheckAlwaysLog( const char *chDesc );
	double getNumberCheckAlwaysLogLim( const char *chDesc, double flim );
	double getNumberDefaultAlwaysLog( const char *chDesc, double fdef );
	double getNumberDefaultNegImplLog( const char *chDesc, double fdef );
	bool lgEOL(void) const
	{
		return m_lgEOL;
	}
	void setEOL(bool val)
	{
		m_lgEOL = val;
	}
	NORETURN void NoNumb(const char *chDesc) const;
private:
	int nMatch1(const char *chKey) const
	{
		const char *p=chKey;

		while (isspace(*p))
			++p;

		for (const char *q=p; *q; ++q)
			ASSERT(!islower(*q));

		if ( !isBoundaryChar(*p))
		{
			const char *q = ::nWord(p, m_card);
			if (NULL == q)
				return 0;
			else
				return q-m_card+1;
		}
		else
		{
			// If the keyword starts with a member of the boundary character
			// set, can't require it to be preceded by one so revert to explicit
			// matching
			return ::nMatch(chKey, m_card);
		}
	}		
public:
	bool nMatch(const char *chKey) const
	{
		return nMatch1(chKey) != 0;
	}
	bool GetParam(const char *chKey, double *val)
	{
		int i = nMatch1(chKey);
		if (i > 0) {
			m_off = i-1;
			*val = FFmtRead();
		}
		return i>0;
	}
	bool GetRange(const char *chKey, double *val1, double *val2)
	{
		int i = nMatch1(chKey);
		if (i > 0) {
			m_off = i-1;
			*val1 = FFmtRead();
			*val2 = FFmtRead();
		}
		return i>0;
	}
	bool nMatchErase(const char *chKey)
	{
		const char *p=chKey;
		while (isspace(*p))
			++p;
		int i = nMatch1(p);
		bool found = (i != 0);
		if(found) 
		{
			char *ptr = m_card+i-1;
			const long len = strlen(p);
			/* erase this keyword, it upsets FFmtRead */
			for (long i=0; i<len; ++i)
			{
				ptr[i] = ' ';
			}
		}		
		return found;
	}
	bool hasCommand(const char *s2);
	bool Command(const char *name, OptionParser doOpts)
	{
		bool lgFound = hasCommand(name);
		if ( lgFound )
			(*doOpts)(*this);
		return lgFound;
	}
	bool isComment(void) const;
	bool isVar(void) const;
	std::string getVarName(void);
	void doSetVar(void);
	void echo(void) const;
	bool last(void) const;
	int PrintLine(FILE *fp) const
	{
		return fprintf( fp, " ==%-.80s==\n", m_card_raw);
	}
	NORETURN void CommandError( void ) const;
	NORETURN void Error( const char *msg ) const;
	NORETURN void StringError( ) const;
	int GetQuote( string& chLabel );
	const char *StandardEnergyUnit(void) const;
	string StandardFluxUnit(void) const;
	string getFirstChunk(long i);
	string getRawTail()
		{
			return string(m_card_raw+m_off);
		}
	void help(FILE *fp) const;
	double getWave();
	double getWaveOpt();
	LineID getLineID();
	Symbol getSymbol();
	void getPairs(vector<double>& a, vector<double> & b);
	void readList(vector<string>& list, const char *chName);
	void readLaw(DepthTable& table);
};

/** Links text string to an action on a specified argument */
template <typename V>
class KeyAction {
	const char * const m_keyword;
	V m_action;
public:
	KeyAction(const char *keyword, const V &action) :
		m_keyword(keyword), m_action(action) {}
		
	const char *key(void) const
	{
		return m_keyword;
	}
	void operator()(realnum *v) const
	{
		m_action(v);
	}
};

/** Helper template to make it easier to generate KeyActions */
template <typename V>
inline KeyAction<V> MakeKeyAction(const char *keyword, const V &action)
{
	return KeyAction<V>(keyword, action);
}

/** Generator for functors which convert the unit of their argument */
class UnitConverter
{
	const realnum m_unit;
public:
	UnitConverter ( double unit ) : m_unit((realnum)unit) {}
		
	void operator()( realnum *t ) const
	{
		*t *= m_unit;
	}
};

/** Interate through a list of KeyActions: apply the first which
	 matches and then quit */
template <typename T, typename V>
bool parserProcess(Parser &p, T *list, unsigned long nlist, V *value)
{
	bool lgFound = false;
	for (unsigned long option=0; option < nlist; ++option)
	{
		if( p.nWord( list[option].key() ) )
		{
			list[option]( value );
			lgFound = true;
			break;
		}
	}
	return lgFound;
}

/**ParseCosmicRays parse the cosmic rays command 
\param *chCard
*/
void ParseCosmicRays( Parser &p );

/**ParseCosmology parse the cosmology command 
\param *chCard
*/
void ParseCosmology( Parser &p );

/**ParseAbundances parse and read in composition as set by abundances command 
\param *chCard
\param lgDSet
*/
void ParseAbundancesNonSolar(Parser &p);

void ParseAbundances(Parser &p);

/**ParseDont parse the dont command */
void ParseDont(Parser &p);

/**ParseSave parse the save command 
\param *chCard
*/
void ParseSave(Parser &p);

void parse_save_line(Parser &p, 
		     /* true, return rel intensity, false, log of luminosity or intensity I */
		     bool lgLog3,
			  ostringstream& chHeader,
							long int ipPun
);

void parse_save_average(Parser &p,
			/* the file we will write to */
			long int ipPun, 
			ostringstream& chHeader);

void parse_save_colden(Parser &p,
		       /* the header for the file, a list of identifications */
		       ostringstream& chHeader);

void Parse_Save_Line_RT(Parser &p);

/**ParseAge - parse the age command */
void ParseAge(Parser &p);

/**ParseAgn parse parameters for the AGN continuum shape command 
\param *chCard
*/
void ParseAgn(Parser &p);

/**ParseState save or recover previous state of the code 
\param *chCard
*/
void ParseState(Parser &p);

/** parse the blackbody command 
\param *chCard input command line, already changed to caps
\param *nqh counter for which continuum source this is
\param *ar1 optional area that might be set here
*/
void ParseBlackbody(Parser &p);					

/**ParseCompile compile werner or kurucz model atmospheres into cloudy format, by K Volk 
\param *chCard
*/
void ParseCompile(Parser &p );

/**ParseConstant parse the constant ... command */
void ParseConstant(Parser &p);

/**ParseDLaw parse parameters on the dlaw command so set some density vs depth law
\param *chCard
*/
void ParseDLaw(Parser &p );

/**ParseTLaw parse parameters on the tlaw command to set some temperature vs depth 
\param *chCard
*/
void ParseTLaw(Parser &p);

/**ParseDrive parse the drive command - drive calls to various subs 
\param *chCard
*/
void ParseDrive(Parser &p );

/**ParseGrain parse parameters on Peter's version of the grains command 
\param *chCard
\param *lgDset
*/
void ParseGrain(Parser &p);

/**ParseFluc parse the fluctuations command */
void ParseFluc(Parser &p);

/**ParseHDEN parse the HDEN command */
void ParseHDEN(Parser &p);

/**ParseDatabaseISO parse the atom XX-like command, to set options for iso sequences 
\param ipISO
\param *chCard
*/
void ParseDatabaseISO(long ipISO, Parser &p);

/**ParseDatabaseH2 parse information from the rotor command line 
\param *chCard
*/
void ParseDatabaseH2(Parser &p );

/**ParseGrid parse the grid command 
\param *chCard
*/
void ParseGrid(Parser &p);

/**ParseInit parse the init command */
void ParseInit(Parser &p);

/**ParseInterp parse parameters on interpolate command 
\param *chCard
\param *lgEOF
*/
void ParseInterp(Parser &p);

/**ParseIonParI parse the ionization parameter command (IONI variant)
\param *nqh
\param *chCard
\param *chType
*/
void ParseIonParI(Parser &p);

/**ParseIonParX parse the ionization parameter command (XI variant)
\param *nqh
\param *chCard
\param *chType
*/

void ParseIonParX(Parser &p);
/**ParseIonPar parse the ionization parameter command 
\param *nqh
\param *chCard
\param *chType
*/
void ParseIonPar(Parser &p,
					  char chType);

/**ParseNorm parse parameters on the normalize command 
\param *chCard
*/
void ParseNorm(Parser &p);

/**ParseOptimize parse the optimize command 
\param *chCard
*/
void ParseOptimize(Parser &p);

/**ParsePrint parse the print command  
\param *chCard
*/
void ParsePrint(Parser &p );

/**ParseRadius parse the radius command */
void ParseRadius(Parser &p);

/**ParseSet parse the set command */
void ParseSet(Parser &p);

/**ParseTable parse the table read command 
\param *nqh
\param *chCard
\param *ar1
*/
void ParseTable(Parser &p);

/**ParseTrace parse the trace command */
void ParseTrace(Parser &p);

/*ParseExtinguish parse the extinguish command */
void ParseExtinguish( Parser &p );

/*ParseIlluminate parse the illuminate command */
void ParseIlluminate( Parser &p );

/*ParseCaseB - parse the Case B command */
void ParseCaseB(Parser &p );

/**ParseTest parse the test command */
void ParseTest(Parser &p);

/**ParseAbsMag parse the absolute magnitude command */
void ParseAbsMag(Parser &p);

/**ParseBackgrd parse the background continuum command */
void ParseBackgrd(Parser &p);

/**ParseCoronal parse the cronal equilibrum command */
void ParseCoronal(Parser &p);

/**ParseElement parse options on element command */
void ParseElement(Parser &p);

/**ParseCMB parse parameters from fireball command 
\param z
\param *nqh
\param *ar1
*/
void ParseCMB(double z, 
  long int *nqh);

/**ParseF_nu parse intensity command parameters 
\param *chCard
\param *nqh
\param *ar1
\param *chType
\param lgNU2
*/
void ParseF_nu(
  Parser &p, 
  const char *chType, 
  bool lgNU2);

/**ParseGlobule parse parameters off the globule command 
\param *chCard
*/
void ParseGlobule(Parser &p);

/**ParseRangeOption parse the range option on the luminosity command */
void ParseRangeOption(Parser &p);

/**ParseMap parse map command to produce map of heating and cooling */
void ParseMap(Parser &p);

/**ParseMetal parse parameters on metal command */
void ParseMetal(Parser &p);

void ParsePrtLineSum(Parser &p);

/**ParsePlot parse the plot command */
void ParsePlot(Parser &p);

/**ParsePowerlawContinuum parse the power law continuum command */
void ParsePowerlawContinuum(Parser &p);

/**ParseRatio parse the ratio command */
void ParseRatio(Parser &p);

/**ParseSphere parse the sphere command */
void ParseSphere(Parser &p);

/**ParseStop parse the stop command */
void ParseStop(Parser &p);

/**ParseCrashDo any of several tests to check that the code can crash 
\param *chCard
*/
void ParseCrashDo(Parser &p);
 
class Option {
public:
	enum Opts { BOOL, LONG, REAL, OPTION, STRING, INVALID } opttype;
	enum Quoted { NOTQUOTED, QUOTED};
	union {
		bool   l;
		long   i;
		double r;
	};
	string s;
	Option(bool val)
	{
		opttype=BOOL;
		l = val;
	}
	Option(long val)
	{
		opttype=LONG;
		i = val;
	}
	Option(double val)
	{
		opttype=REAL;
		r = val;
	}
	Option(const string& val, enum Quoted q)
	{
		if (q == QUOTED)
			opttype=STRING;
		else
			opttype=OPTION;
		s = val;
	}
};

class Properties
{
	bool m_lgDone;
public:
	vector< pair<string, count_ptr<Option> > > p;
   Properties() : m_lgDone(false) {}
	void setDone() 
	{
		m_lgDone = true;
	}
	bool isDone() const
	{
		return m_lgDone;
	}
};

#endif // _PARSER_H_
