/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "parser.h"
#include "called.h"
#include "flux.h"
#include "input.h"
#include "elementnames.h"
#include "service.h"
#include "depth_table.h"
#include "lines.h"

#include <deque>
namespace
{
	class Token
	{
	public:
		enum symType { symNull, symNumber, symOp, symVar };
		string s;
		symType t;
		explicit Token(enum symType type) : s(""), t(type) {}
		explicit Token() : s(""), t(symNull) {}
	};
}

typedef std::map<string,double> symtab;
STATIC bool ParseExpr(deque<Token> &chTokens, vector<double> &valstack,
	const symtab &tab);

const char *Parser::nWord(const char *chKey) const
{
	return ::nWord(chKey, m_card);
}

bool Parser::at_end() const
{
	return current() == '\0';
}
char Parser::current_raw() const
{
	return m_card_raw[m_off];
}
void Parser::skip_whitespace()
{
	while (!at_end() && isspace(current()))
		++m_off;
}
void Parser::newlineProcess(void)
{
	strncpy(m_card,m_card_raw,INPUT_LINE_LENGTH);
	// transform(m_card.begin(),m_card.end(),m_card.begin(),toupper);
	::caps(m_card);
	m_len = strlen(m_card);
	m_off = 0;
	m_lgEOL = false;
}

/*nWord determine whether match to a keyword occurs on command line,
 * return value is 0 if no match, and position of match within string if hit */
const char *nWord(const char *chKey, 
	    const char *chCard)
{
	DEBUG_ENTRY( "nWord()" );

	// Ignore leading space in chKey -- logic below is designed
	// to avoid the need to include this in the first place
	while (isspace(*chKey))
	{
		++chKey;
	}

	const long lenkey = strlen(chKey);
	ASSERT( lenkey > 0 );

	bool atBoundary = true, inQuote=false;
	for (const char *ptr = chCard; *ptr; ++ptr)
	{
		if (!inQuote)
		{
			if (*ptr == '\"')
			{
				inQuote = true;
			}
			else
			{		
				if ( atBoundary && strncmp( ptr, chKey, lenkey) == 0 )
				{
					return ptr;
				}
				
				atBoundary = isBoundaryChar(*ptr);
			}
		}
		else
		{
			if (*ptr == '\"')
			{
				inQuote = false;
			}		
		}
	}

	return NULL;
}

/* check if there is something left at the end of the line after skipping whitespace */
bool Parser::lgReachedEnd()
{
	skip_whitespace();
	return at_end();
}

void Parser::showLocation(FILE *io) const
{
	fprintf( io, "  %s\n", m_card_comment );
	fprintf( io, "  " );
	// make sure tabs are treated correctly....
	for( long i=0; i < m_off; ++i )
	{
		if( m_card_comment[i] == '\t' )
			fputc( '\t', io );
		else
			fputc( ' ', io );
	}
	fprintf( io, "^\n" );
}

bool isBoundaryChar(char c)
{
	const bool lgAnyWhitespacePrecedesWord = false;
	
	if (lgAnyWhitespacePrecedesWord)
		return isspace(c) ? true : false ;
	else 	// Words are strings starting with A-Z, a-z or _
		return (! isalpha(c) ) && c != '_';
}

bool Parser::isComment(void) const
{
	return lgInputComment(m_card_comment);
}
bool Parser::isVar(void) const
{
	return ( current()=='$' );
}
std::string Parser::getVarName(void)
{
	std::string name("");
	while (!at_end())
	{
		char c = current();
		if (!(isalnum(c) || c == '_'))
			break;
		name += c;
		++m_off;
	}
	return name;
}
void Parser::doSetVar(void)
{
	DEBUG_ENTRY( "Parser::doSetVar()" );
	char c='\0';
	++m_off;
	std::string name = getVarName();
	while (!at_end())
	{
		c = current();
		++m_off;
		if (c == '=')
			break;
	}
	if (at_end())
	{
		fprintf(ioQQQ,"Expected '=' in variable definition\n");
		cdEXIT(EXIT_FAILURE);
	}
	while (!at_end())
	{
		c = current();
		if (c != ' ')
			break;
		++m_off;
	}
	m_symtab[name] = FFmtRead();
}

void Parser::echo(void) const
{
	/* >>chng 04 jan 21, add HIDE option, mostly for print quiet command */
	if( called.lgTalk && !::nMatch("HIDE",m_card) )
		fprintf( ioQQQ, "%23c* %-80s*\n", ' ', m_card_comment );
}

bool Parser::last(void) const
{
	// using m_card_raw here could lead to spurious EOF detection if the line
	// is a full-line comment, in which case m_card_raw would be an empty line
	return m_lgEOF || lgInputEOF(m_card_comment);
}

NORETURN void Parser::StringError( ) const
{
	Error(
		" A filename or label must be specified within double quotes, but no quotes were encountered on this command.\n"
		" Name must be surrounded by exactly two double quotes, like \"name.txt\". \n" 
		);
}

/** find string between pair of double quotes.  Returns 0 for success, 1 for failure
 * typical use is to follow failure by p.StringError(); */
int Parser::GetQuote( string& chLabel )
{
	DEBUG_ENTRY( "Parser::GetQuote()" );

	string str_raw( m_card_raw );

	/* find first quote start of string, string begins and ends with quotes */
	size_t p0 = str_raw.find( '\"' );
	/* get pointer to next quote and read label */
	size_t p1 = ( p0 != string::npos ) ? GetString( str_raw, p0, chLabel ) : string::npos;

	/* check that pointers are valid */
	if( p0 == string::npos || p1 == string::npos )
	{		
		/* this branch, ok if not present, return null string in that case */
		chLabel = "";
		/* return value of 1 indicates did not find double quotes */
		return 1;
	}

	/* get pointer to first quote in local copy of line image in calling routine */
	char *i0 = m_card_raw+p0;
	char *i1 = m_card_raw+p1;
	char *iLoc = m_card+p0;

	// blank out label once finished, to not be picked up later
	// erase quotes as well, so that we can find second label, by PvH
	while( i0 < i1 )
	{
		*i0++ = ' ';
		*iLoc++ = ' ';
	}
	/* return condition of 0 indicates success */
	return 0;
}

NORETURN void Parser::Error(const char* msg) const
{
	DEBUG_ENTRY( "Parser::Error()" );
	fprintf( ioQQQ, "  Parser failure\n");
	if (msg)
		fprintf(ioQQQ,"%s",msg);
	fprintf( ioQQQ, " The line image was\n");
	PrintLine(ioQQQ);
	fprintf( ioQQQ, " Sorry.\n" );
	cdEXIT(EXIT_FAILURE);
}
NORETURN void Parser::CommandError(void) const
{
	DEBUG_ENTRY( "Parser::CommandError()" );
	fprintf( ioQQQ, "  Unrecognized command. Key=\"%4.4s\".  This is routine ParseCommands.\n", 
				m_card );
	fprintf( ioQQQ, " The line image was\n");
	PrintLine(ioQQQ);
	fprintf( ioQQQ, " Sorry.\n" );
	if( lgIsExpungedCommentSeq(m_card) )
		fprintf( ioQQQ, " This looks like an invalid old-style comment."
			 " Please use scripts/ccc.pl to convert this input script to convert to #.\n" );
	cdEXIT(EXIT_FAILURE);
}
bool Parser::getline(void)
{
	input.readarray(m_card_raw,m_card_comment,&m_lgEOF);
	newlineProcess();
	if (m_lgEOF)
		return false;
	else
	{
		input.lgVisible[input.nRead] = this->nMatch("HIDE") ? false : input.lgVisibilityStatus;
		return true;
	}
}
	
const char *Parser::StandardEnergyUnit(void) const
{
	return ::StandardEnergyUnit(m_card);
}
string Parser::StandardFluxUnit(void) const
{
	return ::StandardFluxUnit(m_card);
}
void Parser::help(FILE *fp) const
{
	DEBUG_ENTRY( "Parser::help()" );
	fprintf(fp,"Available commands are:\n\n");
	long int i=0, l=0, len;
	while (1)
	{
		len = strlen(m_Commands[i].name);
		if (l+len+2 > 80)
		{
			fprintf(fp,"\n");
			l = 0;
		}
		l += len+2;
		fprintf(fp,"%s",m_Commands[i].name);
		++i;
		if (m_Commands[i].name == NULL)
			break;
		fprintf(fp,", ");
	}
				  
	fprintf(fp,"\n\nSorry, no further help available yet -- try Hazy.\n\n");
	cdEXIT(EXIT_SUCCESS);
}

/*GetElem scans line image, finds element. returns atomic number j, 
 * on C scale, -1 if no hit.  chCARD_CAPS must be in CAPS to hit element */
long int Parser::GetElem(void ) const
{
	int i;

	DEBUG_ENTRY( "Parser::GetElem()" );

	/* find which element */

	/* >>>chng 99 apr 17, lower limit to loop had been 1, so search started with helium,
	 * change to 0 so we can pick up hydrogen.  needed for parseasserts command */
	/* find match with element name, start with helium */
	for( i=0; i<(int)LIMELM; ++i )
	{
		if( nMatch( elementnames.chElementNameShort[i] ) )
		{
			/* return value is in C counting, hydrogen would be 0*/
			return i;
		}
	}
	/* fall through, did not hit, return -1 as error condition */
	return (-1 );
}

/*NoNumb general error handler for no numbers on input line */
NORETURN void Parser::NoNumb(const char * chDesc) const
{
	DEBUG_ENTRY( "Parser::NoNumb()" );

	/* general catch-all for no number when there should have been */
	fprintf( ioQQQ, " There is a problem on the following command line:\n" );
	fprintf( ioQQQ, " %s\n", m_card_raw );
	fprintf( ioQQQ, " A value for %s should have been on this line.\n   Sorry.\n",chDesc );
	cdEXIT(EXIT_FAILURE);
 }

double Parser::getWaveOpt()
{
	double val = FFmtRead();
	/* check for optional micron or cm units, else interpret as Angstroms */
	if( current() == 'A' )
	{
		/* Angstrom */
		++m_off;
	}
	else if( current() == 'M' )
	{
		/* microns */
		val *= 1e4;
		++m_off;
	}
	else if( current() == 'C' )
	{
		/* centimeters */
		val *= 1e8;
		++m_off;
	}
	return val;
}
double Parser::getWave()
{
	double val = getWaveOpt();
	if( lgEOL() )
	{
		NoNumb("wavelength");
	}
	return val;
}
double Parser::getNumberPlain( const char * )
{
	return FFmtRead();
}
double Parser::getNumberCheck( const char *chDesc )
{
	double val = FFmtRead();
	if( lgEOL() )
	{
		NoNumb(chDesc);
	}
	return val;
}
double Parser::getNumberDefault( const char *, double fdef )
{
	double val = FFmtRead();
	if( lgEOL() )
	{
		val = fdef;
	}
	return val;
}
double Parser::getNumberCheckLogLinNegImplLog( const char *chDesc )
{
	double val = getNumberCheck(chDesc);
	if( nMatch(" LOG") )
	{
		val = exp10(val);
	}		
	else if(! nMatch("LINE") )
	{
		/* log, linear not specified, neg so log */
		if( val <= 0. )
		{
			val = exp10(val);
		}
	}
	return val;
}
double Parser::getNumberCheckAlwaysLog( const char *chDesc )
{
	double val = getNumberCheck(chDesc);
	val = exp10( val);
	return val;
}
double Parser::getNumberCheckAlwaysLogLim( const char *chDesc, double flim )
{
	double val = getNumberCheck(chDesc);
	if ( val > flim )
	{
		fprintf(ioQQQ,"WARNING - the log of %s is too "
				  "large, I shall probably crash.  The value was %.2e\n",
				  chDesc, val );
		fflush(ioQQQ);		
	}
	val = exp10( val);
	return val;
}
double Parser::getNumberDefaultAlwaysLog( const char *, double fdef )
{
	double val = exp10(FFmtRead());
	if ( lgEOL() )
	{
		val = fdef;
	}
	return val;
}
double Parser::getNumberDefaultNegImplLog( const char *, double fdef )
{
	double val = FFmtRead();
	if ( lgEOL() )
	{
		val = fdef;
	}
	if (val < 0.0)
	{
		val = exp10(val);
	}
	return val;
}

/*FFmtRead scan input line for free format number */


double Parser::FFmtRead(void)
{

	DEBUG_ENTRY( "Parser::FFmtRead()" );
	
	// Look for start of next expression
	while( m_off < m_len) 
	{
		if ( current() == '$' )
			break;
		int loff = m_off;
		char lchr = current();
		if( lchr == '-' || lchr == '+' )
		{
			++loff;
			lchr = m_card[loff];
		}
		if( lchr == '.' )
		{
			++loff;
			lchr = m_card[loff];
		}
		if( isdigit(lchr) )
			break;
		++m_off;
	}

	if( m_off == m_len )
	{
		m_lgEOL = true;
		return 0.;
	}

	// Lexer for expression
	deque<Token> chTokens(0);
	for(char chr = current();
		 m_off < m_len &&
			 (isdigit(chr) || chr == '.' || chr == '-' || chr == '+' 
			  || chr == 'e' || chr == 'E' || chr == '^' || chr == '*' || chr == '/' 
			  || chr == '$')
			 ;	 chr = current())
	{
		if (chr == '^' || chr == '*' || chr == '/' ) 
		{
			chTokens.push_back(Token(Token::symOp));
			chTokens.back().s += chr;
			++m_off;
		}
		else if (chr == '$')
		{
			chTokens.push_back(Token(Token::symVar));
			++m_off;
			chTokens.back().s += getVarName();
		}
		else
		{
			if (chTokens.size() == 0 || chTokens.back().t != Token::symNumber)
				chTokens.push_back(Token(Token::symNumber));
			chTokens.back().s += chr;
			++m_off;
		}
	}

	ASSERT (chTokens.size() != 0);

	// Parse tokens
	vector<double> valstack;
	const bool lgParseOK = ParseExpr(chTokens, valstack, m_symtab);
	if (!lgParseOK || 1 != valstack.size())
	{
		fprintf(ioQQQ," PROBLEM - syntax error in number or expression on line\n");
		fprintf(ioQQQ, "== %-80s ==\n",m_card);
		m_lgEOL = true;
		return 0.;
	}

	double value = valstack[0];

	m_lgEOL = false;
	return value;
}

string Parser::getFirstChunk(long nchar)
{
	DEBUG_ENTRY( "Parser::getFirstChunk()" );
	if (m_len < nchar)
	{
		fprintf(ioQQQ,
				  "PROBLEM --"
				  " input line too short to provide %ld character label\n"
				  "== %-80s ==\n", nchar,m_card);
		cdEXIT(EXIT_FAILURE);
	}
	m_off = nchar;
	return string(m_card).substr(0,nchar);
}

LineID Parser::getLineID()
{
	DEBUG_ENTRY( "Parser::getLineID()" );
	LineID line;
	if( m_card_raw[0] == '\"' )
	{
		if( GetQuote( line.chLabel ) == 0 )
			trimTrailingWhiteSpace( line.chLabel );
		else
		{
			fprintf( ioQQQ, "getLineID found invalid quoted string:\n" );
			showLocation();
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{	 
		/* order on line is label (col 1-4), wavelength */
		line.chLabel = getFirstChunk(4);
	}

	// Normalize common error "H 1 " or "H 1" for "H  1"
	if ( line.chLabel.size() == 3 || line.chLabel.size() == 4 )
	{
		if ( line.chLabel[1] == ' ' &&
			  ( line.chLabel.size() == 3 || line.chLabel[3] == ' ') )
		{
			fprintf(ioQQQ,"WARNING: read \"%s\" as spectrum\n",line.chLabel.c_str());
			if (line.chLabel.size() == 3)
				line.chLabel += line.chLabel[2];
			else
				line.chLabel[3] = line.chLabel[2];
			line.chLabel[2] = ' ';
			fprintf(ioQQQ,"Assuming required spectrum is \"%s\"\n",line.chLabel.c_str());
		}
	}

	/* now get wavelength */
	line.wave = (realnum)getWave();
	return line;
}

// Simple recursive descent parser for expressions
//
// for discussion, see e.g. http://www.ddj.com/architect/184406384
//
// for a possibly more efficient alternative, see
// http://eli.thegreenplace.net/2010/01/02/top-down-operator-precedence-parsing/

STATIC bool ParseNumber(deque<Token> &chTokens, vector<double> &valstack,
	const symtab &tab)
{
	DEBUG_ENTRY( "ParseNumber()" );
	if ( chTokens.size() < 1)
		return false;

	if (Token::symNumber == chTokens[0].t)
	{
		valstack.push_back(atof(chTokens[0].s.c_str()));
		chTokens.pop_front();
		return true;
	}
	if (Token::symVar == chTokens[0].t)
	{
		symtab::const_iterator var = tab.find(chTokens[0].s);
		if (var == tab.end())
		{
			fprintf(ioQQQ,"ERROR: No value found for variable $%s\n",
					  chTokens[0].s.c_str());
			cdEXIT(EXIT_FAILURE);
		}
		valstack.push_back(var->second);
		chTokens.pop_front();
		return true;
	}

	return false;
}

STATIC bool doop(vector<double> &valstack, const string &op)
{
	DEBUG_ENTRY( "doop()" );
	const double v2 = valstack.back();
	valstack.pop_back();
	const double v1 = valstack.back();
	valstack.pop_back();
	double result;
	if (op == "^")
	{
		result = pow(v1,v2);
	}
	else if (op == "*")
	{
		result = v1*v2;
	}
	else if (op == "/")
	{
		result = v1/v2;
	}
	else
	{
		fprintf(ioQQQ,"Unknown operator '%s'\n",op.c_str());
		return false;
	}
	valstack.push_back(result);
	return true;
}

STATIC bool ParseExp(deque<Token> &chTokens, vector<double> &valstack,
	const symtab& tab)
{
	DEBUG_ENTRY( "ParseExp()" );
	// Right-associative -- need to buffer into stack
	vector<string> opstack;
	if (!ParseNumber(chTokens, valstack, tab))
		return false;

	while (1)
	{
		if ( chTokens.size() == 0 )
			break;
		
		if ( chTokens.size() < 2 )
			return false;
		
		if ( Token::symOp != chTokens[0].t || "^" != chTokens[0].s )
			break;
		
		opstack.push_back(chTokens[0].s);
		chTokens.pop_front();
	
		if (!ParseNumber(chTokens, valstack, tab))
			return false;
	}

	while (!opstack.empty())
	{	
		if (!doop(valstack, opstack.back()))
			return false;
		opstack.pop_back();
	}
	return true;
}

STATIC bool ParseProduct(deque<Token> &chTokens, vector<double> &valstack,
	const symtab& tab)
{
	DEBUG_ENTRY( "ParseProduct()" );
	// Left-associative
	if (!ParseExp(chTokens, valstack, tab))
		return false;

	while ( chTokens.size() > 0 && 
			  Token::symOp == chTokens[0].t &&
			  ( "*" == chTokens[0].s || "/" ==  chTokens[0].s ) ) 
	{
		string op = chTokens[0].s;
		chTokens.pop_front();
		
		if (!ParseExp(chTokens, valstack, tab))
			return false;
		
		if (!doop(valstack, op))
			return false;
	}
	return true;
}

STATIC bool ParseExpr(deque<Token> &chTokens, vector<double> &valstack,
	const symtab& tab)
{
	DEBUG_ENTRY( "ParseExpr()" );
	if (ParseProduct(chTokens, valstack,tab))
		return true;
	return false;
}

bool Parser::hasCommand(const char *s2)
{
	DEBUG_ENTRY( "Parser::hasCommand()" );
	size_t len = strlen(s2);
	if (::strncmp(m_card, s2, len) != 0)
		return false;

	while (isspace(m_card[len-1]))
	{
		--len;
	}

	// Fast forward over any tail to command, i.e. ignore non-space
	// characters after len but before next space
	m_off = len;
	while (!at_end() && (isalnum(current()) || current() == '_'))
		++m_off;
	skip_whitespace();

	return true;
}

void Parser::getPairs(vector<double>& a, vector<double> & b)
{
	DEBUG_ENTRY( "Parser::getPairs()" );
	a.resize(0);
	b.resize(0);
	for(;;)
	{
		getline();
		
		if( m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading element list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		
		if( isComment())
			continue;

		/* this would be a line starting with END to say end on list */
		if( hasCommand("END") )
		{
			return;
		}
		a.push_back(FFmtRead());
		if( lgEOL() )
		{
			fprintf( ioQQQ, " There must be a two numbers on this line, or END.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		b.push_back(FFmtRead());
		if( lgEOL() )
		{
			fprintf( ioQQQ, " There must be a two numbers on this line, or END.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
}


inline Symbol maybeNumber(bool numOK, const Symbol &s)
{
	DEBUG_ENTRY( "maybeNumber()" );
	if (numOK)
	{
		return s;
	}
	else
	{
		fprintf(ioQQQ,"Parser error, incomplete number"
				  " at end of input term, error after '%s'\n",s.value.c_str());
		return Symbol(Symbol::ERROR,s.value);
	}
}

Symbol Parser::getSymbol()
{
	DEBUG_ENTRY( "Parser::getSymbol()" );
	// Eat leading space
	while (!at_end())
	{
		char c = current();
		if (c != ' ' && c != '\t')
			break;
		++m_off;
	}

	if ( at_end() || current() == '\n')
		return Symbol(Symbol::EOSTAT,"");

	if (isdigit(current()) || current() == '-' || current() == '.')
	{
		Symbol s(Symbol::NUMBER,"");
		bool numOK=false;
		if (current() == '-')
		{
			s.value += current();
			++m_off;
			if (at_end())
				return maybeNumber(numOK,s);
		}
		if (isdigit(current()))
		{
			numOK = true;
			do
			{
				s.value += current();
				++m_off;
			}
			while (!at_end() && isdigit(current()));
			if (at_end() )
				return maybeNumber(numOK,s);
		}
		if (current() == '.')
		{
			s.value += current();
			++m_off;
		}
		while (!at_end() && isdigit(current()))
		{
			numOK = true;
			s.value += current();
			++m_off;
		}
		if ( at_end() || current() != 'E' || !numOK)
		{
			return maybeNumber(numOK,s);
		}
		s.value += current();
		++m_off;
		numOK = false;
		if (current() == '-' || current() == '+')
		{
			s.value += current();
			++m_off;
			if (at_end())
				return maybeNumber(numOK,s);
		}
		while (!at_end() && isdigit(current()))
		{
			numOK = true;
			s.value += current();
			++m_off;
		}
		return maybeNumber(numOK,s);
	}

	if (isalpha(current()))
	{
		Symbol s(Symbol::NAME,"");
		do
		{
			s.value += current();
			++m_off;
		}
		while (!at_end() && (isalnum(current()) || current() == '_'));
		return s;
	}

	if ( current() == '"' )
	{
		Symbol s(Symbol::STRING,"");
		++m_off;
		while (!at_end() && current() != '\"')
		{
			if (current() == '\\')
			{
				++m_off;
				if (at_end())
				{
					fprintf(ioQQQ,"Parser error, escape character '\\'"
							  " at end of input term\n");
					return Symbol(Symbol::ERROR,s.value);
				}
			}
			s.value += current_raw();
			++m_off;
		}

		if (at_end())
		{
			fprintf(ioQQQ,"Parser error, unterminated string\n");
			return Symbol(Symbol::ERROR,s.value);
		}
		++m_off;
		return s;
	}

	if ( current() == ',' || current() == '(' || current() == ')' || current() == '=' )
	{
		Symbol s(Symbol::OPERATOR,"");
		s.value += current();
		++m_off;
		return s;
	}

	fprintf(ioQQQ,"Parser error, character not recognized '%c'\n",
		current());
	return Symbol(Symbol::ERROR,"");
}
void Parser::readList(vector<string>& list, const char* chName)
{
	DEBUG_ENTRY( "Parser::readList()" );
	list.clear();
	while( 1 )
	{
		getline();
		if( m_lgEOF )
		{
			fprintf( ioQQQ, 
				" Save %s hit EOF while reading list; use END to end list.\n" ,chName);
			fprintf( ioQQQ,
				" This command requires either a species within quotes or the keyword ALL.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		if( isComment())
			continue;
		if (hasCommand("END" ))
			break;
		string chTerm;
		if (!GetQuote(chTerm))
			list.push_back(chTerm);
		else
			list.push_back(m_card_raw);
	}	
}

void Parser::readLaw(DepthTable& table)
{
	DEBUG_ENTRY( "Parser::readLaw()" );
	if( nMatch("DEPT") )
	{
		table.lgDepth = true;
	}
	else
	{
		table.lgDepth = false;
	}
	if (table.nvals != 0)
	{
		fprintf( ioQQQ, " Warning: over-writing existing table\n" );
		table.clear();
	}

	getline();
	table.dist.push_back(FFmtRead());
	table.val.push_back(FFmtRead());
	if( lgEOL() )
	{
		fprintf( ioQQQ, " No pairs entered - can\'t interpolate.\n Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	table.nvals = 2;
	bool lgEnd = false;
	
	/* read pairs of numbers until we find line starting with END */
	/* >>chng 04 jan 27, loop to LIMTABDLAW from LIMTABD, as per
	 * var definitions, caught by Will Henney */
	while( !lgEnd )
	{
		getline();
		lgEnd = m_lgEOF;
		if( !lgEnd )
		{
			lgEnd = hasCommand("END");
		}
		
		if( !lgEnd )
		{
			double dist = FFmtRead();
			double val = FFmtRead();
			if (lgEOL())
				NoNumb("radius, value pair on each line");
			table.dist.push_back( dist );
			table.val.push_back( val );
			table.nvals += 1;
		}
	}
	--table.nvals;

	for( long i=1; i < table.nvals; i++ )
	{
		/* the radius values are assumed to be strictly increasing */
		if( table.dist[i] <= table.dist[i-1] )
		{
			fprintf( ioQQQ, " Radii must be in increasing order.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
}
