/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "MObjBuilder.h"
#include "MMath.h"
#include "MEvaluate.h"
#include "MFunctions.h"
#include "assert.h"
#include <map>
#include <string>
#include <math.h>
using namespace std;

map<string, FUNCPTR>		FNC;	// 1-D functions
map<string, FUNC2PTR>		FNC2;	// 2-D functions
map<string, FUNCNPTR>		FNCN;	// N-D functions
map<string, FUNCMATPTR>		FMAT;	// Matrix functions
map<string, double>			CNT;	// predefined constants

void init_function_lists()
{
	// 1-parameter functions
	FNC["cos"  ] = cos;
	FNC["sin"  ] = sin;
	FNC["tan"  ] = tan;
	FNC["csc"  ] = csc;
	FNC["sec"  ] = sec;
	FNC["cot"  ] = cot;
	FNC["acos" ] = acos;
	FNC["asin" ] = asin;
	FNC["atan" ] = atan;
	FNC["cosh" ] = cosh;
	FNC["sinh" ] = sinh;
	FNC["tanh" ] = tanh;
	FNC["acosh"] = acosh;
	FNC["ln"   ] = log;
	FNC["log"  ] = log10;
	FNC["exp"  ] = exp;
	FNC["sqrt" ] = sqrt;
	FNC["abs"  ] = fabs;
	FNC["sgn"  ] = sgn;
	FNC["H"    ] = heaviside;
	FNC["step" ] = unit_step;

#ifdef WIN32
	FNC["J0"   ] = _j0;
	FNC["J1"   ] = _j1;
	FNC["Y0"   ] = _y0;
	FNC["Y1"   ] = _y1;
#endif
	FNC["sinc" ] = sinc;
	FNC["fl"   ] = fl;
	FNC["fac"  ] = fac;
	FNC["erf"  ] = erf;
	FNC["erfc" ] = erfc;
	FNC["tgamma"] = tgamma;

	// 2-parameter functions
#ifdef WIN32
	FNC2["Jn" ]   = jn;
	FNC2["Yn" ]   = yn;
#endif
	FNC2["atan2"] = atan2;
	FNC2["pow"  ] = pow;
	FNC2["mod"  ] = fmod;
	FNC2["Tn"]    = chebyshev;
	FNC2["binomial"] = binomial;

	// n-parameter functions
	FNCN["max"] = fmax;
	FNCN["min"] = fmin;
	FNCN["avg"] = avg;

	// Matrix functions
	FMAT["transpose"] = matrix_transpose;
	FMAT["trace"    ] = matrix_trace;
	FMAT["det"      ] = matrix_determinant;
	FMAT["inverse"  ] = matrix_inverse;

	// predefined constants
	CNT["pi"   ] = 3.1415926535897932385;	// "pi"
	CNT["e"    ] = 2.7182818284590452354;	// exp(1)
	CNT["gamma"] = 0.57721566490153286061;	// Euler-Mascheroni constant
	CNT["_inf_"] = 1e308;	// "infinity"
}

//-----------------------------------------------------------------------------
MObjBuilder::MObjBuilder()
{
	static bool bfirst = true;
	if (bfirst) init_function_lists();
	bfirst = false;

	m_autoVars = true;
}

//-----------------------------------------------------------------------------
char* find_next(char* sz)
{
	char* ch = sz;
	int l = 0;
	while (*ch)
	{
		if      (*ch == '(') l++;
		else if (*ch == ')') l--;
		else if ((*ch == ',') && (l==0)) return ch;
		++ch;
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool string_replace(string& s, const std::string& in, const std::string& out)
{
	size_t n = s.find(in);
	if (n != string::npos)
	{
		s.replace(n, in.size(), out);
		return true;
	}
	else return false;
}

//-----------------------------------------------------------------------------
// process string substitutions
// Math expressions can be defined with string literals, which will be substituted
// in the last statement.
// e.g.
// s = a + b; 2*s
// will be replaced by: s*(a+b)

std::string MObjBuilder::processStrings(const std::string& ex)
{
	// remove white space
	string tmp;
	for (size_t i = 0; i < ex.size(); ++i)
	{
		char c = ex[i];
		if (isspace(c) == 0) tmp.push_back(c);
	}

	// keep track of string subs
	std::map<string, string> subs;

	// build string sub list
	size_t lp = 0;
	size_t l1 = tmp.find(';', lp);
	while (l1 != std::string::npos)
	{
		size_t equalSign = tmp.find('=', lp);
		string leftVal = tmp.substr(lp, equalSign - lp);
		string rightVal = tmp.substr(equalSign + 1, l1 - equalSign - 1);
		subs[leftVal] = rightVal;

		lp = l1 + 1;
		l1 = tmp.find(';', lp);
	}

	string out = tmp.substr(lp);

	bool b = true;
	while (b)
	{
		std::map<string, string>::iterator it = subs.begin();
		b = false;
		for (; it != subs.end(); ++it)
		{
			string rep = "(" + it->second + ")";
			bool ret = string_replace(out, it->first, rep);
			if (ret) b = true;
		}
	}

	return out;
}


//-----------------------------------------------------------------------------
bool MObjBuilder::Create(MSimpleExpression* mo, const std::string& expression, bool beval)
{
	// process the strins for string substitutions
	std::string ex = processStrings(expression);

	// make a copy of the original string
	char szcopy[512] = { 0 };
	strcpy(szcopy, ex.c_str());
	m_szexpr = m_szorg = szcopy;

	// keep pointer to object being constructed
	m_po = mo;

	// let's build a simple expression
	try {
		MITEM i = create();
		if (beval) i = MEvaluate(i);
		mo->SetExpression(i);
	}
	catch (MathError e)
	{
		fprintf(stderr, "Error evaluating math expression: %s (position %d)\n", e.GetErrorStr(), e.GetPosition());
		return false;
	}
	catch (...)
	{
		fprintf(stderr, "Error evaluating math expression: (unknown)\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
// Convert a string to a math object
MathObject* MObjBuilder::Create(const std::string& ex, bool beval)
{
	// make a copy of the original string
	char szcopy[512] = {0};
	strcpy(szcopy, ex.c_str());
	m_szexpr = m_szorg = szcopy;

	m_po = 0;

	// let's build a simple expression
	MSimpleExpression* pe = new MSimpleExpression;
	m_po = pe;
	MITEM i = create();
	if (beval) i = MEvaluate(i);
	pe->SetExpression(i);

	return m_po;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::sequence()
{
	MSequence* ps = new MSequence;
	do
	{
		MItem* pi = create();
		ps->add(pi);
	}
	while (curr_tok == COMMA);

	return ps;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::create()
{
	MItem* pi = expr();
	for (;;)
	{
		switch(curr_tok)
		{
		case EQUAL:
			pi = new MEquation(pi, expr());
			break;
		case EQUALITY:
			pi = new MEquality(pi, expr());
			break;
/*		case COMMA:
			{
				MSequence* pe = dynamic_cast<MSequence*>(pi);
				if (pe == 0)
				{
					pe = new MSequence;
					pe->add(pi);
				}
				pe->add(expr());
				pi = pe;
			}
			break;
*/		default:
			return pi;
		}
	}
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::create_sequence()
{
	MItem* pi = expr();
	for (;;)
	{
		switch(curr_tok)
		{
		case EQUAL:
			pi = new MEquation(pi, expr());
			break;
		case COMMA:
			{
				MSequence* pe = dynamic_cast<MSequence*>(pi);
				if (pe == 0)
				{
					pe = new MSequence;
					pe->add(pi);
				}
				pe->add(expr());
				pi = pe;
			}
			break;
		default:
			return pi;
		}
	}
}


//-----------------------------------------------------------------------------
MItem* MObjBuilder::expr()
{
	MItem* pi = term();
	for (;;)
	{
		switch(curr_tok)
		{
		case PLUS:
			pi = new MAdd(pi, term());
			break;
		case MINUS:
			pi = new MSub(pi, term());
			break;
		default:
			return pi;
		}
	}
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::term()
{
	MItem* pi = power();
	for (;;)
		switch(curr_tok)
		{
		case MUL:
			pi = new MMul(pi, power());
			break;
		case DIV:
			pi = new MDiv(pi, power());
			break;
		case CONTRACT:
			{
				MMatrix* pl = mmatrix(pi);
				if (pl == 0) throw MathError(Position(), "invalid left term");
				MMatrix* pr = mmatrix(power());
				if (pr == 0) throw MathError(Position(), "invalid right term");
				pi = new MFuncMat2(matrix_contract, "contract", pl, pr);
			}
			break;
		default:
			return pi;
		}
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::power()
{
	MItem* pi = prim();
	for (;;)
		switch(curr_tok)
		{
		case POW:
			pi = new MPow(pi, prim());
			break;
		default:
			return pi;
		}
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::prim()
{
	MItem* pi = 0;
	get_token();

	switch (curr_tok)
	{
	case NUMBER:
		{
			pi = new MConstant(number_value);
			get_token();
		}
		break;
	case NAME:
		{
			string s(string_value);
			get_token();

			// if the next token is a left bracket, let's assume it's a function
			if (curr_tok == LP)
			{
				pi = func();
				if (pi == 0) throw MathError(Position(), "Unknown function"); 
			}
			else
			{
				pi = var();
				if (is_matrix(pi))
				{
					// we might be extracting the component of a matrix
					if (curr_tok == LB)
					{
						MSequence* pe = dynamic_cast<MSequence*>(sequence());
						read_token(RB);

						if (pe == 0) { delete pe; throw MathError(Position(), "comma expected"); }
						MSequence& e = *pe;
						if (e.size() != 2) throw MathError(Position(), "syntax error");
						const MConstant* pc0 = mconst(e[0]); if ((pc0 == 0) || (is_int(pc0->value()) == false)) throw MathError(Position(), "syntax error");
						const MConstant* pc1 = mconst(e[1]); if ((pc1 == 0) || (is_int(pc1->value()) == false)) throw MathError(Position(), "syntax error");
						int i = (int)(pc0->value());
						int j = (int)(pc1->value());

						const MMatrix& m = *mmatrix(pi);
						int nr = m.rows();
						int nc = m.columns();
						if ((i<0)||(i>=nr)) throw MathError(Position(), "invalid matrix component");
						if ((j<0)||(j>=nc)) throw MathError(Position(), "invalid matrix component");

						pi = (m[i][j]->copy());
						delete pe;
					}
				}
			}
		}
		break;
	case MINUS:
		{
			pi = new MNeg(power());
		}
		break;
	case LP:
		{
			pi = create();
			read_token(RP); // eat ')'
		}
		break;
	case LC:
		{
			pi = sequence();
			if (pi == 0) { throw MathError(Position(), "comma expected"); }
			read_token(RC);
		}
		break;
	case LB:
		{
			MSequence* pe = dynamic_cast<MSequence*>(sequence());
			if (pe == 0) { delete pe; throw MathError(Position(), "comma expected"); }
			MSequence& e = *pe;
			read_token(RB);
			if (is_matrix(e[0]))
			{
				int nrow = e.size();
				int ncol = mmatrix(e[0])->columns();
				int N = e.size();
				for (int i=0; i<N; ++i)
				{
					const MMatrix* pm = mmatrix(e[i]);
					if (pm == 0) throw MathError(Position(), "invalid matrix definition");
					if (pm->columns() != ncol) throw MathError(Position(), "invalid matrix definition");
					if (pm->rows() != 1) throw MathError(Position(), "invalid matrix definition");
				}

				MMatrix* pm = new MMatrix();
				pm->Create(nrow, ncol);
				for (int i=0; i<nrow; ++i) 
					for (int j=0; j<ncol; ++j)
					{
						const MMatrix* pk = mmatrix(e[i]);
						(*pm)[i][j] = (*pk)[0][j]->copy();
					}
				pi = pm;
			}
			else
			{
				MMatrix* pm = new MMatrix();
				int ncol = e.size();
				pm->Create(1, ncol);
				for (int i=0; i<ncol; ++i) (*pm)[0][i] = e[i]->copy();
				pi = pm;
			}
			delete pe;
		}
		break;
	case AB:
		{
			pi = new MFunc1D(fabs, "abs", create());
			if (curr_tok != AB) throw MathError(Position(), "'|' expected");
			get_token(); // eat '|'
		}
		break;
	default:
		throw MathError(Position(), "Token expected");
	}
	assert(pi);


	// look for post-operators
	do 
	{
		if (curr_tok == FACT)
		{
			pi = new MFunc1D(fac, "fac", pi);
			get_token(); // eat '!'
		}
		else if (curr_tok == TRANS)
		{
			pi = new MFuncMat(matrix_transpose, "transpose", pi);
			get_token(); // eat '''
		}
		else break;
	}
	while (true);

	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::var()
{
	string s(string_value);
	MItem* pi = 0;
	// see if it is a predefined variable
	map<string, double>::iterator pc = CNT.find(s);
	if (pc != CNT.end())
	{
		pi = new MNamedCt(pc->second, pc->first);
	}
	else
	{
		// see if it is a user variable
		MVariable* pvar = m_po->FindVariable(s);
		if (pvar == 0)
		{
			if (m_autoVars == false) throw MathError(Position(), "Unknown variable");
			// create a new variable
			pvar = new MVariable(s);
			m_po->AddVariable(pvar);
		}
		pi = new MVarRef(pvar);
	}
	return pi;
}

//-----------------------------------------------------------------------------
// parse a function
MItem* MObjBuilder::func()
{
	MItem* pi = 0;
	pi = fncnd();
	if (pi == 0) pi = fnc1d();
	if (pi == 0) pi = fnc2d();
	if (pi == 0) pi = fmat();
	if (pi == 0) pi = fsym();
	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::fnc1d()
{
	MItem* pi = 0;
	string s(string_value);
	map<string, FUNCPTR>::iterator pf = FNC.find(s);
	if (pf != FNC.end())
	{
		FUNCPTR fnc = pf->second;
		if (curr_tok != LP) throw MathError(Position(), "'(' expected");
		pi = new MFunc1D(fnc, pf->first, create());
		if (curr_tok != RP) throw MathError(Position(), "')' expected");
		get_token();
	}
	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::fnc2d()
{
	MItem* pi = 0;
	string s(string_value);
	map<string, FUNC2PTR>::iterator pf = FNC2.find(s);
	if (pf != FNC2.end())
	{
		FUNC2PTR fnc = pf->second;
		MItem *p1, *p2;
		p1 = create();
		if (curr_tok != COMMA) throw MathError(Position(), "',' expected");
		p2 = create();
		pi = new MFunc2D(fnc, pf->first, p1, p2);
		if (curr_tok != RP) throw MathError(Position(), "')' expected");
		get_token();
	}
	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::fncnd()
{
	MItem* pi = 0;
	string s(string_value);
	map<string, FUNCNPTR>::iterator pf = FNCN.find(s);
	if (pf != FNCN.end())
	{
		FUNCNPTR fnc = pf->second;
		MSequence pl;
		do
		{
			pl.add(create());
		}
		while (curr_tok == COMMA);
		if (curr_tok != RP) throw MathError(Position(), "')' expected");
		get_token();
		pi = new MFuncND(fnc, pf->first, pl);
	}
	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::fmat()
{
	MItem* pi = 0;
	string s(string_value);
	map<string, FUNCMATPTR>::iterator pf = FMAT.find(s);
	if (pf != FMAT.end())
	{
		FUNCMATPTR fnc = pf->second;
		if (curr_tok != LP) throw MathError(Position(), "'(' expected");
		pi = new MFuncMat(fnc, pf->first, create());
		if (curr_tok != RP) throw MathError(Position(), "')' expected");
		get_token();
	}
	return pi;
}

//-----------------------------------------------------------------------------
// Read a symbolic function.
MItem*	MObjBuilder::fsym()
{
	string s(string_value);
	if      (s.compare("derive"   ) == 0) return derive   ();
	else if (s.compare("replace"  ) == 0) return replace  ();
	else if (s.compare("taylor"   ) == 0) return taylor   ();
	else if (s.compare("integrate") == 0) return integrate();
	else if (s.compare("expand"   ) == 0) return expand   ();
	else if (s.compare("simplify" ) == 0) return simplify ();
	else if (s.compare("solve"    ) == 0) return solve    ();
	else if (s.compare("collect"  ) == 0) return collect  ();
	else if (s.compare("identity" ) == 0) return identity ();
	else if (s.compare("mdiag"    ) == 0) return mdiag    ();
	return 0;
}

//-----------------------------------------------------------------------------
MObjBuilder::Token_value MObjBuilder::get_token()
{
	// remove leading whitespace
	while ((*m_szexpr==' ') || (*m_szexpr=='\t')) m_szexpr++;

	// get the first character
	char ch = *m_szexpr++;

	switch(ch)
	{
	case 0:
		return curr_tok = END;
	case '^':
	case '*':
	case '/':
	case '+':
	case '-':
	case '(':
	case ')':
	case '[':
	case ']':
	case '{':
	case '}':
	case '%':
	case ',':
	case '|':
	case '!':
	case '\'':
	case ':':
		return curr_tok = Token_value(ch);
	case '=':
		if (*m_szexpr=='=') { m_szexpr++; return curr_tok = EQUALITY; } else return curr_tok = EQUAL;
		break;
	case '0': case '1': case '2': case '3': case '4':
	case '5': case '6': case '7': case '8': case '9':
	case '.':
		m_szexpr--;
		number_value = get_number();
		return curr_tok = NUMBER;
	case '$':
		{
		if (*m_szexpr != '{') throw MathError(Position(), "'{' expected");
		m_szexpr++;
		int n = 0;
		string_value[0] = 0;
		while ((*m_szexpr) && (*m_szexpr != '}')) string_value[n++] = *m_szexpr++;
		string_value[n] = 0;
		if (*m_szexpr != '}') throw MathError(Position(), "'}' expected");
		m_szexpr++;
		}
		return curr_tok = NAME;
	default:
		if (isalpha(ch)||(ch=='_'))
		{
			m_szexpr--;
			get_name(string_value);
			return curr_tok = NAME;
		}
		throw MathError(Position(), "Bad token");
		return curr_tok=PRINT;
	}
}

double MObjBuilder::get_number()
{	
	const char* ch = m_szexpr;

	// read first digits
	while (isdigit(*ch)) ch++;
	if (*ch == '.')
	{
		// read more digits after decimal point
		ch++;
		while (isdigit(*ch)) ch++;
	}

	// see if we're using exp notation
	if ((*ch == 'E') || (*ch == 'e'))
	{
		ch++;
		if ((*ch == '-') || (*ch == '+')) ch++;

		// read digits
		while (isdigit(*ch)) ch++;
	}

	double val = atof(m_szexpr);
	m_szexpr = ch;

	return val;
}

void MObjBuilder::get_name(char* str)
{
	int n = 0;
	bool binside = false;
	while (isalnum(*m_szexpr) || (*m_szexpr == '_') || (*m_szexpr == '.') || (*m_szexpr == '{') || ((*m_szexpr == '}')&&binside) || binside)
	{
		if (*m_szexpr == '{') binside = true;
		if (*m_szexpr == '}') binside = false;
		str[n++] = *m_szexpr++;
	}
	str[n] = 0;
}

//-----------------------------------------------------------------------------
// Check's to see if the current token is
void MObjBuilder::read_token(MObjBuilder::Token_value n, bool bnext)
{
	switch (n)
	{
	case COMMA : if (curr_tok != COMMA ) throw MathError(Position(), "',' expected"); break;
	case NAME  : if (curr_tok != NAME  ) throw MathError(Position(), "syntax error"); break;
	case NUMBER: if (curr_tok != NUMBER) throw MathError(Position(), "syntax error"); break;
	case RP    : if (curr_tok != RP    ) throw MathError(Position(), "')' expected"); break;
	case LC    : if (curr_tok != LC    ) throw MathError(Position(), "'{' expected"); break;
	case RC    : if (curr_tok != RC    ) throw MathError(Position(), "'}' expected"); break;
	case RB    : if (curr_tok != RB    ) throw MathError(Position(), "']' expected"); break;
	}

	if (bnext) get_token();
}

//-----------------------------------------------------------------------------
// Read a variable name; assuming the current token is indeed a variable
MVariable* MObjBuilder::read_var(bool badd)
{
	// read the name token
	read_token(NAME);

	// get the variable
	MVariable* pv = m_po->FindVariable(string_value);
	if (pv == 0)
	{
		if (badd)
		{
			pv = new MVariable(string_value);
			m_po->AddVariable(pv);
		}
		else throw MathError(Position(), "unknown variable"); 
	}
	return pv;
}

//-----------------------------------------------------------------------------
// Read an integer number
int MObjBuilder::read_int()
{
	read_token(NUMBER);
	double f = number_value;
	if (f - floor(f) != 0) throw MathError(Position(), "integer expected");
	return (int) f;
}

//-----------------------------------------------------------------------------
// Read a double number
double MObjBuilder::read_double()
{
	read_token(NUMBER);
	return number_value;
}

//-----------------------------------------------------------------------------
// read a mathematical expression
MItem* MObjBuilder::read_math()
{
	return create();
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::derive()
{
	// read the argument
	MITEM e = read_math();

	// read the comma
	read_token(COMMA, false);

	// find the variable
	MITEM v = read_math();

	MItem* pi = 0;
	if (is_sequence(v))
	{
		// read the parentheses
		read_token(RP);

		// calculate derivatives
		const MSequence& s = *msequence(v);
		pi = MDerive(e, s).copy();
	}
	else
	{
		// make sure v is a variable
		if (is_var(v) == false) throw MathError(Position(), "This is not a variable");
		const MVariable& x = *(mvar(v)->GetVariable());
		
		// read optional arguments
		if (curr_tok == COMMA)
		{
			// skip the comma
			get_token();

			if (curr_tok == NUMBER)
			{
				// read the number
				int n = read_int();

				// make sure it is at least one
				if (n <= 0) throw MathError(Position(), "integer must be at least one");

				read_token(RP);

				pi = MDerive(e, x, n).copy();
			}
			else MathError(Position(), "number unexpected");
		}
		else
		{
			read_token(RP);
			pi = MDerive(e, x).copy();
		}
	}

	// return the result
	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::replace()
{
	// read the expression
	MITEM e = read_math();

	// read the comma
	read_token(COMMA, false);

	// read the expression to replace
	MITEM x = read_math();

	if (is_sequence(x))
	{
		const MSequence& v = *msequence(x);

		// read the comma
		read_token(COMMA);

		// read the left curly bracket
		read_token(LC, false);

		// read the list of new variables
		const MSequence& s = *msequence(sequence());

		// read the right curly bracket
		read_token(RC);

		// read the right parenthesis
		read_token(RP);

		// calculate the new expression
		MItem* pi = MReplace(e, v, s).copy();

		// clean-up
		delete &s;

		// done
		return pi;
	}
	else
	{
		// read comma
		if (curr_tok == COMMA)
		{
			// read the second expression
			MITEM s = read_math();

			// read the right parenthesis
			read_token(RP);

			// call replace function
			return MReplace(e, x, s).copy();
		}
		else
		{
			read_token(RP);
			return MReplace(e,x).copy();
		}
	}
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::taylor()
{
	// read the expression
	MITEM e = read_math();

	// read comma
	read_token(COMMA);

	// read the variable name
	MVariable* pv = read_var();

	// read comma
	read_token(COMMA);

	// read the point at which to expand the expression
	double z = read_double();

	// read comma
	read_token(COMMA);

	// read the number of terms in taylor expansion
	int n = read_int();
	if (n < 0) throw MathError(Position(), "integer must be larger than zero");

	// read the right parenthesis
	read_token(RP);

	return MTaylor(e, *pv, z, n).copy();
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::integrate()
{
	// read the expression
	MITEM e = read_math();

	// read comma
	read_token(COMMA);

	// read the variable
	MVariable* pv = read_var(true);

	MItem* pi = 0;
	if (curr_tok == COMMA)
	{
		// if a comma is there we calculate the definite integral
		MITEM a = read_math();

		read_token(COMMA, false);

		MITEM b = create();

		read_token(RP);

		pi = MIntegral(e, *pv, a, b).copy();
	}
	else
	{
		// else we calculate the indefinite integral
		read_token(RP);

		pi = MIntegral(e, *pv).copy();
	}

	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::expand()
{
	// read the expression
	MITEM e = read_math();

	MItem* pi = 0;
	if (curr_tok == COMMA)
	{
		MITEM s = read_math();
		read_token(RP);
		pi = MExpand(e, s).copy();
	}
	else 
	{
		read_token(RP);
		pi = MExpand(e).copy();
	}
	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::simplify()
{
	// read the expression
	MITEM e = read_math();
	read_token(RP);

	MItem* pi = MExpand(e).copy();
	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::solve()
{
	// read the argument
	MITEM e = read_math();

	// read the comma
	read_token(COMMA, false);

	// find the variable
	MITEM v = read_math();

	read_token(RP);

	MItem* pi = MSolve(e, v).copy();

	// return the result
	return pi;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::collect()
{
	// read the expression
	MITEM e = read_math();
	read_token(COMMA, false);
	MITEM a = read_math();
	read_token(RP);
	return MCollect(e, a).copy();
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::identity()
{
	// read the expression
	MITEM e = read_math();
	read_token(RP);
	if (is_int(e))
	{
		return matrix_identity((int)e.value());
	}
	else throw MathError(Position(), "constant expected");
	return 0;
}

//-----------------------------------------------------------------------------
MItem* MObjBuilder::mdiag()
{
	// read the expression
	MSequence* pe = dynamic_cast<MSequence*>(sequence());
	read_token(RP);
	MSequence& e = *pe;

	int n = pe->size();
	MMatrix* pm = new MMatrix();
	pm->Create(n, n);
	for (int i=0; i<n; ++i) 
		for (int j=0; j<n; ++j)
			{
				if (i == j)
				{
					const MItem* pi = e[i];
					(*pm)[i][j] = pi->copy();
				}
				else
				{
					(*pm)[i][j] = new MConstant(0.0);
				}
			}
	delete pe;
	return pm;
}
