#pragma once

#include "MathObject.h"

class MathError
{
public:
	MathError(int n, const char* sz) : m_npos(n), m_szerr(sz){}
	int GetPosition() { return m_npos; }
	const char* GetErrorStr() { return m_szerr; }

protected:
	int	m_npos;
	const char*	m_szerr;
};


class MObjBuilder
{
private:
	enum Token_value {
		NAME, NUMBER, END,
		EQUAL='=', PLUS='+', MINUS='-', MUL='*', DIV='/', MOD = '%', POW='^', FACT='!', TRANS = '\'',
		LP='(',	RP=')', AB = '|', LB = '[', RB = ']', LC='{', RC='}',
		CONTRACT = ':',
		COMMA = ',', PRINT, EQUALITY
	};

public:
	MObjBuilder();

	bool Create(MSimpleExpression* mo, const std::string& ex, bool eval);
	MathObject* Create(const std::string& ex , bool eval);

	static void Clear();

	void setAutoVars(bool b) { m_autoVars = b; }

protected:
	MItem*	create();
	MItem*	create_sequence();
	MItem*	expr();
	MItem*	term ();
	MItem*	power();
	MItem*	prim ();
	MItem*  var  ();
	MItem*  func ();
	MItem*	fnc1d();
	MItem*	fnc2d();
	MItem*	fncnd();
	MItem*	fmat ();
	MItem*	fsym();
	MItem*	sequence();

	double get_number();
	void get_name(char* str);

	int Position() { return (int)(m_szexpr - m_szorg); }

protected:
	MItem*	derive   ();
	MItem*	replace  ();
	MItem*	taylor   ();
	MItem*	integrate();
	MItem*	expand   ();
	MItem*	simplify ();
	MItem*	solve    ();
	MItem*	collect  ();
	MItem*  identity ();
	MItem*  mdiag    ();

protected:
	MathObject*	m_po;

	Token_value get_token();
	Token_value	curr_tok;

	// read an expression
	MItem* read_math();

	// read a required token
	void read_token(Token_value n, bool bnext = true);

	// read a variable name
	MVariable* read_var(bool badd = false);

	// read an int
	int read_int();

	// read a double
	double read_double();

protected:
	const char*	m_szexpr, *m_szorg;
	double		number_value;
	char		string_value[256];

	bool		m_autoVars;	// add new variables automatically
};
