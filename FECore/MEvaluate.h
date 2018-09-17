#pragma once

#include "MathObject.h"
#include "MMatrix.h"
#include <stack>
#include <list>
using namespace std;

//-----------------------------------------------------------------------------
void MEvaluate(MathObject* po);

//-----------------------------------------------------------------------------
MITEM MEvaluate(MITEM& e);

//-----------------------------------------------------------------------------
MITEM MEvaluate(MMatrix& A);

//-----------------------------------------------------------------------------
MITEM MMultiply(MITEM& l, MITEM& r);

//-----------------------------------------------------------------------------
MITEM MDivide(MITEM& n, MITEM& d);

//-----------------------------------------------------------------------------
MITEM MAddition(MITEM& l, MITEM& r);

//-----------------------------------------------------------------------------
// This class describes a multi-factor product. It is used to simplify 
// expressions involving factors for which the binary operations might be
// too difficult to process
class MProduct
{
public:
	MProduct(MITEM& a);

public:
	void Multiply(MITEM& a);
	MITEM operator / (MProduct& D);

	MITEM Item();

	bool operator == (MProduct& a);

	bool contains(MITEM& i);

protected:
	list<MITEM>		m_p;
};

//-----------------------------------------------------------------------------
// This class describes a general sum
class MSum
{
	class MTerm
	{
	public:
		MITEM		m_a;	// term
		FRACTION	m_s;	// multiplier (can be negative!)

	public:
		MTerm(MITEM& i);
		MTerm(const MTerm& i) { m_a = i.m_a; m_s = i.m_s; }
		void operator = (const MTerm& i) { m_a = i.m_a; m_s = i.m_s; }
	};

public:
	MSum(MITEM& a);
	void Add(MITEM& a);
	void Sub(MITEM& a);

	MITEM Item();

protected:
	list<MTerm>	m_t;
	FRACTION	m_c;	// accumulator for constants
};
