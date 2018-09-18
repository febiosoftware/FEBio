#pragma once

#include "MathObject.h"
#include "MMatrix.h"
#include <stack>
#include <list>
using namespace std;

//-----------------------------------------------------------------------------
void MEvaluate(MathObject* po);

//-----------------------------------------------------------------------------
MITEM MEvaluate(const MITEM& e);

//-----------------------------------------------------------------------------
MITEM MEvaluate(const MMatrix& A);

//-----------------------------------------------------------------------------
MITEM MMultiply(const MITEM& l, const MITEM& r);

//-----------------------------------------------------------------------------
MITEM MDivide(const MITEM& n, const MITEM& d);

//-----------------------------------------------------------------------------
MITEM MAddition(const MITEM& l, const MITEM& r);

//-----------------------------------------------------------------------------
// This class describes a multi-factor product. It is used to simplify 
// expressions involving factors for which the binary operations might be
// too difficult to process
class MProduct
{
public:
	MProduct(const MITEM& a);

public:
	void Multiply(const MITEM& a);
	MITEM operator / (const MProduct& D);

	MITEM Item();

	bool operator == (const MProduct& a);

	bool contains(const MITEM& i) const;

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
		MTerm(const MITEM& i);
		MTerm(const MTerm& i) { m_a = i.m_a; m_s = i.m_s; }
		void operator = (const MTerm& i) { m_a = i.m_a; m_s = i.m_s; }
	};

public:
	MSum(const MITEM& a);
	void Add(const MITEM& a);
	void Sub(const MITEM& a);

	MITEM Item();

protected:
	list<MTerm>	m_t;
	FRACTION	m_c;	// accumulator for constants
};
