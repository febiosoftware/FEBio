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



#pragma once

#include "MathObject.h"
#include "MMatrix.h"
#include <list>

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
	std::list<MITEM>		m_p;
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
	std::list<MTerm>	m_t;
	FRACTION	m_c;	// accumulator for constants
};
