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
#include "MMath.h"
#include "MEvaluate.h"
using namespace std;

//-----------------------------------------------------------------------------
//! Calculate the definite integral of an expression
MITEM MIntegral(const MITEM& i, const MVariable& x, const MITEM& a, const MITEM& b)
{
	MITEM e = MEvaluate(i);
	MITEM Ie = MIntegral(e, x);

	return MReplace(Ie, x, b) - MReplace(Ie, x, a);
}

//-----------------------------------------------------------------------------
//! Calculate the indefinite integral of an expression.
//! Note that the integration constant is not included in the result
MITEM MIntegral(const MITEM& i, const MVariable& x)
{
	// simplify the expression
	MITEM e = MEvaluate(i);

	// see if there is a dependency on x
	if (is_dependent(e, x) == false) 
	{
		// if not, multiply by x and return
		return e*MITEM(x); 
	}

	// if we get here, there is a dependency on x so we
	// need to find an appropriate integration rule
	switch (e.Type())
	{
	case MVAR: 
		if (e == x) 
		{
			MITEM X(x);
			return Fraction(1.0, 2.0)*(X ^ 2.0); 
		}
		break;
	case MNEG: return -MIntegral(e.Item(), x); 
	case MADD: return MIntegral(e.Left(), x) + MIntegral(e.Right(), x); 
	case MSUB: return MIntegral(e.Left(), x) - MIntegral(e.Right(), x); 
	case MMUL:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_dependent(l, x) == false) return l*MIntegral(r, x);
			if (is_dependent(r, x) == false) return r*MIntegral(l, x);
			return MIntegral(MExpand(e), x);
		}
		break;
	case MDIV:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_dependent(r, x) == false) return MIntegral(l, x) / r;
		}
		break;
	case MPOW:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_var(l) && isConst(r))
			{
				if (l == x) 
				{
					if (r.value() != -1.0)
					{
						MITEM np1( r.value() + 1.0);
						return (l^np1)/np1;
					}
					else return Log(Abs(l));
				}
				else return e*MITEM(x);
			}
			if (is_int(r) && (is_add(l) || is_sub(l))) return MIntegral(MExpand(e),x);
			if (is_number(l) && (is_dependent(r,x)))
			{
				if ((is_int(l) == false) || (mnumber(l)->value() > 1.0))
				{
					if (r == x)
					{
						return e/Log(l);
					}
					else if (is_mul(r))
					{
						MITEM rl = r.Left();
						MITEM rr = r.Right();
						if (is_number(rl) && (rr == x))
						{
							return e/(rl*Log(l));
						}
					}
				}
			}
		}
		break;
	case MF1D:
		{
			const string& s = mfnc1d(e)->Name();
			MITEM p = mfnc1d(e)->Item()->copy();
			if (p == x)
			{
				if (s.compare("cos") == 0) return Sin(p);
				if (s.compare("sin") == 0) return -Cos(p);
				if (s.compare("tan") == 0) return -Log(Abs(Cos(p)));
				if (s.compare("cot") == 0) return Log(Abs(Sin(p)));
				if (s.compare("sec") == 0) return Log(Abs(Sec(p) + Tan(p)));
				if (s.compare("csc") == 0) return -Log(Abs(Csc(p) + Cot(p)));
				if (s.compare("sinh") == 0) return Cosh(p);
				if (s.compare("cosh") == 0) return Sinh(p);
				if (s.compare("tanh") == 0) return Log(Cosh(p));
				if (s.compare("sech") == 0) return Atan(Sinh(p));
				if (s.compare("exp") == 0) return Exp(p);
			}
		}
		break;
	}

	return new MOpIntegral(e.copy(), new MVarRef(&x));
}
