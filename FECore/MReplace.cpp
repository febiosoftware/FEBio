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
#include "MathObject.h"
#include <map>
#include <string>
using namespace std;

//-----------------------------------------------------------------------------
MITEM MReplace(const MITEM& e, const MITEM& x)
{
	if (is_equation(x))
	{
		MITEM l = x.Left();
		MITEM r = x.Right();
		return MReplace(e, l, r);
	}
	else return e;
}

//-----------------------------------------------------------------------------
MITEM MReplace(const MITEM& e, const MVariable& x, const MITEM& s)
{
	MITEM v(x);
	return MReplace(e,v,s);
}

//-----------------------------------------------------------------------------
MITEM MReplace(const MITEM& e, const MITEM& x, const MITEM& s)
{
	// see if they are equal
	if (e == x) return s;

	switch (e.Type())
	{
	case MCONST:
	case MFRAC:
	case MNAMED:
	case MVAR:
		return e;
	case MNEG:
		return MEvaluate(-(MReplace(-e, x, s)));
	case MADD:
		{
			MITEM l = MReplace(e.Left (), x, s);
			MITEM r = MReplace(e.Right(), x, s);
			return MEvaluate(l+r);
		}
	case MSUB:
		{
			MITEM l = MReplace(e.Left (), x, s);
			MITEM r = MReplace(e.Right(), x, s);
			return MEvaluate(l-r);
		}
	case MMUL:
		{
			MITEM l = MReplace(e.Left (), x, s);
			MITEM r = MReplace(e.Right(), x, s);
			return MEvaluate(l*r);
		}
	case MDIV:
		{
			MITEM l = MReplace(e.Left (), x, s);
			MITEM r = MReplace(e.Right(), x, s);
			return MEvaluate(l/r);
		}
	case MPOW:
		{
			MITEM l = MReplace(e.Left (), x, s);
			MITEM r = MReplace(e.Right(), x, s);
			return MEvaluate(l^r);
		}
	case MEQUATION:
		{
			MITEM l = MReplace(e.Left (), x, s);
			MITEM r = MReplace(e.Right(), x, s);
			return MITEM(new MEquation(l.copy(), r.copy()));
		}
	case MF1D:
		{
			const MFunc1D* pf = mfnc1d(e);
			MITEM p = pf->Item()->copy();
			return new MFunc1D(pf->funcptr(), pf->Name(), MReplace(p, x, s).copy());
		}
	case MF2D:
		{
			const MFunc2D* pf = mfnc2d(e);
			MITEM l = MReplace(e.Left (), x, s);
			MITEM r = MReplace(e.Right(), x, s);
			return new MFunc2D(pf->funcptr(), pf->Name(), l.copy(), r.copy());
		}
/*	case MSFND:
		{
			MSFuncND* pf = msfncnd(e);
			MITEM v(pf->Item()->copy());
			MITEM vnew = MReplace(v, x, s);

			string s = pf->Name();
			M1DFuncDef* pfd = DEF1D.find(s)->second;
			assert(pfd);
			return pfd->CreateFunction(vnew.ItemPtr());
		}
*/	case MMATRIX:
		{
			const MMatrix* pm = mmatrix(e);
			int nrows = pm->rows();
			int ncols = pm->columns();
			MMatrix* pmnew = new MMatrix;
			pmnew->Create(nrows, ncols);
			for (int i=0; i<nrows; ++i)
				for (int j=0; j<nrows; ++j)
				{
					MITEM mij = pm->Item(i,j)->copy();
					(*pmnew)[i][j] = MReplace(mij, x, s).copy();
				}
			return pmnew;
		}
	}

	assert(false);
	return e;
}

//-----------------------------------------------------------------------------
// replace multiple expressions with other expressions
MITEM MReplace(const MITEM& e, const MSequence& x, const MSequence& s)
{
	if (x.size() != s.size()) return e;

	MITEM r = e;
	const int N = x.size();
	for (int i=0; i<N; ++i)
	{
		MITEM xi = x[i]->copy();
		MITEM si = s[i]->copy();
		r = MReplace(r, xi, si);
	}

	return r;
}
