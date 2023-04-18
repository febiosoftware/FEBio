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
MITEM MExpand(const MITEM& i)
{
	return MExpand(i, MITEM((MItem*)0));
}

//-----------------------------------------------------------------------------
MITEM MExpand(const MITEM& i, const MITEM& s)
{
	if (i == s) return i;

	MITEM e = MEvaluate(i);

	switch (e.Type())
	{
	case MFRAC:
		{
			FRACTION f = mfrac(i)->fraction();
			if (f.d == 1.0) return f.n;
		}
		break;
	case MNEG: return -MExpand(e.Item(),s);
	case MMUL:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_add(r)&&(s != r)) 
			{
				MITEM a = r.Left();
				MITEM b = r.Right();
				return MExpand(l*a,s) + MExpand(l*b,s);
			}
			if (is_sub(r)&&(s != r))
			{
				MITEM a = r.Left();
				MITEM b = r.Right();
				return MExpand(l*a,s) - MExpand(l*b,s);
			}
			if (is_add(l)&&(s != l)) 
			{
				MITEM a = l.Left();
				MITEM b = l.Right();
				return MExpand(r*a,s) + MExpand(r*b,s);
			}
			if (is_sub(l)&&(s != l))
			{
				MITEM a = l.Left();
				MITEM b = l.Right();
				return MExpand(r*a,s) - MExpand(r*b,s);
			}
			MITEM Ml = MExpand(l,s);
			MITEM Mr = MExpand(r,s);
//			if ((l != Ml) || (r != Mr)) return Ml*Mr;
			if ((l != Ml) || (r != Mr)) return MExpand(Ml*Mr,s); // it looks like this could create an infinite loop
		}
		break;
	case MDIV:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_add(l)&&(s != l))
			{
				MITEM a = l.Left();
				MITEM b = l.Right();
				return MExpand(a/r,s) + MExpand(b/r,s);
			}
			if (is_sub(l)&&(s != l))
			{
				MITEM a = l.Left();
				MITEM b = l.Right();
				return MExpand(a/r,s) - MExpand(b/r,s);
			}
			MITEM Ml = MExpand(l,s);
			MITEM Mr = MExpand(r, s);
			if ((l != Ml) || (r != Mr)) return MExpand(Ml/Mr,s); else return Ml/Mr;
		}
		break;
	case MADD:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			return MExpand(l,s) + MExpand(r,s);
		}
		break;
	case MSUB:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			return MExpand(l,s) - MExpand(r,s);
		}
		break;
	case MPOW:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_int(r) && is_add(l) && (l != s))
			{
				MITEM a = l.Left();
				MITEM b = l.Right();
				int n = (int) r.value();
				if (n == 0) return 1.0;
				MITEM sum(0.0);
				for (int i=0; i<=n; ++i)
				{
					double C = binomial(n,i);
					sum = sum + C*MExpand(a^(n-i),s)*MExpand(b^i,s);
				}
				return sum;
			}
			if (is_int(r) && is_sub(l))
			{
				MITEM a = l.Left();
				MITEM b = l.Right();
				int n = (int) r.value();
				if (n == 0) return 1.0;
				MITEM sum(0.0);
				for (int i=0; i<=n; ++i)
				{
					double C = binomial(n,i);
					MITEM t = C*MExpand(a^(n-i),s)*MExpand(b^i,s);
					if (i%2 == 0) sum = sum + t;
					else sum = sum - t;
				}
				return sum;
			}
			if (is_add(r) && (r != s))
			{
				MITEM a = r.Left();
				MITEM b = r.Right();
				return ((l^a)*(l^b));
			}
			MITEM le = MExpand(l);
			MITEM re = MExpand(r);
			return le^re;
		}
		break;
	case MF1D:
		{
			const string& f = mfnc1d(e)->Name();
			if (f.compare("cos") == 0)
			{
				MITEM p = e.Param();
				if (p.Type() == MADD)
				{
					MITEM a = p.Left();
					MITEM b = p.Right();
					return MExpand(Cos(a)*Cos(b) - Sin(a)*Sin(b),s); 
				}
				if (p.Type() == MSUB)
				{
					MITEM a = p.Left();
					MITEM b = p.Right();
					return MExpand(Cos(a)*Cos(b) + Sin(a)*Sin(b),s); 
				}
			}
			if (f.compare("sin") == 0)
			{
				MITEM p = e.Param();
				if (p.Type() == MADD)
				{
					MITEM a = p.Left();
					MITEM b = p.Right();
					return MExpand(Sin(a)*Cos(b) + Cos(a)*Sin(b),s); 
				}
				if (p.Type() == MSUB)
				{
					MITEM a = p.Left();
					MITEM b = p.Right();
					return MExpand(Sin(a)*Cos(b) - Cos(a)*Sin(b),s); 
				}
			}
			if (f.compare("tan") == 0)
			{
				MITEM p = e.Param();
				if (p.Type() == MADD)
				{
					MITEM a = p.Left();
					MITEM b = p.Right();
					return MExpand((Tan(a)+Tan(b))/(1.0 - Tan(a)*Tan(b)),s); 
				}
				if (p.Type() == MSUB)
				{
					MITEM a = p.Left();
					MITEM b = p.Right();
					return MExpand((Tan(a)-Tan(b))/(1.0 + Tan(a)*Tan(b)),s); 
				}
			}
/*			if (f.compare("ln") == 0)
			{
				MITEM p = e.Param();
				if (is_mul(p))
				{
					MITEM a = MExpand(p.Left());
					MITEM b = MExpand(p.Right());
					return MExpand(Log(a)+Log(b));
				}
				if (is_div(p))
				{
					MITEM a = MExpand(p.Left());
					MITEM b = MExpand(p.Right());
					return MExpand(Log(a)-Log(b));
				}
				if (is_pow(p))
				{
					MITEM a = MExpand(p.Left());
					MITEM b = MExpand(p.Right());
					return MExpand(b*Log(a));
				}
			}
*/			const MFunc1D* pf = mfnc1d(e);
			MITEM v = MExpand(e.Param(), s);
			return new MFunc1D(pf->funcptr(), pf->Name(), v.copy());
		}
		break;
	case MSFNC:
		{
			MITEM f(msfncnd(i)->Value()->copy());
			return MExpand(f, s);
		}
		break;
	case MEQUATION:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			MITEM Ml = MExpand(l,s);
			MITEM Mr = MExpand(r,s);
			return new MEquation(Ml.copy(), Mr.copy());
		}
		break;
	case MSEQUENCE:
		{
			const MSequence& q = *msequence(i);
			MSequence* ps = new MSequence();
			for (int i=0; i<q.size(); ++i)
			{
				MITEM qi = q[i]->copy();
				MITEM vi = MExpand(qi, s);
				ps->add(vi.copy());
			}
			return ps;
		}
		break;
	case MMATRIX:
		{
			const MMatrix& m = *mmatrix(e);
			int ncol = m.columns();
			int nrow = m.rows();

			MMatrix* pdm = new MMatrix;
			pdm->Create(nrow, ncol);
			for (int i=0; i<nrow; ++i)
				for (int j=0; j<ncol; ++j)
				{
					MITEM mij(m.Item(i,j)->copy());
					MITEM aij = MExpand(mij, s);
					(*pdm)[i][j] = aij.copy();
				}
			return MITEM(pdm);
		}
		break;
	}
	return e;
}
