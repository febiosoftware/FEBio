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

#ifndef PI
#define PI	3.141592653589793
#endif

//-----------------------------------------------------------------------------
MITEM MDerive(const MITEM& a, const MVariable& x)
{
	MITEM e = MEvaluate(a);

	if (is_dependent(e, x) == false) return 0.0;

	switch (e.Type())
	{
	case MCONST:
	case MFRAC:
	case MNAMED: return 0.0; 
	case MVAR: if (e == x) return 1.0; else return 0.0;
	case MNEG: return -MDerive(-e, x); 
	case MADD: return (MDerive(e.Left(), x) + MDerive(e.Right(), x));
	case MSUB: 
		{
			MITEM l = MDerive(e.Left(), x);
			MITEM r = MDerive(e.Right(), x);
			return l - r;
		}
	case MMUL:
		{
			MITEM l = e.Left();
			MITEM dl = MDerive(l, x);
			MITEM r = e.Right();
			MITEM dr = MDerive(r, x);

			return (dl*r + l*dr);
		}
	case MDIV:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			MITEM dl = MDerive(l, x)*r;
			MITEM dr = l*MDerive(r, x);
			// be careful here that we are not subtracting two pointers
			MITEM a = dl/(r^2.0);
			MITEM b = dr/(r^2.0);
			return (a-b);
		}
	case MPOW:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_dependent(r, x) == false) return (r*(l^(r-1.0)))*MDerive(l, x);
			else if (is_dependent(l, x) == false) return (Log(l)*e)*MDerive(r, x);
			else
			{
				MITEM dl = MDerive(l, x);
				MITEM dr = MDerive(r, x);
				return (e*(dr*Log(l) + (r*dl)/l));
			}
		}
	case MF1D:
		{
			const string& s = mfnc1d(e)->Name();
			MITEM p = mfnc1d(e)->Item()->copy();
			MITEM dp = MDerive(p,x);
			if (s.compare("cos") == 0) return -Sin(p)*dp;
			if (s.compare("sin") == 0) return  Cos(p)*dp;
			if (s.compare("tan") == 0) return (Sec(p)^2.0)*dp;
			if (s.compare("sec") == 0) return (Sec(p)*Tan(p))*dp;
			if (s.compare("csc") == 0) return (-Csc(p)*Cot(p))*dp;
			if (s.compare("cot") == 0) return -((Csc(p)^2.0)*dp);
			if (s.compare("abs") == 0) return Sgn(p)*dp;
			if (s.compare("ln" ) == 0) return (dp/p);
			if (s.compare("log") == 0) return (dp/p)/Log(MITEM(10.0));
			if (s.compare("asin") == 0) return (dp/Sqrt(1.0 - (p^2)));
			if (s.compare("acos") == 0) return (-dp/Sqrt(1.0 - (p^2)));
			if (s.compare("atan") == 0) return (dp/((p^2) + 1.0));
			if (s.compare("cosh") == 0) return (dp*Sinh(p));
			if (s.compare("sinh") == 0) return (dp*Cosh(p));
			if (s.compare("sqrt") == 0) return (dp/(2.0*Sqrt(p)));
			if (s.compare("acosh") == 0) return (dp/Sqrt((p^2) - 1.0));
#ifdef WIN32
			if (s.compare("J0"  ) == 0) return (-J1(p))*dp;
			if (s.compare("J1"  ) == 0) return dp*(J0(p) - Jn(2, p))/2.0;
			if (s.compare("Y0"  ) == 0) return (-Y1(p))*dp;
			if (s.compare("Y1"  ) == 0) return dp*(Y0(p) - Yn(2, p))/2.0;
#endif
			if (s.compare("erf" ) == 0)
			{
				MITEM Pi = new MNamedCt(PI, "pi");
				return 2/Sqrt(Pi)*Exp(-(p^2))*dp;
			}
			if (s.compare("erfc") == 0)
			{
				MITEM Pi = new MNamedCt(PI, "pi");
				return -(2/Sqrt(Pi))*Exp(-(p^2))*dp;
			}
			if (s.compare("H") == 0)
			{
				return 0.0;
			}
			assert(false);
		}
		break;
	case MF2D:
		{
			const string& s = mfnc2d(e)->Name();
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (s.compare("pow") == 0)
			{
				if (isConst(r) || is_named(r) || is_frac(r)) return (r * (l ^ (r - 1.0))) * MDerive(l, x);
				else if (isConst(l) || is_named(l) || is_frac(l)) return (Log(l) * e) * MDerive(r, x);
				else
				{
					MITEM dl = MDerive(l, x);
					MITEM dr = MDerive(r, x);
					return (e * (dr * Log(l) + (r * dl) / l));
				}
			}
#ifdef WIN32
			else if (s.compare("Jn") == 0)
			{
				MITEM dr = MDerive(r, x);
				if (is_int(l))
				{
					int n = (int) l.value();
					if (n==0) return (-J1(r))*dr;
					else return ((Jn(n-1, r) - Jn(n+1, r))/2.0)*dr;
				}
			}
			else if (s.compare("Yn") == 0)
			{
				MITEM dr = MDerive(r, x);
				if (is_int(l))
				{
					int n = (int) l.value();
					if (n==0) return (-Y1(r))*dr;
					else return ((Yn(n-1, r) - Yn(n+1, r))/2.0)*dr;
				}
			}
#endif
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
					MITEM dmij = MDerive(mij, x);
					(*pdm)[i][j] = dmij.copy();
				}
			return MITEM(pdm);
		}
		break;
	case MSFNC:
		{
			const MSFuncND& f = *msfncnd(e);
			MITEM v(f.Value()->copy());
			return MDerive(v, x);
		}
		break;
	}
	assert(false);
	return e;
}

//-----------------------------------------------------------------------------
MITEM MDerive(const MITEM& e, const MVariable& x, int n)
{
	MITEM d = e;
	for (int i=0; i<n; ++i) d = MDerive(d, x);
	return d;
}

//-----------------------------------------------------------------------------
MITEM MDerive(const MITEM& e, const MSequence& x)
{
	MITEM d = e;
	for (int i=0; i<(int) x.size(); ++i)
	{
		const MVariable& xi = *(mvar(x[i])->GetVariable());
		d = MDerive(d, xi);
	}
	return d;
}

//-----------------------------------------------------------------------------
MITEM MTaylor(const MITEM& e, const MVariable& v, double z, int n)
{
	MITEM a(z);
	MITEM x(v);

	MITEM dx = x - a;
	MITEM s = MReplace(e, v, a);
	MITEM t(e);
	double d = 1;

	for (int i=1; i<=n; ++i)
	{
		t = MDerive(t /d, v);
		d = (double) i;

		MITEM f = MReplace(t, v, a);
		MITEM Di = dx^d;
		MITEM ds = ((f/d)*Di);
		s = s + ds;
	}

	return s;
}
