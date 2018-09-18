#include "MMath.h"
#include "MEvaluate.h"
#include "MathObject.h"
#include <map>
#include <string>
using namespace std;

//-----------------------------------------------------------------------------
MITEM MReplace(MITEM& e, MITEM& x)
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
MITEM MReplace(MITEM& e, MVariable& x, MITEM& s)
{
	MITEM v(x);
	return MReplace(e,v,s);
}

//-----------------------------------------------------------------------------
MITEM MReplace(MITEM& e, MITEM& x, MITEM& s)
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
			MFunc1D* pf = mfnc1d(e);
			MITEM p = pf->Item()->copy();
			return new MFunc1D(pf->funcptr(), pf->Name(), MReplace(p, x, s).copy());
		}
	case MF2D:
		{
			MFunc2D* pf = mfnc2d(e);
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
			MMatrix* pm = mmatrix(e);
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
MITEM MReplace(MITEM& e, MSequence& x, MSequence& s)
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
