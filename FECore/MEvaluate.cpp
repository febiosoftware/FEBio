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
#include "MEvaluate.h"
#include "MMath.h"
#include <map>
using namespace std;

//-----------------------------------------------------------------------------
void MEvaluate(MathObject* po)
{
	if (dynamic_cast<MSimpleExpression*>(po))
	{
		MSimpleExpression* pe = dynamic_cast<MSimpleExpression*>(po);
		MITEM it = MEvaluate(pe->GetExpression());
		pe->SetExpression(it);
	}
}

//-----------------------------------------------------------------------------
MITEM MEvaluate(const MITEM& i)
{
	MITEM a;
	switch (i.Type())
	{
	case MNEG: a = -MEvaluate(i.Item()); break;
	case MADD:
		{
			MITEM l = MEvaluate(i.Left());
			MITEM r = MEvaluate(i.Right());
			a = MAddition(l, r);
		}
		break;
	case MSUB:
		{
			MITEM l = MEvaluate(i.Left());
			MITEM r = MEvaluate(i.Right());
			a = MAddition(l, -r);
		}
		break;
	case MMUL:
		{
			MITEM l = MEvaluate(i.Left());
			MITEM r = MEvaluate(i.Right());
			a = MMultiply(l, r);
		}
		break;
	case MDIV:
		{
			MITEM l = MEvaluate(i.Left());
			MITEM r = MEvaluate(i.Right());
			a = MDivide(l, r);
		}
		break;
	case MPOW:
		{
			MITEM l = MEvaluate(i.Left());
			MITEM r = MEvaluate(i.Right());
			a = (l ^ r);
		}
		break;
	case MEQUATION:
		{
			MITEM l = MEvaluate(i.Left());
			MITEM r = MEvaluate(i.Right());
			a = MITEM(new MEquation(l.copy(), r.copy()));
		}
		break;
	case MEQUALITY:
		{
			MITEM l = MEvaluate(i.Left());
			MITEM r = MEvaluate(i.Right());
			bool b = is_equal(l.ItemPtr(), r.ItemPtr());
			a = MITEM(new MBoolean(b));
		}
		break;
	case MMATRIX:
		{
			const MMatrix& A = *mmatrix(i);
			a = MEvaluate(A);
		}
		break;
	case MF1D:
		{
			const MFunc1D* pf = mfnc1d(i);
			string s = pf->Name();
			MITEM p = MEvaluate(i.Item());
			if (s.compare("sin" ) == 0) return Sin  (p);
			if (s.compare("cos" ) == 0) return Cos  (p);
			if (s.compare("tan" ) == 0) return Tan  (p);
			if (s.compare("cot" ) == 0) return Cot  (p);
			if (s.compare("sec" ) == 0) return Sec  (p);
			if (s.compare("csc" ) == 0) return Csc  (p);
			if (s.compare("atan") == 0) return Atan (p);
			if (s.compare("cosh") == 0) return Cosh (p);
			if (s.compare("sinh") == 0) return Sinh (p);
			if (s.compare("exp" ) == 0) return Exp  (p);
			if (s.compare("ln"  ) == 0) return Log  (p);
			if (s.compare("log" ) == 0) return Log10(p);
			if (s.compare("sqrt") == 0) return Sqrt (p);
			if (s.compare("fl"  ) == 0) return Float(p);
			if (s.compare("abs" ) == 0) return Abs  (p);
			if (s.compare("sgn" ) == 0) return Sgn  (p);
#ifdef WIN32
			if (s.compare("J0"  ) == 0) return J0   (p);
			if (s.compare("J1"  ) == 0) return J1   (p);
			if (s.compare("Y0"  ) == 0) return Y0   (p);
			if (s.compare("Y1"  ) == 0) return Y1   (p);
#endif
			if (s.compare("fac" ) == 0) return Fac  (p);
			if (s.compare("erf" ) == 0) return Erf  (p);
			if (s.compare("erfc") == 0) return Erfc (p);
			return i;
		}
		break;
	case MF2D:
		{
			const MFunc2D* pf = mfnc2d(i);
			string s = pf->Name();
			MITEM a = MEvaluate(i.Left());
			MITEM b = MEvaluate(i.Right());
			if (s.compare("binomial") == 0) return Binomial(a, b);
#ifdef WIN32
			else if (s.compare("Jn")==0) return Jn(a,b);
			else if (s.compare("Yn")==0) return Yn(a,b);
#endif
			else if (s.compare("Tn")==0) return Tn(a,b);
			return i;
		}
		break;
	case MFMAT:
		{
			const MFuncMat* pfm = mfncmat(i);
			MITEM m = MEvaluate(i.Item());
			const MMatrix& A = *mmatrix(m.ItemPtr());
			FUNCMATPTR pf = pfm->funcptr();
			a = pf(A);
		}
		break;
	case MFMAT2:
		{
			const MFuncMat2* pfm = mfncmat2(i);
			MITEM l = MEvaluate(i.Left());
			MITEM r = MEvaluate(i.Right());
			const MMatrix& A = *mmatrix(l.ItemPtr());
			const MMatrix& B = *mmatrix(r.ItemPtr());
			FUNCMAT2PTR pf = pfm->funcptr();
			a = pf(A, B);
		}
		break;
	case MSFNC:
		{
			const MSFuncND* pf = msfncnd(i);
			MITEM v = pf->Value()->copy();
			MITEM e = MEvaluate(v);
//			if (is_rconst(e.ItemPtr())) return e; else return i;
			return e;
		}
		break;
	default:
		return i;
	}

	if (a.Type() != i.Type()) a = MEvaluate(a);
	return a;
}

//-----------------------------------------------------------------------------
MITEM MEvaluate(const MMatrix& A)
{
	MMatrix& B = *(new MMatrix());
	B.Create(A.rows(), A.columns());
	for (int i=0; i<A.rows(); ++i)
		for (int j=0; j<A.columns(); ++j)
		{
			MITEM aij = A[i][j]->copy();
			B[i][j] = MEvaluate(aij).copy();
		}
	return (&B);
}

//-----------------------------------------------------------------------------
MITEM MMultiply(const MITEM& l, const MITEM& r)
{
	MProduct M(l*r);
	return M.Item();
}

//-----------------------------------------------------------------------------
MITEM MDivide(const MITEM& n, const MITEM& d)
{
	if (d == 0.0) throw DivisionByZero();
	if (n == d) return 1.0;
	MProduct L(n), D(d);
	return L / D;
}

//-----------------------------------------------------------------------------
MProduct::MProduct(const MITEM& a)
{
	m_p.push_back(1.0);
	Multiply(a);
}

//-----------------------------------------------------------------------------
void MProduct::Multiply(const MITEM& a)
{
	if (is_mul(a))
	{
		MITEM l = a.Left();
		MITEM r = a.Right();
		Multiply(l);
		Multiply(r);
	}
	else 
	{
		if (isConst(a))
		{
			MITEM& i0 = *m_p.begin();
			if (isConst(i0)) i0 = (a.value()*i0.value()); else m_p.push_front(a);
		}
		else 
		{
			// see if a already appears in the product
			bool binsert = true;
			MITEM b(a), pb(1.0);
			if (is_pow(a))
			{
				b = a.Left();
				pb = a.Right();
			}
			list<MITEM>::iterator ic;
			for (ic=m_p.begin(); ic != m_p.end(); ++ic)
			{
				MITEM c(*ic), pc(1.0);
				if (is_pow(c))
				{
					c = (*ic).Left();
					pc = (*ic).Right();
				}

				if (b == c)
				{
					(*ic) = c^(pc + pb);
					binsert = false;
					break;
				}
			}

			if (binsert) m_p.push_back(a);
		}
	}
}

//-----------------------------------------------------------------------------
MITEM MProduct::Item()
{
	MITEM c(1.0), r(1.0);
	list<MITEM>::iterator it;
	for (it=m_p.begin(); it != m_p.end(); ++it) 
	{
		if (isConst(*it)) c = c*(*it);
		else r = (r*(*it));
	}
	return c*r;
}

//-----------------------------------------------------------------------------
MITEM MProduct::operator / (const MProduct& d)
{
	MProduct tmp(d);
	list<MITEM>::iterator ib;
	list<MITEM>::iterator ia;
	for (ib = tmp.m_p.begin(); ib != tmp.m_p.end();)
	{
		MITEM b(*ib), pb(1.0);
		bool binc = true;
		if (b != 0.0)
		{
			if (is_pow(b))
			{
				b  = (*ib).Left();
				pb = (*ib).Right();
			}
			for (ia = m_p.begin(); ia != m_p.end(); ++ia)
			{
				MITEM a(*ia), pa(1.0);
				if (is_pow(a))
				{
					a  = (*ia).Left();
					pa = (*ia).Right();
				}

				if (a == b)
				{
					(*ia) = a^(pa - pb);
					ib = tmp.m_p.erase(ib);
					binc = false;
					break;
				}
			}
		}
		if (binc) ++ib;
	}

	if (!tmp.m_p.empty())
	{
		if (isConst(*m_p.begin()) && isConst(*tmp.m_p.begin()))
		{
			MITEM& a = *m_p.begin();
			MITEM& b = *tmp.m_p.begin();
			double n = a.value();
			double d = b.value();
			if ((d != 0) && is_int(n) && is_int(d))
			{
				long  c = gcf((long) n, (long) d);
				n /= (double) c;
				d /= (double) c;
				a = n;
				b = d;
			}
		}
	}

	return Item() / tmp.Item();
}

//-----------------------------------------------------------------------------
bool MProduct::operator==(const MProduct& b)
{
	list<MITEM>::iterator pf;
	for (pf = m_p.begin(); pf != m_p.end(); ++pf)
	{
		MITEM& i = *pf;
		if (b.contains(i) == false) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool MProduct::contains(const MITEM& a) const
{
	list<MITEM>::const_iterator pf;
	for (pf = m_p.begin(); pf != m_p.end(); ++pf)
	{
		const MITEM& i = *pf;
		if (i == a) return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
MITEM MAddition(const MITEM& l, const MITEM& r)
{
	MSum S(l+r);
	return S.Item();
}

//-----------------------------------------------------------------------------
MSum::MSum(const MITEM& a)
{
	m_c = 0.0;
	Add(a);
}

//-----------------------------------------------------------------------------
MSum::MTerm::MTerm(const MITEM& i) : m_s(1.0)
{
	// extract the scalar from i
	if (is_mul(i))
	{
		MITEM l = i.Left();
		if (isConst(l))
		{
			m_s = l.value();
			m_a = i.Right();
		}
		else m_a = i;
	}
	else if (is_div(i))
	{
		MITEM n = i.Left();
		MITEM d = i.Right();

		if (is_mul(n))
		{
			MITEM l = n.Left();
			if (isConst(l))
			{
				m_s = l.value();
				n = n.Right();
			}
		}

		if (is_mul(d))
		{
			MITEM l = d.Left();
			if (isConst(l))
			{
				m_s /= l.value();
				d = d.Right();
			}
		}
		else if (isConst(d))
		{
			m_s /= d.value();
			d = 1.0;
		}

		m_a = MDivide(n, d);
	}
	else m_a = i;
}

//-----------------------------------------------------------------------------
void MSum::Add(const MITEM& a)
{
	if (is_add(a))
	{
		Add(a.Left());
		Add(a.Right());
	}
	else if (is_sub(a))
	{
		Add(a.Left());
		Sub(a.Right());
	}
	else if (isConst(a)) m_c += a.value();
	else if (is_frac(a)) m_c += mfrac(a)->fraction();
	else if (m_t.empty()) m_t.push_back(MTerm(a));
	else
	{
		MTerm b(a);
		bool badd = true;
		// see if term already exists
		list<MTerm>::iterator pt;
		for (pt = m_t.begin(); pt != m_t.end(); ++pt)
		{
			if (b.m_a == pt->m_a)
			{
				badd = false;
				pt->m_s += b.m_s;
				break;
			}
		}
		if (badd) m_t.push_back(b);
	}
}

//-----------------------------------------------------------------------------
void MSum::Sub(const MITEM& a)
{
	if (is_add(a))
	{
		Sub(a.Left());
		Sub(a.Right());
	}
	else if (is_sub(a))
	{
		Sub(a.Left());
		Add(a.Right());
	}
	else if (isConst(a)) m_c -= a.value();
	else if (is_frac(a)) m_c -= mfrac(a)->fraction();
	else
	{
		MTerm b(a);

		bool badd = true;
		// see if term already exists
		list<MTerm>::iterator pt;
		for (pt = m_t.begin(); pt != m_t.end(); ++pt)
		{
			if (b.m_a == pt->m_a)
			{
				badd = false;
				pt->m_s -= b.m_s;
				break;
			}
		}
		if (badd) { b.m_s = -b.m_s; m_t.push_back(b); }
	}
}

//-----------------------------------------------------------------------------
MITEM MSum::Item()
{
	if (m_t.empty()) return Fraction(m_c.n, m_c.d);

	list<MTerm>::iterator pt = m_t.begin();
	MITEM g = Fraction(pt->m_s.n, pt->m_s.d);
	MITEM s = g*pt->m_a;
	++pt;
	for (; pt != m_t.end(); ++pt)
	{
		g = Fraction(pt->m_s.n, pt->m_s.d);
		s = s + g*pt->m_a;
	}
	if (m_c.n != 0.0) s = s + Fraction(m_c.n, m_c.d);
	return s;
}
