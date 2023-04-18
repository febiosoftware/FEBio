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
#include "MItem.h"
#include "MMatrix.h"
#include "MEvaluate.h"
#include "MMath.h"
using namespace std;

#ifndef PI
#define PI	3.141592653589793
#endif

//-----------------------------------------------------------------------------
// Some special constants
MITEM MPi (new MNamedCt(PI      , "pi"));
MITEM ME  (new MNamedCt(exp(1.0), "e" ));
MITEM MINF(new MNamedCt(1e308   , "_inf_"));

//-----------------------------------------------------------------------------
MItem* MFuncND::copy() const
{
	return new MFuncND(m_pf, m_name, GetSequence());
}

//-----------------------------------------------------------------------------
MSequence::MSequence(const MSequence& s) : MItem(MSEQUENCE)
{
	int N = s.size();
	for (int i=0; i<N; ++i) add(s[i]->copy());
}

//-----------------------------------------------------------------------------
MSequence::~MSequence()
{
	int M = size();
	for (int i=0; i<M; ++i) delete m_item[i]; 
	m_item.clear();
}

//-----------------------------------------------------------------------------
void MSequence::remove(int i)
{
	delete m_item[i];
	m_item.erase(m_item.begin()+i);
}

//-----------------------------------------------------------------------------
void MSequence::replace(int i, MItem* pi)
{
	delete m_item[i];
	m_item[i] = pi;
}

//-----------------------------------------------------------------------------
MSequence& MSequence::operator = (MSequence& s)
{
	int M = size();
	for (int i=0; i<M; ++i) delete m_item[i]; m_item.clear();
	int N = s.size();
	for (int i=0; i<N; ++i) add(s[i]->copy());
	return (*this);
}


//-----------------------------------------------------------------------------
MITEM Fraction(double n, double d)
{
	double s = (n*d<0?-1:1);
	n = fabs(n);
	d = fabs(d);

	if ((n - floor(n) == 0) && (d - floor(d) == 0))
	{
		if (d != 0.0) 
		{
			if (n==0) return 0.0;
			double f = n/d;
			if (f - floor(f) == 0)
			{
				MITEM c(f);
				if (s<0) return -c; else return c;
			}

			double c = (double) gcf((long) n, (long) d);
			n /= c;
			d /= c;
		}
	}
/*	else if (d - floor(d) == 0.0)
	{
		double v = n / d;
		return v;
	}
*/
	if (d == 1.0) return (s*n);
	MITEM r = new MFraction(n, d);
	if (s<0) return -r; else return r;
}

//=============================================================================
//                        N E G A T I O N   ( - )
//=============================================================================

//-----------------------------------------------------------------------------
MITEM operator - (const MITEM& l) 
{
	// get the Item pointer
	const MItem* pi = l.ItemPtr();

	// remove unnecessary negative signs
	int n = 0;
	while (is_neg(pi)) { pi = mneg(pi)->Item(); n = (n+1)%2; }

	// -(a-b) = b-a
	if (is_sub(pi))
	{
		const MItem* pl = msub(pi)->LeftItem ();
		const MItem* pr = msub(pi)->RightItem();

		if (!is_neg(pl) && !is_neg(pr))
		{
			MITEM a = pl->copy();
			MITEM b = pr->copy();
			return (b - a);
		}
	}

	// remove negative from zero
	if (is_zero(pi)) return 0.0;

	// return the negative of a copy of pi
	MItem* pret = pi->copy();
	if (n==0) pret = new MNeg(pret);
	return pret; 
}

//-----------------------------------------------------------------------------


//=============================================================================
//                            A D D I T I O N   ( + )
//=============================================================================

//-----------------------------------------------------------------------------
MITEM operator + (const MITEM& l, const MITEM& r)
{
	if (is_matrix(l) || is_matrix(r))
	{
		if (is_matrix(l) && is_matrix(r))
		{
			const MMatrix& A = *mmatrix(l);
			const MMatrix& B = *mmatrix(r);
			return A+B;
		}
//		else throw InvalidOperation();
	}
	if (isConst(l)) return (r + l.value());
	if (isConst(r)) return (l + r.value());
	if (is_neg  (r)) return (l - (-r)); // a+(-b) = a - b
	if (is_neg  (l)) return (r - (-l)); // ((-a)+b) = b - a
	if (is_frac(l) && is_frac(r)) 
	{
		FRACTION a = mfrac(l)->fraction();
		FRACTION b = mfrac(r)->fraction();
		return Fraction(a + b);
	}
	return new MAdd(l.copy(), r.copy());
}

//-----------------------------------------------------------------------------
MITEM operator + (const MITEM& l, double r)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (r==0) return l;
	if (isConst(l)) return (l.value() + r);
	if (is_neg  (l)) return (r - l.Item());
	if (is_frac (l))
	{
		FRACTION a = mfrac(l)->fraction();
		return Fraction(a + r);
	}
	return new MAdd(l.copy(), r);
}

//-----------------------------------------------------------------------------
MITEM operator + (double l, const MITEM& r)
{
	if (is_matrix(r)) throw InvalidOperation();
	if (l==0) return r;
	if (isConst(r)) return (l + r.value());
	if (is_neg  (r)) return (l - r.Item());
	if (is_frac (r))
	{
		FRACTION a = mfrac(r)->fraction();
		return Fraction(a + l);
	}
	return new MAdd(r.copy(), l);
}

//=============================================================================
//                         S U B T R A C T I O N   ( - )
//=============================================================================

//-----------------------------------------------------------------------------
MITEM operator - (const MITEM& l, const MITEM& r)
{
	if (is_matrix(l) || is_matrix(r))
	{
		if (is_matrix(l) && is_matrix(r))
		{
			const MMatrix& A = *mmatrix(l);
			const MMatrix& B = *mmatrix(r);
			return A-B;
		}
//		else throw InvalidOperation();
	}
	if (isConst(l)) return (l.value() - r);
	if (isConst(r)) return (l - r.value());
	if (is_neg(l) && is_neg(r)) return (-r)-(-l);  // (-a)-(-b) = b - a;
	if (is_neg(r)) return (l+(-r));
	if (is_neg(l)) return -((-l)+r);
	if (is_frac(l) && is_frac(r)) 
	{
		FRACTION a = mfrac(l)->fraction();
		FRACTION b = mfrac(r)->fraction();
		return Fraction(a - b);
	}
	return new MSub(l.copy(), r.copy());
}

//-----------------------------------------------------------------------------
MITEM operator - (double l, const MITEM& r)
{
	if (is_matrix(r)) throw InvalidOperation();
	if (l==0) return -r;
	if (isConst(r)) return (l - r.value());
	if (is_neg  (r)) return (-r) + l;
	if (is_frac (r))
	{
		FRACTION a = mfrac(r)->fraction();
		return Fraction(l - a);
	}
	return new MSub(l, r.copy());
}

//-----------------------------------------------------------------------------
MITEM operator - (const MITEM& l, double r)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (r == 0) return l;
	if (isConst(l)) return (l.value() - r);
	if (is_neg  (l)) return -((-l) + r);
	if (is_frac (l))
	{
		FRACTION a = mfrac(l)->fraction();
		return Fraction(a - r);
	}
	else return new MSub(l.copy(), r);
}

//=============================================================================
//                      M U L T I P L I C A T I O N   ( * )
//=============================================================================

//-----------------------------------------------------------------------------
MITEM operator * (const MITEM& l, const MITEM& r)
{
	if (is_matrix(l) && ::is_scalar(r))
	{
		const MMatrix& A = *mmatrix(l);
		return A*r;
	}
	if (is_matrix(r) && ::is_scalar(l))
	{
		const MMatrix& A = *mmatrix(r);
		return A*l;
	}
	if (isConst(l)) return (l.value() * r);
	if (isConst(r)) return (r.value() * l);
	if (is_neg(l) && is_neg(r)) return ((-l)*(-r)); 
	if (is_neg(l)) return -((-l)*r);
	if (is_neg(r)) return -(l*(-r));
	if (is_frac(l) && is_frac(r))
	{
		FRACTION a = mfrac(l)->fraction();
		FRACTION b = mfrac(r)->fraction();
		return Fraction(a*b);
	}
	if (is_frac(l))
	{
		FRACTION a = mfrac(l)->fraction();
		return (a.n*r)/a.d;
	}
	if (is_frac(r))
	{
		FRACTION a = mfrac(r)->fraction();
		return (a.n*l)/a.d;
	}
	if (is_div(l) && (is_matrix(r) == false))
	{
		MITEM l1 = l.Left();
		MITEM l2 = l.Right();
		return (l1*r)/l2;
	}
	if (is_div(r) && (is_matrix(r) == false))
	{
		MITEM r1 = r.Left();
		MITEM r2 = r.Right();
		return (l*r1)/r2;
	}
	if (is_matrix(l) && is_matrix(r))
	{
		const MMatrix& A = *mmatrix(l);
		const MMatrix& B = *mmatrix(r);
		return A*B;
	}
	return new MMul(l.copy(), r.copy());
}

//-----------------------------------------------------------------------------
MITEM operator * (double l, const MITEM& r)
{ 
	if (l == 0.0) return 0.0;
	if (l == 1.0) return r;
	if (isConst(r)) return l*r.value();
	if (is_neg  (r)) return -(l * (-r)); 
	if (is_frac (r))
	{
		FRACTION a = mfrac(r)->fraction();
		return Fraction(l*a);
	}
	if (is_div(r))
	{
		MITEM fl = r.Left ();
		MITEM fr = r.Right();
		return (l*fl)/fr;
	}
	return new MMul(l, r.copy()); 
}

//=============================================================================
//                            D I V I S I O N   ( / )
//=============================================================================

MITEM operator / (const MITEM& l, const MITEM& r)
{
	if (r == 0.0) throw DivisionByZero();
	if (r == 1.0) return l;
	if (is_matrix(r)) throw InvalidOperation();
	if (is_matrix(l))
	{
		const MMatrix& A = *mmatrix(l);
		return A/r;
	}
	if (isConst(l)) return (l.value() / r);
	if (isConst(r)) return (l / r.value());
	if (is_neg  (l)) return -((-l)/r);
	if (is_neg  (r)) return -(l/(-r));
	if (is_frac(l) && is_frac(r))
	{
		FRACTION a = mfrac(l)->fraction();
		FRACTION b = mfrac(r)->fraction();
		return Fraction(a/b);
	}
	if (is_div(r))
	{
		MITEM r1 = r.Left();
		MITEM r2 = r.Right();
		return (l*r2)/r1;
	}
	if (is_div(l))
	{
		MITEM l1 = l.Left();
		MITEM l2 = l.Right();
		return (l1/(l2*r));
	}
	//------------------------------------------------------
	return new MDiv(l.copy(), r.copy());
}

//-----------------------------------------------------------------------------
MITEM operator / (double l, const MITEM& r)
{
	if (is_matrix(r)) throw InvalidOperation();
	if ((l==0.0) && (r!=0.0)) return 0.0;
	if (isConst(r)) return Fraction(l, r.value());
	if (is_neg  (r)) return -(l/(-r));
	if (is_frac (r))
	{
		FRACTION a = mfrac(r)->fraction();
		return Fraction(l/a);
	}
	if (is_div(r))
	{
		MITEM a = r.Left();
		MITEM b = r.Right();
		return (l*b)/a;
	}
	return new MDiv(l, r.copy()); 
}

//-----------------------------------------------------------------------------
MITEM operator / (const MITEM& l, double r)
{
	if (r == 1.0) return l;
	if (isConst(l)) return Fraction(l.value(), r);
	if (is_neg  (l)) return -((-l)/r);
	if (is_frac (l))
	{
		FRACTION a = mfrac(l)->fraction();
		return Fraction(a/r);
	}
	if (is_div(l))
	{
		MITEM a = l.Left();
		MITEM b = l.Right();
		return (a/(r*b));
	}
	return new MDiv(l.copy(), r); 
}

//=============================================================================
//                            P O W E R   ( ^ )
//=============================================================================

//-----------------------------------------------------------------------------
MITEM operator ^ (const MITEM& l, const MITEM& r)
{
	if (is_matrix(l) || is_matrix(r)) throw InvalidOperation();
	if (isConst(r)) return (l ^ r.value());
	if (is_neg(r))
	{
		MITEM mr = -r;
		if (isConst(mr)) return 1/ (l^mr);
	}
	if (is_fnc1d(l, "sqrt"))
	{
		MITEM a = l.Item();
		return (a^(r / 2));
	}
	return new MPow(l.copy(), r.copy());
}

//-----------------------------------------------------------------------------
MITEM operator ^ (const MITEM& l, double r)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (r==1) return l;
	if (r==0) return 1.0;
	if (l==1.0) return 1.0;
	if (isConst(l) && is_int(r)) return pow(l.value(), r);
	if (is_neg(l) && is_int(r))
	{
		int n = (int) r;
		if (n%2==0) return (-l)^r;
		else return -((-l)^r);
	}
	if (is_frac(l) && is_int(r))
	{
		FRACTION a = mfrac(l)->fraction();
		return Fraction(pow(a.n,r), pow(a.d,r));
	}
	if (is_pow(l))
	{
		MITEM v = l.Left();
		MITEM p = l.Right();
		return (v^(r*p));
	}
	if (is_mul(l))
	{
		MITEM l1 = l.Left ();
		MITEM r1 = l.Right();
		return ((l1^r)*(r1^r));
	}
	if (is_div(l))
	{
		MITEM l1 = l.Left ();
		MITEM r1 = l.Right();
		return ((l1^r)/(r1^r));
	}
	if (is_fnc1d(l, "sqrt"))
	{
		MITEM a = l.Item();
		return (a^(Fraction(r, 2)));
	}
	return new MPow(l.copy(), r);
}

//=============================================================================
//                            F U N C T I O N S
//=============================================================================

//-----------------------------------------------------------------------------
// Absolute value
MITEM Abs(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	if (is_neg(a)) return Abs(-a);

	// remove redundant abs
	const MItem* pi = a.ItemPtr();
	while (is_neg(pi) || is_abs(pi)) { pi = munary(pi)->Item(); }

	if (isConst(pi)) return fabs(mnumber(pi)->value());
	if (is_named(pi)) return pi->copy(); // This assumes that all named constants are positive!
	return new MFunc1D(fabs, "abs", pi->copy());
}

//-----------------------------------------------------------------------------
// return the sign of the expression
MITEM Sgn(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	if (isConst(a)) return 1.0; // this assumes that all constants are positive
	if (is_neg(a)) return -Sgn(-a);
	if (is_pow(a))
	{
		const MItem *pr = mpow(a)->RightItem();
		if (isConst(pr) && is_int(mnumber(pr)->value()))
		{
			MITEM l = a.Left();
			int n = (int) mnumber(pr)->value();
			if (n%2 == 0) return 1.0;
			else return Sgn(l);
		}
	}
	return new MFunc1D(sgn, "sgn", a.copy());
}

//-----------------------------------------------------------------------------
// This function sees if a is a multiple of pi where the multiplier is a fraction
// of integers. 
bool pi_multiple(const MITEM& a, int& num, int& den)
{
	if (a==0.0) { num = 0; den = 1; return true; }
	if (is_pi(a)) { num = den = 1; return true; }
	else if (is_mul(a))
	{
		const MItem* pl = mmul(a)->LeftItem();
		const MItem* pr = mmul(a)->RightItem();
		if (is_pi(pl) && is_int(pr)) { num = (int)mconst(pr)->value(); den = 1; return true; }
		if (is_pi(pr) && is_int(pl)) { num = (int)mconst(pl)->value(); den = 1; return true; }
	}
	else if (is_div(a))
	{
		const MItem* pl = mdiv(a)->LeftItem();
		const MItem* pr = mdiv(a)->RightItem();
		if (is_pi(pl)&&is_int(pr)) { num = 1; den = (int)mconst(pr)->value(); return true; }
		if (is_mul(pl)&&is_int(pr))
		{
			den = (int)mconst(pr)->value();
			const MItem* pl1 = mmul(pl)->LeftItem();
			const MItem* pl2 = mmul(pl)->RightItem();
			if (is_pi(pl1) && (is_int(pl2))) { num = (int)mconst(pl2)->value(); return true; }
			if (is_pi(pl2) && (is_int(pl1))) { num = (int)mconst(pl1)->value(); return true; }
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
MITEM Sin(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	// sin(-x) = -sin(x)
	if (is_neg(a)) return -Sin(-a);

	// check for multiples of pi
	int num, den;
	if (pi_multiple(a, num, den))
	{
		if (num == 0) return 0.0;
		if (den == 1) return 0.0;
		if ((num == 1) && (den == 2)) return 1.0;
		if ((num == 1) && (den == 3)) return Sqrt(MITEM(3.0))/2.0;
		if ((num == 1) && (den == 4)) return Sqrt(MITEM(2.0))/2.0;
		if ((num == 1) && (den == 6)) return 1.0 / MITEM(2.0);
		if ((num == 1) && (den == 12)) return Sqrt(MITEM(2.0))/4*(Sqrt(MITEM(3.0)) - 1.0);
		if ((num == 5) && (den == 12)) return Sqrt(MITEM(2.0))/4*(Sqrt(MITEM(3.0)) + 1.0);
		if ((num == 7) && (den == 12)) return Sqrt(MITEM(2.0))/4*(Sqrt(MITEM(3.0)) + 1.0);
		if ((num == 2) && (den == 3 )) return Sqrt(MITEM(3.0))/2;
		if ((num == 3) && (den == 4 )) return Sqrt(MITEM(2.0))/2;
		if ((num == 5) && (den == 6 )) return 1.0 / MITEM(2.0);
		if ((num == 11) && (den == 12)) return Sqrt(MITEM(2.0))/4*(Sqrt(MITEM(3.0)) - 1.0);
	}
	return new MFunc1D(sin, "sin", a.copy());
}

//-----------------------------------------------------------------------------
MITEM Cos(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	// cos(-x) = cos(x)
	if (is_neg(a)) return Cos(-a);

	// check for multiples of pi
	int num, den;
	if (pi_multiple(a, num, den))
	{
		if (num == 0) return  1.0;
		if (den == 1) return (num%2==0?1.0:-1.0);
		if ((num == 1) && (den == 2)) return 0.0;
		if ((num == 1) && (den == 3)) return 1.0 / MITEM(2.0);
		if ((num == 1) && (den == 4)) return Sqrt(MITEM(2.0))/2.0;
		if ((num == 1) && (den == 6)) return Sqrt(MITEM(3.0))/2.0;
		if ((num == 1) && (den == 12)) return Sqrt(MITEM(2.0))/4*(Sqrt(MITEM(3.0)) + 1.0);
		if ((num == 5) && (den == 12)) return Sqrt(MITEM(2.0))/4*(Sqrt(MITEM(3.0)) - 1.0);
		if ((num == 7) && (den == 12)) return -Sqrt(MITEM(2.0))/4*(Sqrt(MITEM(3.0)) - 1.0);
		if ((num == 2) && (den == 3 )) return -1.0 / MITEM(2.0);
		if ((num == 3) && (den == 4 )) return -Sqrt(MITEM(2.0))/2;
		if ((num == 5) && (den == 6 )) return -Sqrt(MITEM(3.0))/2;
		if ((num == 11) && (den == 12)) return -Sqrt(MITEM(2.0))/4*(Sqrt(MITEM(3.0)) + 1.0);
	}
	return new MFunc1D(cos, "cos", a.copy());
}

//-----------------------------------------------------------------------------
MITEM Sec(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();

	// sec(-x) = sec(x)
	if (is_neg(a)) return Sec(-a);

	// check for multiples of pi
	int num, den;
	if (pi_multiple(a, num, den))
	{
		if (num == 0) return  1.0;
		if (den == 1) return (num%2==0?1.0:-1.0);
//		if ((num == 1) && (den == 2)) return INFINITY;
		if ((num == 1) && (den == 3)) return 2.0;
		if ((num == 1) && (den == 4)) return Sqrt(MITEM(2.0));
		if ((num == 1) && (den == 6)) return 2.0*Sqrt(MITEM(3.0))/3.0;
		if ((num == 1) && (den == 12)) return Sqrt(MITEM(2.0))*(Sqrt(MITEM(3.0)) - 1.0);
		if ((num == 5) && (den == 12)) return Sqrt(MITEM(2.0))*(Sqrt(MITEM(3.0)) + 1.0);
		if ((num == 7) && (den == 12)) return -Sqrt(MITEM(2.0))*(Sqrt(MITEM(3.0)) + 1.0);
		if ((num == 2) && (den == 3 )) return -2.0;
		if ((num == 3) && (den == 4 )) return -Sqrt(MITEM(2.0));
		if ((num == 5) && (den == 6 )) return -(2.0*Sqrt(MITEM(3.0))/3);
		if ((num == 11) && (den == 12)) return -Sqrt(MITEM(2.0))*(Sqrt(MITEM(3.0)) - 1.0);
	}

	return new MFunc1D(sec, "sec", a.copy());
}

//-----------------------------------------------------------------------------
MITEM Csc(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	// csc(-x) = -csc(x)
	if (is_neg(a)) return -Csc(-a);

	// check for multiples of pi
	int num, den;
	if (pi_multiple(a, num, den))
	{
//		if (num == 0) return INFINITY;
//		if (den == 1) return INFINITY;
		if ((num == 1) && (den == 2)) return 1.0;
		if ((num == 1) && (den == 3)) return 2.0*Sqrt(MITEM(3.0))/3.0;
		if ((num == 1) && (den == 4)) return Sqrt(MITEM(2.0));
		if ((num == 1) && (den == 6)) return 2.0;
		if ((num == 1) && (den == 12)) return Sqrt(MITEM(2.0))*(Sqrt(MITEM(3.0)) + 1.0);
		if ((num == 5) && (den == 12)) return Sqrt(MITEM(2.0))*(Sqrt(MITEM(3.0)) - 1.0);
		if ((num == 7) && (den == 12)) return Sqrt(MITEM(2.0))*(Sqrt(MITEM(3.0)) - 1.0);
		if ((num == 2) && (den == 3 )) return 2.0*Sqrt(MITEM(3.0))/3.0;
		if ((num == 3) && (den == 4 )) return Sqrt(MITEM(2.0));
		if ((num == 5) && (den == 6 )) return 2.0;
		if ((num == 11) && (den == 12)) return Sqrt(MITEM(2.0))*(Sqrt(MITEM(3.0)) + 1.0);
	}

	return new MFunc1D(csc, "csc", a.copy());
}

//-----------------------------------------------------------------------------
MITEM Tan(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();

	// tan(-x) = -tan(x)
	if (is_neg(a)) return -Tan(-a);

	// check for multiples of pi
	int num, den;
	if (pi_multiple(a, num, den))
	{
		if (num == 0) return 0.0;
		if (den == 1) return 0.0;
//		if ((num == 1) && (den == 2)) +/- return INFINITY;
		if ((num == 1) && (den == 3)) return Sqrt(MITEM(3.0));
		if ((num == 1) && (den == 4)) return 1.0;
		if ((num == 1) && (den == 6)) return Sqrt(MITEM(3.0))/3.0;
		if ((num == 1) && (den == 12)) return MITEM(2.0) - Sqrt(MITEM(3.0));
		if ((num == 5) && (den == 12)) return MITEM(2.0) + Sqrt(MITEM(3.0));
		if ((num == 7) && (den == 12)) return -(MITEM(2.0) + Sqrt(MITEM(3.0)));
		if ((num == 2) && (den == 3 )) return -Sqrt(MITEM(3.0));
		if ((num == 3) && (den == 4 )) return -1.0;
		if ((num == 5) && (den == 6 )) return -Sqrt(MITEM(3.0))/3.0;
		if ((num == 11) && (den == 12)) return -(MITEM(2.0) - Sqrt(MITEM(3.0)));
	}

	return new MFunc1D(tan, "tan", a.copy());
}

//-----------------------------------------------------------------------------
MITEM Cot(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();

	// cot(-x) = -cot(x)
	if (is_neg(a)) return -Cot(-a);

	// check for multiples of pi
	int num, den;
	if (pi_multiple(a, num, den))
	{
//		if (num == 0) return +/- INFINITY;
//		if (den == 1) return +/- INFINITY;
		if ((num == 1) && (den == 2)) return 0.0;
		if ((num == 1) && (den == 3)) return Sqrt(MITEM(3.0))/3.0;
		if ((num == 1) && (den == 4)) return 1.0;
		if ((num == 1) && (den == 6)) return Sqrt(MITEM(3.0));
		if ((num == 1) && (den == 12)) return MITEM(2.0) + Sqrt(MITEM(3.0));
		if ((num == 5) && (den == 12)) return MITEM(2.0) - Sqrt(MITEM(3.0));
		if ((num == 7) && (den == 12)) return -(MITEM(2.0) - Sqrt(MITEM(3.0)));
		if ((num == 2) && (den == 3 )) return -Sqrt(MITEM(3.0))/3.0;
		if ((num == 3) && (den == 4 )) return -1.0;
		if ((num == 5) && (den == 6 )) return -Sqrt(MITEM(3.0));
		if ((num == 11) && (den == 12)) return -(MITEM(2.0) + Sqrt(MITEM(3.0)));
	}
	return new MFunc1D(cot, "cot", a.copy());
}

//-----------------------------------------------------------------------------
MITEM Atan(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	if (is_neg(a)) return -Atan(-a);
	if (a == 0.0) return 0.0;
	if (a == 1.0) return MPi/4.0;
	return new MFunc1D(atan, "atan", a.copy());
}

//-----------------------------------------------------------------------------
MITEM Cosh(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (l == 0) return 1.0;
	if (is_neg(l)) return Cosh(-l);
	return new MFunc1D(cosh, "cosh", l.copy());
}

//-----------------------------------------------------------------------------
MITEM Sinh(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (l == 0) return 0.0;
	if (is_neg(l)) return -Sinh(-l);
	return new MFunc1D(sinh, "sinh", l.copy());
}

//-----------------------------------------------------------------------------
MITEM Exp(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	return ME^a;
}

//-----------------------------------------------------------------------------
// Natural logarithm (base e)
MITEM Log(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	if (a==0.0) return -MINF;
	if (a==1.0) return 0.0;
	if (is_named(a))
	{
		string sz = mnamed(a)->Name();
		if (sz.compare("e") == 0) return 1.0;
	}
	if (is_pow(a))
	{
		const MItem* pl = mpow(a)->LeftItem();
		const MItem* pr = mpow(a)->RightItem();
		if (is_named(pl))
		{
			string sz = mnamed(pl)->Name();
			if (sz.compare("e") == 0)
			{
				return pr->copy();
			}
		}
	}
	return new MFunc1D(log, "ln", a.copy());
}

//-----------------------------------------------------------------------------
// base 10 log
MITEM Log10(const MITEM& a)
{
	if (is_matrix(a)) throw InvalidOperation();
	if (isConst(a))
	{
		double w = a.value();
		if (w > 0)
		{
			double p = log10(w);
			double pi = (double) ((int) p);
			if (pi == p) return p; 
		}
	}
	return new MFunc1D(log10, "log", a.copy());
}

//-----------------------------------------------------------------------------
MITEM Sqrt(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (l == 1.0) return 1.0;
	if (is_pow(l))
	{
		MITEM a = l.Left();
		MITEM b = l.Right();
		return a^(b/2);
	}
	if (isConst(l))
	{
		double a = l.value();
		double f = sqrt(a);
		if (is_int(f)) return f;

		return new MFunc1D(sqrt, "sqrt", l.copy());
	}
	if (is_var(l) || is_named(l))
	{
		return new MFunc1D(sqrt, "sqrt", l.copy());
	}
	return l ^ Fraction(1,2);
}

//-----------------------------------------------------------------------------
MITEM Erf(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (l == 0.0) return 0.0;
	return new MFunc1D(erf, "erf", l.copy());
}

//-----------------------------------------------------------------------------
MITEM Erfc(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (l == 0.0) return 1.0;
	return new MFunc1D(erfc, "erfc", l.copy());
}

//-----------------------------------------------------------------------------
MITEM Tn(int n, const MITEM& r)
{
	if (is_matrix(r)) throw InvalidOperation();
	if (n >= 0)
	{
		if (n == 0) return MITEM(1.0);
		if (n == 1) return r;
		return MEvaluate(MExpand(2*r*Tn(n-1,r) - Tn(n-2,r)));
	}
	else return new MFunc2D(chebyshev, "Tn", new MConstant((double)n), r.copy());
}

//-----------------------------------------------------------------------------
MITEM Tn(const MITEM& l, const MITEM& r)
{
	if (is_matrix(l) || is_matrix(r)) throw InvalidOperation();
	if (is_int(l) && (l.value() >= 0.0)) return Tn((int) l.value(), r);
	return new MFunc2D(chebyshev, "Tn", l.copy(), r.copy());
}

//-----------------------------------------------------------------------------
#ifdef WIN32
MITEM J0(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (l == 0.0) return 1.0;
	return new MFunc1D(_j0, "J0", l.copy());
}
#endif

//-----------------------------------------------------------------------------
#ifdef WIN32
MITEM J1(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (l == 0.0) return 0.0;
	return new MFunc1D(_j1, "J1", l.copy());
}
#endif

//-----------------------------------------------------------------------------
#ifdef WIN32
MITEM Jn(int n, const MITEM& r)
{
	if (is_matrix(r)) throw InvalidOperation();
	if (n >= 0)
	{
		if (r == 0.0) return ((n == 0)?1.0:0.0);
	}
	return new MFunc2D(jn, "Jn", new MConstant((double)n), r.copy());
}
#endif

//-----------------------------------------------------------------------------
#ifdef WIN32
MITEM Jn(const MITEM& l, const MITEM& r)
{
	if (is_matrix(l) || is_matrix(r)) throw InvalidOperation();
	if (is_int(l)) return Jn((int)l.value(), r);
	return new MFunc2D(jn, "Jn", l.copy(), r.copy());
}
#endif

//-----------------------------------------------------------------------------
#ifdef WIN32
MITEM Y0(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	return new MFunc1D(_y0, "Y0", l.copy());
}
#endif

//-----------------------------------------------------------------------------
#ifdef WIN32
MITEM Y1(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	return new MFunc1D(_y1, "Y1", l.copy());
}
#endif

//-----------------------------------------------------------------------------
#ifdef WIN32
MITEM Yn(int n, const MITEM& r)
{
	if (is_matrix(r)) throw InvalidOperation();
	return new MFunc2D(yn, "Yn", new MConstant((double)n), r.copy());
}
#endif

//-----------------------------------------------------------------------------
#ifdef WIN32
MITEM Yn(const MITEM& l, const MITEM& r)
{
	if (is_matrix(l) || is_matrix(r)) throw InvalidOperation();
	if (is_int(l)) return Yn((int)l.value(), r);
	return new MFunc2D(yn, "Yn", l.copy(), r.copy());
}
#endif

//-----------------------------------------------------------------------------
MITEM Fac(const MITEM& l)
{
	if (is_matrix(l)) throw InvalidOperation();
	if (is_int(l) && (l.value() >= 0)) return fac(l.value());
	return new MFunc1D(fac, "fac", l.copy());
}

//-----------------------------------------------------------------------------
MITEM Binomial(const MITEM& l, const MITEM& r)
{
	if (is_matrix(l) || is_matrix(r)) throw InvalidOperation();
	if (is_rconst(l.ItemPtr()) && is_rconst(r.ItemPtr())) return binomial(l.value(), r.value());
	return new MFunc2D(binomial, "binomial", l.copy(), r.copy());
}

//-----------------------------------------------------------------------------
MITEM Identity(const MITEM& n)
{
	if (is_int(n))
	{
		int m = (int) n.value();
		if (m <= 0) throw InvalidOperation();
		return matrix_identity(m);
	}
	else throw InvalidOperation();
}

//-----------------------------------------------------------------------------
MITEM Float(const MITEM& a)
{
	if (is_rconst(a.ItemPtr()))
	{
		if (is_number(a.ItemPtr())) return a.value();
		if (is_neg(a)) return -Float(-a);
		if (is_add(a)) return Float(a.Left()).value() + Float(a.Right()).value();
		if (is_sub(a)) return Float(a.Left()).value() - Float(a.Right()).value();
		if (is_mul(a)) return Float(a.Left()).value() * Float(a.Right()).value();
		if (is_div(a)) return Float(a.Left()).value() / Float(a.Right()).value();
		if (is_pow(a))
		{
			if (is_neg(mpow(a)->RightItem())) return pow(Float(a.Left()).value(), -(Float(-(a.Right())).value()));
			else return pow(Float(a.Left()).value(), Float(a.Right()).value());
		}
		if (is_func1d(a))
		{
			FUNCPTR pf = mfnc1d(a)->funcptr();
			MITEM v(mfnc1d(a)->Item()->copy());
			return pf(Float(v).value());
		}
		return a.value();
	}
	return a;
}

//-----------------------------------------------------------------------------
bool is_equal(const MItem* pl, const MItem* pr)
{
	if ((pl == 0)&&(pr == 0)) return true;
	if ((pl == 0)||(pr == 0)) return false;

	Item_Type t = pl->Type();
	if (t != pr->Type()) return false;
	switch (t)
	{
	case MCONST: 
	case MFRAC:
	case MNAMED: return (mnumber(pl)->value() == mnumber(pr)->value());
	case MVAR:
		{
			const string& vl = mvar(pl)->Name();
			const string& vr = mvar(pr)->Name();
			return (vl.compare(vr) == 0);
		}
	case MNEG:
		{
			const MItem* pa = munary(pl)->Item();
			const MItem* pb = munary(pr)->Item();
			return is_equal(pa, pb);
		}
	case MADD:
		{
			const MItem* pl1 = mbinary(pl)->LeftItem(), *pr1 = mbinary(pl)->RightItem();
			const MItem* pl2 = mbinary(pr)->LeftItem(), *pr2 = mbinary(pr)->RightItem();

			// TODO: since this only looks one level deep, this may not always return the correct answer
			return ((is_equal(pl1, pl2) && is_equal(pr1, pr2)) || (is_equal(pl1, pr2) && is_equal(pr1, pl2)));
		}
	case MMUL:
		{
			MITEM l(pl->copy());
			MITEM r(pr->copy());
			MProduct a(l), b(r);
			return (a==b);
		}
	case MSUB:
	case MDIV:
	case MPOW:
		{
			const MItem* pl1 = mbinary(pl)->LeftItem(), *pr1 = mbinary(pl)->RightItem();
			const MItem* pl2 = mbinary(pr)->LeftItem(), *pr2 = mbinary(pr)->RightItem();
			return (is_equal(pl1, pl2) && is_equal(pr1, pr2));
		}
	case MF1D:
		{
			if (mfnc1d(pl)->funcptr() == mfnc1d(pr)->funcptr())
			{
				const MItem* pa = munary(pl)->Item();
				const MItem* pb = munary(pr)->Item();
				return is_equal(pa, pb);
			}
			else return false;
		}
	case MSFNC:
		{
			const string& vl = msfncnd(pl)->Name();
			const string& vr = msfncnd(pr)->Name();
			if (vl.compare(vr) != 0) return false;
			const MItem* pa = msfncnd(pl)->Value();
			const MItem* pb = msfncnd(pr)->Value();
			return is_equal(pa, pb);
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
// See if the expression is constant. That is, there are no variables in
// this expression.
bool is_rconst(const MItem* pi)
{
	if (is_var(pi)) return false;
	if (isConst(pi) || is_named(pi) || is_frac(pi)) return true;
	if (is_unary(pi)) return is_rconst(munary(pi)->Item());
	if (is_binary(pi))
	{
		const MItem* pl = mbinary(pi)->LeftItem ();
		const MItem* pr = mbinary(pi)->RightItem();
		return (is_rconst(pl) && is_rconst(pr));
	}
	if (is_nary(pi))
	{
		const MNary* pn = mnary(pi);
		int n = pn->Params();
		for (int i=0; i<n; ++i)
			if (is_rconst(pn->Param(i)) == false) return false;
		return true;
	}
	if (is_matrix(pi))
	{
		const MMatrix& A = *mmatrix(pi);
		for (int i=0; i<A.rows(); ++i)
			for (int j=0; j<A.columns(); ++j)
			{
				MItem* pij = A[i][j];
				if (is_rconst(pij) == false) return false;
			}
		return true;
	}
	assert(false);
	return false;
}

//-----------------------------------------------------------------------------
// See if the expression is dependent on the variable x
bool is_dependent(const MItem* pi, const MVariable& x)
{
	if (is_var(pi)) return (mvar(pi)->Name().compare(x.Name()) == 0);
	if (isConst(pi) || is_named(pi) || is_frac(pi)) return false;
	if (is_unary(pi)) return is_dependent(munary(pi)->Item(), x);
	if (is_binary(pi))
	{
		const MItem* pl = mbinary(pi)->LeftItem ();
		const MItem* pr = mbinary(pi)->RightItem();
		return (is_dependent(pl, x) || is_dependent(pr, x));
	}
	if (is_nary(pi))
	{
		const MNary* pn = mnary(pi);
		int n = pn->Params();
		bool b = false;
		for (int i=0; i<n; ++i) b = (b || is_dependent(pn->Param(i), x));
		return b;
	}
	if (is_matrix(pi))
	{
		const MMatrix& m = *mmatrix(pi);
		int nr = m.rows();
		int nc = m.columns();
		for (int i=0; i<nr; ++i)
			for (int j=0; j<nc; ++j)
			{
				const MItem* mij = m(i,j);
				if (is_dependent(mij, x)) return true;
			}
		return false;
	}
	assert(false);
	return false;
}

//-----------------------------------------------------------------------------
// See if the expression contains the expression px
bool is_dependent(const MItem* pi, const MItem* px)
{
	if (is_equal(pi, px)) return true;
	if (is_unary(pi)) return is_dependent(munary(pi)->Item(), px);
	if (is_binary(pi))
	{
		const MItem* pl = mbinary(pi)->LeftItem ();
		const MItem* pr = mbinary(pi)->RightItem();
		return (is_dependent(pl, px) || is_dependent(pr, px));
	}
	if (is_nary(pi))
	{
		const MNary* pn = mnary(pi);
		int n = pn->Params();
		bool b = false;
		for (int i=0; i<n; ++i) b = (b || is_dependent(pn->Param(i), px));
		return b;
	}
	if (is_matrix(pi))
	{
		const MMatrix& m = *mmatrix(pi);
		int nr = m.rows();
		int nc = m.columns();
		for (int i=0; i<nr; ++i)
			for (int j=0; j<nc; ++j)
			{
				const MItem* mij = m(i,j);
				if (is_dependent(mij, px)) return true;
			}
		return false;
	}
	return false;
}

//-----------------------------------------------------------------------------
// see if the item is the named constant pi
bool is_pi(const MItem* pi)
{
	if (is_named(pi))
	{
		const string& sz = mnamed(pi)->Name();
		if (sz.compare("pi") == 0) return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// see if the expression is a scalar expression (e.g. does not contain matrices)
bool is_scalar(const MItem* pi)
{
	if (is_number(pi)) return true;
	if (is_unary(pi)) return ::is_scalar(munary(pi)->Item());
	if (is_binary(pi))
	{
		const MItem* pl = mbinary(pi)->LeftItem ();
		const MItem* pr = mbinary(pi)->RightItem();
		return (::is_scalar(pl) && ::is_scalar(pr));
	}
	if (is_nary(pi))
	{
		const MNary* pn = mnary(pi);
		int n = pn->Params();
		bool b = true;
		for (int i=0; i<n; ++i) b = (b && ::is_scalar(pn->Param(i)));
		return b;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Calculate the number of operations for the expression.
int op_count(const MItem* pi)
{
	int n = 0;
	if (is_unary(pi)) n = op_count(munary(pi)->Item()) + 1;
	else if (is_binary(pi))
	{
		int n1 = op_count(mbinary(pi)->LeftItem());
		int n2 = op_count(mbinary(pi)->RightItem());
		n = n1 + n2 + 1;
	}
	else if (is_nary(pi))
	{
		n = 1;
		const MNary* pn = mnary(pi);
		int N = pn->Params();
		for (int i=0; i<N; ++i) n += op_count(pn->Param(i));
	}
	return n;
}

//-----------------------------------------------------------------------------
const char* read_format(const MItem* pe, const char* sz)
{
	switch (sz[0])
	{
	case '+':
		{
			const MBinary* pb = mbinary(pe);
			if (pe->Type() != MADD) return nullptr;
			sz = read_format(pb->LeftItem(), sz+1);
			read_format(pb->RightItem(), sz);			
		}
		break;
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool is_format(const MItem* pe, const char* sz)
{
	return (read_format(pe, sz) != 0);
}
