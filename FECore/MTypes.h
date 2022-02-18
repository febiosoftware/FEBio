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

//-----------------------------------------------------------------------------
// class that describes a fraction of two real numbers
class FRACTION
{
public:
	double n, d;		// nominator, denominator

public:
	FRACTION() : n(0), d(1) {}
	FRACTION(double a) : n(a), d(1) {}
	FRACTION(double a, double b) : n(a), d(b) {}
	FRACTION(const FRACTION& f) { n = f.n; d = f.d; }

	FRACTION& operator = (double a) { n = a; d = 1; return (*this); }
	FRACTION& operator = (const FRACTION& f) { n = f.n; d = f.d; return (*this); }

	FRACTION& operator += (double a) { n += a*d; return (*this); }
	FRACTION& operator += (const FRACTION& a) { n = n*a.d + d*a.n; d *= a.d; return (*this); }

	FRACTION& operator -= (double a) { n -= a*d; return (*this); }
	FRACTION& operator -= (const FRACTION& a) { n = n*a.d - d*a.n; d *= a.d; return (*this); }

	FRACTION& operator /= (double a) { d *= a; return (*this); }

	operator double() { return n / d; }

	void normalize();
};

//-----------------------------------------------------------------------------
// operators for FRACTION
inline FRACTION operator + (FRACTION& l, FRACTION& r) { return FRACTION(l.n*r.d + r.n*l.d, l.d*r.d); }
inline FRACTION operator + (FRACTION& l, double    r) { return FRACTION(l.n + r*l.d, l.d); }
inline FRACTION operator + (double    l, FRACTION& r) { return FRACTION(l*r.d + r.n, r.d); }

inline FRACTION operator - (FRACTION& l, FRACTION& r) { return FRACTION(l.n*r.d - r.n*l.d, l.d*r.d); }
inline FRACTION operator - (FRACTION& l, double    r) { return FRACTION(l.n - r*l.d, l.d); }
inline FRACTION operator - (double    l, FRACTION& r) { return FRACTION(l*r.d - r.n, r.d); }

inline FRACTION operator * (FRACTION& l, FRACTION& r) { return FRACTION(l.n*r.n, l.d*r.d); }
inline FRACTION operator * (FRACTION& l, double    r) { return FRACTION(l.n*r, l.d); }
inline FRACTION operator * (double    l, FRACTION& r) { return FRACTION(r.n*l, r.d); }

inline FRACTION operator / (FRACTION& l, FRACTION& r) { return FRACTION(l.n*r.d, l.d*r.n); }
inline FRACTION operator / (FRACTION& l, double    r) { return FRACTION(l.n, l.d*r); }
inline FRACTION operator / (double    l, FRACTION& r) { return FRACTION(l*r.d, r.n); }

inline FRACTION operator - (FRACTION& a) { return FRACTION(-a.n, a.d); }
//-----------------------------------------------------------------------------
// Calculates the greatest common factor
long gcf(long a, long b);

