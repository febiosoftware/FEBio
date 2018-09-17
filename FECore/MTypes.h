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

