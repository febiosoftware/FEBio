#include "MFunctions.h"
#include <assert.h>

double fl  (double x) { return x; }
double csc (double x) { return 1.0 / sin(x); }
double sec (double x) { return 1.0 / cos(x); }
double cot (double x) { return 1.0 / tan(x); }
double sinc(double x) { return (x == 0? 1 : sin(x)/x); }
double sgn (double x) { return (x<0?-1.0:1.0); }

#ifdef WIN32
double jn(double x, double y)
{
	return _jn((int) x, y);
}

double yn(double x, double y)
{
	return _yn((int) x, y);
}
#endif

double fmax(double* x, int n)
{
	double v = x[0];
	for (int i=1; i<n; ++i)
		if (x[i] > v) v = x[i];
	return v;
}

double fmin(double* x, int n)
{
	double v = x[0];
	for (int i=1; i<n; ++i)
		if (x[i] < v) v = x[i];
	return v;
}

double avg(double* x, int n)
{
	double v = x[0];
	for (int i=1; i<n; ++i) v += x[i];
	return (v/n);
}

double chebyshev(double f, double x)
{
	int n = (int)(f);
	if (n<=0) return 1;
	if (n==1) return x;

	double T0 = 0;
	double T1 = 1;
	double Tn = x;
	for (int i=2; i<=n; ++i)
	{
		T0 = T1;
		T1 = Tn;
		Tn = 2*x*T1 - T0;
	}
	return Tn;
}

double fac(double n)
{
	assert(n >= 0.0);
	if (n <= 1) return 1.0;
	double p = 1.0;
	while (n > 1.0) p *= (n--);
	return p;
}

double prod(double a, double b)
{
	double p = 1;
	if (a >= b)
	{
		while (a >= b) p *= a--;
	}
	else
	{
		while (b >= a) p *= b--;
	}
	return p;
}

double binomial(double n, double r)
{
	assert(r <= n);
	if (n == r) return 1.0;
	if (r == 0.0) return 1.0;
	double b = 0;
	if (r >= n/2)
	{
		b = prod(n, r+1)/fac(n - r);
	}
	else b = prod(n, n-r+1) / fac(r);
	return b;
}

//-----------------------------------------------------------------------------
// approximation of gamma function using Lanczos approximation
double gamma(double z)
{
	const int g = 7;
	const double p[] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
     771.32342877765313, -176.61502916214059, 12.507343278686905,
	 -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
	const double pi = 4.0*atan(1.0);
 
    if (z < 0.5)
	{
        return pi / (sin(pi*z) * gamma(1-z));
	}
    else
	{
        z -= 1;
        double x = p[0];
        for (int i=1; i<g+2; ++i) x += p[i]/(z+i);
		double t = z + g + 0.5;
        return sqrt(2*pi) * pow(t, z+0.5) * exp(-t) * x;
	}
}
