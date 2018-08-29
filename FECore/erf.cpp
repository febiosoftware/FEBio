#include "stdafx.h"
#include "erf.h"
#include <math.h>

#ifdef WIN32
//-----------------------------------------------------------------------------
// approximation to error function using rational expansion.
double erf(double x)
{
	const double p = 0.3275911;
	const double a1 =  0.254829592;
	const double a2 = -0.284496736;
	const double a3 =  1.421413741;
	const double a4 = -1.453152027;
	const double a5 =  1.061405429;
	double t = 1.0/(1.0 + p*fabs(x));
	double e = t*(a1 + t*(a2 + t*(a3 + t*(a4 + a5*t))))*exp(-x*x);
	return (x > 0 ? 1.0 - e : e - 1.0);
}

//-----------------------------------------------------------------------------
// complementary error function (erfc(x) = 1 - erf(x))
double erfc(double x)
{
	return 1.0 - erf(x);
}
#endif
