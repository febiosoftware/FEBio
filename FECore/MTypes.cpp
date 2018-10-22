#include "stdafx.h"
#include "MTypes.h"
#include <math.h>

//-----------------------------------------------------------------------------
// Calculates the greatest common factor of two integers
long gcf(long a, long b)
{
	if ((a==0) || (b==0)) return 1;
	if (a<0) a = -a;
	if (b<0) b = -b;

	if (a > b) { a ^= b; b ^= a; a ^= b; }

	long q, r;
	do
	{
		q = b/a;
		r = b%a;

		if (r > 0)
		{
			b = (b-r)/q;
			a = r;
		}
	}
	while (r>0);

	return a;
}

//-----------------------------------------------------------------------------
// reduces the fraction to it simplest form
void FRACTION::normalize()
{
	double s = (n*d<0?-1:1);
	n = fabs(n);
	d = fabs(d);

	if (d != 0)
	{
		long in = (long ) n;
		long id = (long ) d;
		if ((n==(double) in) && (d == (double) id))
		{
			double c = (double) gcf(in, id);
			n /= c;
			d /= c;
		}
	}

	if (s<0) n = -n;
}
