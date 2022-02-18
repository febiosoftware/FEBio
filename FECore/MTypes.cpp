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
