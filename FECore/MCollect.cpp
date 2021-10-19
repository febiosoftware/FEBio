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

//---------------------------------------------------
void MCollect(const MITEM& e, const MITEM& x, MITEM& a, MITEM& b);

//---------------------------------------------------
// collect terms in x
MITEM MCollect(const MITEM& e, const MITEM& x)
{
	MITEM a(0.0);
	MITEM b(0.0);
	MCollect(e, x, a, b);
	return MEvaluate(a*x + b);
}

//---------------------------------------------------
void MCollect(const MITEM& e, const MITEM& x, MITEM& a, MITEM& b)
{
	// check for equality
	if (e == x)
	{
		a = a + 1.0;
		return;
	}

	// check for dependancy
	if (is_dependent(e, x) == false)
	{
		b = MEvaluate(b + e);
		return;
	}

	// process operators
	switch (e.Type())
	{
	case MNEG:
		{
			MITEM c(0.0), d(0.0);
			MCollect(e.Item(), x, c, d);
			a = a - c;
			b = b - d;
		}
		break;
	case MADD:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_dependent(l, x) == false)
			{
				b = MEvaluate(b + l);
				MCollect(r, x, a, b);
			}
			else if (is_dependent(r, x) == false)
			{
				b = MEvaluate(b + r);
				MCollect(l, x, a, b);
			}
			else
			{
				MITEM al(0.0), ar(0.0), bl(0.0), br(0.0);
				MCollect(l, x, al, bl);
				MCollect(r, x, ar, br);
				a = a + al + ar;
				b = b + bl + br;
			}
		}
		break;
	case MSUB:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_dependent(l, x) == false)
			{
				b = MEvaluate(b + l);
				MCollect(-r, x, a, b);
			}
			else if (is_dependent(r, x) == false)
			{
				b = MEvaluate(b - r);
				MCollect(l, x, a, b);
			}
			else
			{
				MITEM al(0.0), ar(0.0), bl(0.0), br(0.0);
				MCollect(l, x, al, bl);
				MCollect(r, x, ar, br);
				a = a + al - ar;
				b = b + bl - br;
			}
		}
		break;
	case MMUL:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (l == x) a = a + r;
			else if (r == x) a = a + l;
			else if (is_dependent(l, x) == false)
			{
				MCollect(r, x, a, b);
				a = MEvaluate(l*a);
				b = MEvaluate(l*b);
			}
			else if (is_dependent(r, x) == false)
			{
				MCollect(l, x, a, b);
				a = MEvaluate(r*a);
				b = MEvaluate(r*b);
			}
			else b = b + e;
		}
		break;
	case MDIV:
		{
			MITEM l = e.Left();
			MITEM r = e.Right();
			if (is_dependent(l, x))
			{
				MCollect(l, x, a, b);
				a = MEvaluate(a/r);
				b = MEvaluate(b/r);
			}
			else b = b + e;
		}
		break;
	default:
		b = b + e;
	}
}
