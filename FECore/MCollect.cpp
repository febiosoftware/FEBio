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
