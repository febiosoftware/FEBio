#include "MMath.h"
#include "MEvaluate.h"

//-----------------------------------------------------------------------------
MITEM MSolve(MITEM& e, MVariable& x);
MITEM MSolve(MITEM& l, MITEM& r, MVariable& x);
MITEM MSolve(MSequence& e, MSequence& x);
MITEM MSolve(MEquation& e, MVariable& x);

//-----------------------------------------------------------------------------
// Solve an expression e for the variable(s) v.
// The expression can be an equation linear in v or a system of linear equations.
MITEM MSolve(MITEM& e, MITEM& v)
{
	if (is_sequence(e))
	{
		// solve a linear system of equations
		// make sure v contains a list of variables
		if (is_sequence(v) == false) throw InvalidOperation();

		MSequence& se = *msequence(e);
		MSequence& sv = *msequence(v);

		// make sure we have as many equations as unknowns
		if (se.size() != sv.size()) throw InvalidOperation();

		// make sure all entries in e are equations
		for (int i=0; i<se.size(); ++i) if (is_equation(se[i]) == false) throw InvalidOperation();

		// make sure all items in v are variables
		for (int i=0; i<sv.size(); ++i) if (is_var(sv[i]) == false) throw InvalidOperation();

		// solve the system of equations
		return MSolve(se, sv);
	}
	else
	{
		if (is_var(v))
		{
			MVariable& x = *mvar(v)->GetVariable();
			return MSolve(e, x);
		}
		else throw InvalidOperation();
	}
	assert(false);
	return e;
}

//-----------------------------------------------------------------------------
// solve a system of linear equations
MITEM MSolve(MSequence& e, MSequence& v)
{
	MSequence sol = e;
	MSequence var = v;
	int neq = e.size();
	for (int i=0; i<neq; ++i)
	{
		MEquation& eqi = *mequation(sol[i]);

		// find a variable this equation depends on
		int nvar = var.size();
		MVariable* pv = 0;
		for (int j=0; j<nvar; ++j)
		{
			MVariable* pvj = mvar(var[j])->GetVariable();
			if (is_dependent(sol[i], *pvj))
			{ 
				pv = pvj; 
				var.remove(j);
				break; 
			}
		}
		if (pv == 0) throw InvalidOperation();

		MVariable& x = *pv;

		// solve the i-th equation for x
		MITEM si = MSolve(eqi, x);

		// make sure the LHS is the variable
		if (is_equation(si) == false) throw InvalidOperation();
		MEquation& ei = *mequation(si);
		if (is_equal(ei.LeftItem(), x) == false) throw InvalidOperation();

		// replace i-th equation with this solution
		sol.replace(i, si.copy());

		// substitute the solution in all other equations
		for (int j=0; j<neq; ++j)
		{
			if (i != j)
			{
				MITEM ej = sol[j]->copy();
				MITEM newj = MReplace(ej, si);
				sol.replace(j, newj.copy());
			}
		}
	}
	return sol.copy();
}

//-----------------------------------------------------------------------------
MITEM MSolve(MEquation& e, MVariable& x)
{
	MITEM l = e.LeftItem()->copy();
	MITEM r = e.RightItem()->copy();
	return MSolve(l, r, x);
}

//-----------------------------------------------------------------------------
MITEM MSolve(MITEM& e, MVariable& x)
{
	if (is_equation(e))
	{
		MEquation& eq = *mequation(e);
		return MSolve(eq, x);
	}
	else
	{
		MITEM i(new MEquation(e.copy(), new MConstant(0)));
		return MSolve(i, x);
	}
}

//-----------------------------------------------------------------------------
MITEM MSolve(MITEM& L, MITEM& r, MVariable& x)
{
	// if the RHS depends on x, bring it to the left side
	if (is_dependent(r,x))
	{
		MITEM c(0.0);
		return MEvaluate(MSolve(L - r, c, x));
	}

	// simplify LHS as much as possible
	MITEM l = MEvaluate(L);

	// check LHS for dependancy on x and try to solve
	switch (l.Type())
	{
	case MNEG: return MEvaluate(MSolve(l.Item(), -r, x)); break;
	case MADD:
		{
			MITEM a = l.Left();
			MITEM b = l.Right();
			if (is_dependent(a, x) == false) return MEvaluate(MSolve(b, r - a, x));
			if (is_dependent(b, x) == false) return MEvaluate(MSolve(a, r - b, x));
			else
			{
				MITEM v(x);
				MITEM a = MCollect(l, v);
				return MSolve(a, r, x);
			}
		}
		break;
	case MSUB:
		{
			MITEM a = l.Left();
			MITEM b = l.Right();
			if (is_dependent(a, x) == false) return MEvaluate(MSolve(-b, r - a, x));
			if (is_dependent(b, x) == false) return MEvaluate(MSolve(a, r + b, x));
			else
			{
				MITEM v(x);
				MITEM a = MCollect(l, v);
				return MSolve(a, r, x);
			}
		}
		break;
	case MMUL:
		{
			MITEM a = l.Left();
			MITEM b = l.Right();
			if (is_dependent(a, x) == false) return MEvaluate(MSolve(b, r/a, x));
			if (is_dependent(b, x) == false) return MEvaluate(MSolve(a, r/b, x));
		}
		break;
	case MDIV:
		{
			MITEM a = l.Left();
			MITEM b = l.Right();
			return MEvaluate(MSolve(a, r*b, x));
		}
		break;
	}

	return MITEM(new MEquation(l.copy(), r.copy()));
}
