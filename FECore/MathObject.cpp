#include "MathObject.h"
#include "MMath.h"
#include <float.h>
using namespace std;

//-----------------------------------------------------------------------------
bool operator > (const string& a, const string& b) { return (strcmp(a.c_str(), b.c_str()) > 0); }

void MathObject::AddVariable(MVariable* pv)
{
	if (m_Var.empty() == false)
	{
		MVarList::iterator it = m_Var.begin();
		for (it = m_Var.end(); it != m_Var.end(); ++it)
			if (*it == pv) return;
	}
	m_Var.push_back(pv);
}

//-----------------------------------------------------------------------------
MVariable* MathObject::FindVariable(const string& s)
{
	for (MVarList::iterator it = m_Var.begin(); it != m_Var.end(); ++it)
		if ((*it)->Name() == s) return *it;
	return 0;
}

//-----------------------------------------------------------------------------
// check an item for variables recursively
void MathObject::BuildVarList(MItem* pi)
{
	if (is_var(pi))
	{
		// see if we already have this variable
		MVarRef& var = *mvar(pi);
		MVariable* pv = FindVariable(var.Name());
		if (pv != 0) return;

		// if not, add it to the list
		m_Var.push_back(var.GetVariable());
	}
	else if (is_unary(pi))
	{
		MUnary* pe = munary(pi);
		BuildVarList(pe->Item());
	}
	else if (is_binary(pi))
	{
		MBinary* pe = mbinary(pi);
		BuildVarList(pe->LeftItem());
		BuildVarList(pe->RightItem());
	}
	else if (is_nary(pi))
	{
		MNary* pe = mnary(pi);
		for (int i=0; i<pe->Params(); ++i) BuildVarList(pe->Param(i));
	}
}

//-----------------------------------------------------------------------------
int MSimpleExpression::Items()
{
	if (m_item.Type() == MSEQUENCE)
	{
		MSequence* pe = dynamic_cast<MSequence*>(m_item.ItemPtr());
		return pe->size();
	}
	else return 1;
}

//-----------------------------------------------------------------------------
double MSimpleExpression::value(MItem* pi)
{
	switch (pi->Type())
	{
	case MCONST: return (mnumber(pi)->value());
	case MFRAC : return (mnumber(pi)->value());
	case MNAMED: return (mnumber(pi)->value());
	case MVAR  : return (mnumber(pi)->value());
	case MNEG: return -value(munary(pi)->Item());
	case MADD: return value(mbinary(pi)->LeftItem()) + value(mbinary(pi)->RightItem());
	case MSUB: return value(mbinary(pi)->LeftItem()) - value(mbinary(pi)->RightItem());
	case MMUL: return value(mbinary(pi)->LeftItem()) * value(mbinary(pi)->RightItem());
	case MDIV: return value(mbinary(pi)->LeftItem()) / value(mbinary(pi)->RightItem());
	case MPOW: return pow(value(mbinary(pi)->LeftItem()), value(mbinary(pi)->RightItem()));
	case MF1D:
		{
			double a = value(munary(pi)->Item());
			return (mfnc1d(pi)->funcptr())(a);
		}
		break;
	case MF2D:
		{
			double a = value(mbinary(pi)->LeftItem());
			double b = value(mbinary(pi)->RightItem());
			return (mfnc2d(pi)->funcptr())(a, b);
		}
		break;
	case MSFNC:
		{
			return value(msfncnd(pi)->Value());
		};		
	default:
		assert(false);
		return 0;
	}
}

//-----------------------------------------------------------------------------
MSFuncND* MFuncDef::CreateFunction(MSequence& s)
{
	assert((int)s.size() == GetVariables());
	MSequence v;
	int N = GetVariables();
	for (int i=0; i<N; ++i)
	{
		MITEM var(new MVarRef(m_pv[i]));
		v.add(new MVarRef(m_pv[i]));
	}
	MITEM res = MReplace(m_item, v, s);
	return new MSFuncND(m_name, res.copy(), &s);
}
