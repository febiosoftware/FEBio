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
#include "MathObject.h"
#include "MMath.h"
#include "MObjBuilder.h"
#include <float.h>
using namespace std;

//-----------------------------------------------------------------------------
MathObject::MathObject()
{

}

MathObject::MathObject(const MathObject& mo)
{
	// copy variable list
	for (int i = 0; i < mo.Variables(); ++i)
	{
		AddVariable(new MVariable(*mo.Variable(i)));
	}
}

void MathObject::Clear()
{
	for (int i = 0; i < m_Var.size(); ++i) delete m_Var[i];
	m_Var.clear();
}

//-----------------------------------------------------------------------------
void MathObject::operator = (const MathObject& mo)
{
	// copy variable list
	for (int i = 0; i < mo.Variables(); ++i)
	{
		AddVariable(new MVariable(*mo.Variable(i)));
	}
}

//-----------------------------------------------------------------------------
MathObject::~MathObject()
{
	for (int i = 0; i < m_Var.size(); ++i) delete m_Var[i];
	m_Var.clear();
}

//-----------------------------------------------------------------------------
MVariable* MathObject::AddVariable(const std::string& var, double initVal)
{
	if (m_Var.empty() == false)
	{
		MVarList::iterator it = m_Var.begin();
		for (it = m_Var.end(); it != m_Var.end(); ++it)
			if ((*it)->Name() == var) return *it;
	}

	MVariable* pv = new MVariable(var, initVal);
	pv->setIndex((int)m_Var.size());
	m_Var.push_back(pv);
	return pv;
}

void MathObject::AddVariable(MVariable* pv)
{
	if (m_Var.empty() == false)
	{
		MVarList::iterator it = m_Var.begin();
		for (it = m_Var.end(); it != m_Var.end(); ++it)
			if (*it == pv) return;
	}
	pv->setIndex((int)m_Var.size());
	m_Var.push_back(pv);
}

void MathObject::AddVariables(const vector<std::string>& varList)
{
	for (const std::string& s : varList)
	{
		AddVariable(s);
	}
}

//-----------------------------------------------------------------------------
MVariable* MathObject::FindVariable(const string& s)
{
	for (MVarList::iterator it = m_Var.begin(); it != m_Var.end(); ++it)
		if ((*it)->Name() == s) return *it;
	return 0;
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
double MSimpleExpression::value(const std::string& s)
{
	Create(s);
	return value();
}

//-----------------------------------------------------------------------------
double MSimpleExpression::value(const MItem* pi) const
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
double MSimpleExpression::value(const MItem* pi, const std::vector<double>& var) const
{
	switch (pi->Type())
	{
	case MCONST: return (mnumber(pi)->value());
	case MFRAC : return (mnumber(pi)->value());
	case MNAMED: return (mnumber(pi)->value());
	case MVAR  : return (var[mvar(pi)->index()]);
	case MNEG: return -value(munary(pi)->Item(), var);
	case MADD: return value(mbinary(pi)->LeftItem(), var) + value(mbinary(pi)->RightItem(), var);
	case MSUB: return value(mbinary(pi)->LeftItem(), var) - value(mbinary(pi)->RightItem(), var);
	case MMUL: return value(mbinary(pi)->LeftItem(), var) * value(mbinary(pi)->RightItem(), var);
	case MDIV: return value(mbinary(pi)->LeftItem(), var) / value(mbinary(pi)->RightItem(), var);
	case MPOW: return pow(value(mbinary(pi)->LeftItem(), var), value(mbinary(pi)->RightItem(), var));
	case MF1D:
		{
			double a = value(munary(pi)->Item(), var);
			return (mfnc1d(pi)->funcptr())(a);
		}
		break;
	case MF2D:
		{
			double a = value(mbinary(pi)->LeftItem(), var);
			double b = value(mbinary(pi)->RightItem(), var);
			return (mfnc2d(pi)->funcptr())(a, b);
		}
		break;
	case MSFNC:
		{
			return value(msfncnd(pi)->Value(), var);
		};		
	default:
		assert(false);
		return 0;
	}
}

//-----------------------------------------------------------------------------
MSimpleExpression::MSimpleExpression(const MSimpleExpression& mo) : MathObject(mo), m_item(mo.m_item)
{
	// The copy c'tor of MathObject copied the variables, but any MVarRefs still point to the mo object, not this object's var list.
	// Calling the following function fixes this
	fixVariableRefs(m_item.ItemPtr());
}

//-----------------------------------------------------------------------------
void MSimpleExpression::operator=(const MSimpleExpression& mo)
{
	// copy base object
	MathObject::operator=(mo);

	// copy the item
	m_item = mo.m_item;

	// The = operator of MathObject copied the variables, but any MVarRefs still point to the mo object, not this object's var list.
	// Calling the following function fixes this
	fixVariableRefs(m_item.ItemPtr());
}

//-----------------------------------------------------------------------------
void MSimpleExpression::fixVariableRefs(MItem* pi)
{
	if (pi->Type() == MVAR)
	{
		MVarRef* var = static_cast<MVarRef*>(pi);
		int index = var->GetVariable()->index();
		var->SetVariable(Variable(index));
	}
	else if (is_unary(pi))
	{
		MUnary* uno = static_cast<MUnary*>(pi);
		fixVariableRefs(uno->Item());
	}
	else if (is_binary(pi))
	{
		MBinary* bin = static_cast<MBinary*>(pi);
		fixVariableRefs(bin->LeftItem());
		fixVariableRefs(bin->RightItem());
	}
	else if (is_nary(pi))
	{
		MNary* any = static_cast<MNary*>(pi);
		for (int i = 0; i < any->Params(); ++i) fixVariableRefs(any->Param(i));
	}
}

//-----------------------------------------------------------------------------
// Create a simple expression object from a string
bool MSimpleExpression::Create(const std::string& expr, bool autoVars)
{
	MObjBuilder mob;
	mob.setAutoVars(autoVars);
	if (mob.Create(this, expr, false) == false) return false;
	return true;
}
