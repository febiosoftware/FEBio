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
#include "MObj2String.h"
#include <sstream>
using namespace std;

//-----------------------------------------------------------------------------
// Convert a math object to a string
string MObj2String::Convert(const MathObject& o)
{
	const MathObject* po = &o;
	if (dynamic_cast<const MSimpleExpression*>(po))
	{
		const MSimpleExpression* pe = dynamic_cast<const MSimpleExpression*>(po);
		const MITEM& i = pe->GetExpression();
		return Convert(i.ItemPtr());
	}
	assert(false);
	return "";
}

//-----------------------------------------------------------------------------
// Convert a math item to string
string MObj2String::Convert(const MItem* pi)
{
	switch(pi->Type())
	{
	case MCONST   : return Constant (dynamic_cast<const MConstant*> (pi));
	case MFRAC    : return Fraction (dynamic_cast<const MFraction*> (pi));
	case MNAMED   : return NamedCt  (dynamic_cast<const MNamedCt* > (pi));
	case MVAR     : return Variable (dynamic_cast<const MVarRef*  > (pi));
	case MNEG     : return OpNeg    (dynamic_cast<const MNeg*     > (pi));
	case MADD     : return OpAdd    (dynamic_cast<const MAdd*     > (pi));
	case MSUB     : return OpSub    (dynamic_cast<const MSub*     > (pi));
	case MMUL     : return OpMul    (dynamic_cast<const MMul*     > (pi));
	case MDIV     : return OpDiv    (dynamic_cast<const MDiv*     > (pi));
	case MPOW     : return OpPow    (dynamic_cast<const MPow*     > (pi));
	case MEQUATION: return OpEqual  (dynamic_cast<const MEquation*> (pi));
	case MF1D     : return OpFnc1D  (dynamic_cast<const MFunc1D*  > (pi));
	case MF2D     : return OpFnc2D  (dynamic_cast<const MFunc2D*  > (pi));
	case MFND     : return OpFncND  (dynamic_cast<const MFuncND*  > (pi));
	case MSFNC    : return OpSFnc   (dynamic_cast<const MSFuncND* > (pi));
	}
	assert(false);
	return "";
}

//-----------------------------------------------------------------------------
// convert a constant to a string
string MObj2String::Constant(const MConstant *pc)
{
	stringstream ss;
	ss << pc->value();
	return ss.str();
}

//-----------------------------------------------------------------------------
// convert a fraction to a string
string MObj2String::Fraction(const MFraction *pc)
{
	FRACTION a = pc->fraction();
	stringstream ss;
	ss << a.n << "/" << a.d;
	return ss.str();
}

//-----------------------------------------------------------------------------
// convert a constant to a string
string MObj2String::NamedCt(const MNamedCt *pc)
{
	return pc->Name();
}

//-----------------------------------------------------------------------------
// convert a variable to a string
string MObj2String::Variable(const MVarRef* pv)
{
	const MVariable* pvar = pv->GetVariable();
	return string(pvar->Name());
}

//-----------------------------------------------------------------------------
string MObj2String::OpNeg(const MNeg* po)
{
	const MItem* pi = po->Item();
	string s = Convert(pi);
#ifdef _DEBUG
	s = "-[" + s + "]";
#else
	if (is_add(pi) || is_sub(pi) || is_neg(pi)) s = "-(" + s + ")"; else s = "-" + s;
#endif
	return s;
}

//-----------------------------------------------------------------------------
string MObj2String::OpAdd(const MAdd* po)
{
	string s;
	const MItem* pl = po->LeftItem ();
	const MItem* pr = po->RightItem();
	string sl = Convert(pl);
	string sr = Convert(pr);
#ifdef _DEBUG
	if (is_binary(pl)) s = "(" + sl + ")"; else s = sl;
	s += '+';
	if (is_binary(pr)) s += "(" + sr + ")"; else s += sr;
#else
	if (is_neg(pl)) s = '(' + sl + ')'; else s = sl; 
	s += '+';
	if (is_neg(pr)) s += '(' + sr + ')'; else s += sr; 
#endif
	return s;
}

//-----------------------------------------------------------------------------
string MObj2String::OpSub(const MSub* po)
{
	string s;
	const MItem* pl = po->LeftItem ();
	const MItem* pr = po->RightItem();
	string sl = Convert(pl);
	string sr = Convert(pr);
#ifdef _DEBUG
	if (is_binary(pl)) s = "(" + sl + ")"; else s = sl;
	s += '-';
	if (is_binary(pr)) s += "(" + sr + ")"; else s += sr;
#else
	if (is_neg(pl)) s = '(' + sl + ')'; else s = sl; 
	s += '-';
	if (is_add(pr) || is_sub(pr) || is_neg(pr)) s += '(' + sr + ')'; else s += sr; 
#endif
	return s;
}

//-----------------------------------------------------------------------------
string MObj2String::OpMul(const MMul* po)
{
	string s;
	const MItem* pl = po->LeftItem ();
	const MItem* pr = po->RightItem();
	string sl = Convert(pl);
	string sr = Convert(pr);
#ifdef _DEBUG
	if (is_binary(pl)) s = "(" + sl + ")"; else s = sl;
	s += '*';
	if (is_binary(pr)) s += "(" + sr + ")"; else s += sr;
#else
	if (is_add(pl) || is_sub(pl) || is_neg(pl) || is_div(pl)) s = '(' + sl + ')'; else s = sl;
	s += '*';
	if (is_add(pr) || is_sub(pr) || is_neg(pr) || is_div(pr)) s += '(' + sr + ')'; else s += sr;
#endif
	return s;
}

//-----------------------------------------------------------------------------
string MObj2String::OpDiv(const MDiv* po)
{
	string s;
	const MItem* pl = po->LeftItem ();
	const MItem* pr = po->RightItem();
	string sl = Convert(pl);
	string sr = Convert(pr);
#ifdef _DEBUG
	if (is_binary(pl)) s = "(" + sl + ")"; else s = sl;
	s += '/';
	if (is_binary(pr)) s += "(" + sr + ")"; else s += sr;
#else
	if (is_add(pl) || is_sub(pl) || is_neg(pl) || is_mul(pl) || is_div(pl) || is_pow(pr)) s = '(' + sl + ')'; else s = sl;
	s += '/';
	if (is_add(pr) || is_sub(pr) || is_neg(pr) || is_mul(pr) || is_div(pr)) s += '(' + sr + ')'; else s += sr;
#endif
	return s;
}

//-----------------------------------------------------------------------------
string MObj2String::OpPow(const MPow* po)
{
	string s;
	const MItem* pl = po->LeftItem ();
	const MItem* pr = po->RightItem();
	string sl = Convert(pl);
	string sr = Convert(pr);
#ifdef _DEBUG
	if (is_binary(pl)) s = "(" + sl + ")"; else s = sl;
	s += '^';
	if (is_binary(pr)) s += "(" + sr + ")"; else s += sr;
#else
	if (is_add(pl) || is_sub(pl) || is_neg(pl) || is_pow(pl) || is_mul(pl) || is_div(pl)) s = '(' + sl + ')'; else s = sl; 
	s += '^';
	if (is_add(pr) || is_sub(pr) || is_neg(pr) || is_pow(pr) || is_mul(pr) || is_div(pr)) s += '(' + sr + ')'; else s += sr;
#endif
	return s;
}

//-----------------------------------------------------------------------------
string MObj2String::OpEqual(const MEquation* po)
{
	string s;
	const MItem* pl = po->LeftItem ();
	const MItem* pr = po->RightItem();
	string sl = Convert(pl);
	string sr = Convert(pr);
#ifdef _DEBUG
	if (is_binary(pl)) s = "(" + sl + ")"; else s = sl;
	s += '=';
	if (is_binary(pr)) s += "(" + sr + ")"; else s += sr;
#else
	s += '=';
#endif
	return s;
}

//-----------------------------------------------------------------------------
string MObj2String::OpFnc1D(const MFunc1D *po)
{
	return string(po->Name() + "(" + Convert(po->Item()) + ")");
}

//-----------------------------------------------------------------------------
string MObj2String::OpFnc2D(const MFunc2D *po)
{
	return string(po->Name() + "(" + Convert(po->LeftItem()) + "," + Convert(po->RightItem()) + ")");
}

//-----------------------------------------------------------------------------
string MObj2String::OpFncND(const MFuncND *po)
{
	int N = po->Params();
	string s = po->Name() + "(";
	for (int i=0; i<N; ++i)
	{
		const MItem* pi = po->Param(i);
		s += Convert(pi);
		if (i != N-1) s += ",";
	}
	return s + ")";
}

//-----------------------------------------------------------------------------
string MObj2String::OpSFnc(const MSFuncND *po)
{
	int N = po->Params();
	string s = po->Name() + "(";
	for (int i=0; i<N; ++i)
	{
		const MItem* pi = po->Param(i);
		s += Convert(pi);
		if (i != N-1) s += ",";
	}
	return s + ")";
}
