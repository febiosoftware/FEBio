#include "MObj2String.h"
#include <sstream>

//-----------------------------------------------------------------------------
// Convert a list of math objects to a string
string MObj2String::Convert(MObjList* pol)
{
	string s;
	for (int i=0; i<pol->Objects(); ++i)
	{
		MathObject* po = pol->Object(i);
		s += Convert(*po);
		if (i != pol->Objects()-1) s += ",";
	}
	return s;
}

//-----------------------------------------------------------------------------
// Convert a math object to a string
string MObj2String::Convert(MathObject& o)
{
	MathObject* po = &o;
	if (dynamic_cast<MSimpleExpression*>(po))
	{
		MSimpleExpression* pe = dynamic_cast<MSimpleExpression*>(po);
		MITEM& i = pe->GetExpression();
		return Convert(i.ItemPtr());
	}
	assert(false);
	return "";
}

//-----------------------------------------------------------------------------
// Convert a math item to string
string MObj2String::Convert(MItem* pi)
{
	switch(pi->Type())
	{
	case MCONST   : return Constant (dynamic_cast<MConstant*> (pi));
	case MFRAC    : return Fraction (dynamic_cast<MFraction*> (pi));
	case MNAMED   : return NamedCt  (dynamic_cast<MNamedCt* > (pi));
	case MVAR     : return Variable (dynamic_cast<MVarRef*  > (pi));
	case MNEG     : return OpNeg    (dynamic_cast<MNeg*     > (pi));
	case MADD     : return OpAdd    (dynamic_cast<MAdd*     > (pi));
	case MSUB     : return OpSub    (dynamic_cast<MSub*     > (pi));
	case MMUL     : return OpMul    (dynamic_cast<MMul*     > (pi));
	case MDIV     : return OpDiv    (dynamic_cast<MDiv*     > (pi));
	case MPOW     : return OpPow    (dynamic_cast<MPow*     > (pi));
	case MEQUATION: return OpEqual  (dynamic_cast<MEquation*   > (pi));
	case MF1D     : return OpFnc1D  (dynamic_cast<MFunc1D*  > (pi));
	case MF2D     : return OpFnc2D  (dynamic_cast<MFunc2D*  > (pi));
	case MFND     : return OpFncND  (dynamic_cast<MFuncND*  > (pi));
	case MSFNC    : return OpSFnc   (msfncnd(pi));
	}
	assert(false);
	return "";
}

//-----------------------------------------------------------------------------
// convert a constant to a string
string MObj2String::Constant(MConstant *pc)
{
	stringstream ss;
	ss << pc->value();
	return ss.str();
}

//-----------------------------------------------------------------------------
// convert a fraction to a string
string MObj2String::Fraction(MFraction *pc)
{
	FRACTION a = pc->fraction();
	stringstream ss;
	ss << a.n << "/" << a.d;
	return ss.str();
}

//-----------------------------------------------------------------------------
// convert a constant to a string
string MObj2String::NamedCt(MNamedCt *pc)
{
	return pc->Name();
}

//-----------------------------------------------------------------------------
// convert a variable to a string
string MObj2String::Variable(MVarRef* pv)
{
	MVariable* pvar = pv->GetVariable();
	return string(pvar->Name());
}

//-----------------------------------------------------------------------------
string MObj2String::OpNeg(MNeg* po)
{
	MItem* pi = po->Item();
	string s = Convert(pi);
#ifdef _DEBUG
	s = "-[" + s + "]";
#else
	if (is_add(pi) || is_sub(pi) || is_neg(pi)) s = "-(" + s + ")"; else s = "-" + s;
#endif
	return s;
}

//-----------------------------------------------------------------------------
string MObj2String::OpAdd(MAdd* po)
{
	string s;
	MItem* pl = po->LeftItem ();
	MItem* pr = po->RightItem();
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
string MObj2String::OpSub(MSub* po)
{
	string s;
	MItem* pl = po->LeftItem ();
	MItem* pr = po->RightItem();
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
string MObj2String::OpMul(MMul* po)
{
	string s;
	MItem* pl = po->LeftItem ();
	MItem* pr = po->RightItem();
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
string MObj2String::OpDiv(MDiv* po)
{
	string s;
	MItem* pl = po->LeftItem ();
	MItem* pr = po->RightItem();
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
string MObj2String::OpPow(MPow* po)
{
	string s;
	MItem* pl = po->LeftItem ();
	MItem* pr = po->RightItem();
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
string MObj2String::OpEqual(MEquation* po)
{
	string s;
	MItem* pl = po->LeftItem ();
	MItem* pr = po->RightItem();
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
string MObj2String::OpFnc1D(MFunc1D *po)
{
	return string(po->Name() + "(" + Convert(po->Item()) + ")");
}

//-----------------------------------------------------------------------------
string MObj2String::OpFnc2D(MFunc2D *po)
{
	return string(po->Name() + "(" + Convert(po->LeftItem()) + "," + Convert(po->RightItem()) + ")");
}

//-----------------------------------------------------------------------------
string MObj2String::OpFncND(MFuncND *po)
{
	int N = po->Params();
	string s = po->Name() + "(";
	for (int i=0; i<N; ++i)
	{
		MItem* pi = po->Param(i);
		s += Convert(pi);
		if (i != N-1) s += ",";
	}
	return s + ")";
}

//-----------------------------------------------------------------------------
string MObj2String::OpSFnc(MSFuncND *po)
{
	int N = po->Params();
	string s = po->Name() + "(";
	for (int i=0; i<N; ++i)
	{
		MItem* pi = po->Param(i);
		s += Convert(pi);
		if (i != N-1) s += ",";
	}
	return s + ")";
}
