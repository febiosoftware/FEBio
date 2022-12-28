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



#pragma once
#include "MFunctions.h"
#include "MTypes.h"
#include <vector>
#include <assert.h>
#include <string>
#include <string.h>

//-----------------------------------------------------------------------------
// This exception will get thrown when an invalid operation is performed.
class InvalidOperation {};
class DivisionByZero   {};

//-----------------------------------------------------------------------------
// The following classes define the componenents of a mathematical expression
enum Item_Type {
	MCONST, MFRAC, MNAMED, MVAR, MBOOL,
	MMATRIX,
	MNEG,
	MADD, MSUB, MMUL, MDIV, MPOW,
	MEQUATION, MEQUALITY,
	MF1D, MF2D, MFND, MFMAT, MFMAT2,
	MSFNC,
	MINTEGRAL,
	MSEQUENCE
};

//-----------------------------------------------------------------------------
// stores the value of a variable
class MVariable
{
public:
	MVariable(const std::string& s, double initVal = 0.0) : m_v(initVal), m_name(s), m_index(-1) {}
	MVariable(const MVariable& v) : m_v(0), m_name(v.m_name), m_index(-1) {}
	void operator = (const MVariable& v) { m_name = v.m_name; m_index = v.m_index; }
	double value() const { return m_v; }
	void value(double v) { m_v = v; }
	const std::string& Name() const { return m_name; }
	void setIndex(int n) { m_index = n; }
	int index() const { return m_index; }

protected:
	double		m_v;
	int			m_index;
	std::string	m_name;
};

//-----------------------------------------------------------------------------
// base class for Math items
// a math item can be a constant, a variable, a matrix, an operator or a function
class MItem
{
public:
	MItem(Item_Type ntype) : m_ntype(ntype) {}
	virtual ~MItem(){}

	// return the item's type
	Item_Type Type() const { return m_ntype; }

	// create a copy of this item
	virtual MItem* copy() const = 0;

protected:
	Item_Type m_ntype;
};

//-----------------------------------------------------------------------------
// base class for numbers.
// numbers can be either constant or variables.
class MNumber : public MItem
{
public: 
	MNumber(Item_Type ntype) : MItem(ntype){}

	// return the float value for this number
	virtual double value() const = 0;
};

//-----------------------------------------------------------------------------
// defines a constant number
class MConstant : public MNumber
{
public:
	MConstant(double x) : MNumber(MCONST), m_v(x) { assert(m_v>=0); }
	void value(double v) { m_v = v; }

public:
	double value() const override { return m_v; }
	MItem* copy() const override { return new MConstant(m_v); }

private:
	double	m_v;
};

//-----------------------------------------------------------------------------
// defines a constant fraction
class MFraction : public MNumber
{
public:
	MFraction(double a, double b) : MNumber(MFRAC), m_v(a,b) {}
	MFraction(const FRACTION& a) : MNumber(MFRAC), m_v(a) {}
	FRACTION fraction() const { return m_v; }

public:
	double value() const override { return m_v.n / m_v.d; }
	MItem* copy() const override { return new MFraction(m_v); }

private:
	FRACTION	m_v;
};

//-----------------------------------------------------------------------------
// defines a named constant number
class MNamedCt : public MNumber
{
public:
	MNamedCt(double x, const std::string& sz) : MNumber(MNAMED), m_v(x), m_name(sz) {}
	const std::string& Name() const { return m_name; }

public:
	double value() const override { return m_v; }
	MItem* copy() const override { return new MNamedCt(m_v, m_name); }

private:
	double			m_v;
	std::string		m_name;
};

//-----------------------------------------------------------------------------
// defines a reference to a variable
class MVarRef : public MNumber
{
public:
	MVarRef(const MVariable* pv) : MNumber(MVAR), m_pv(pv) {}
	const MVariable* GetVariable() const { return m_pv; }
	const std::string& Name() const { return m_pv->Name(); }
	int index() const { return m_pv->index(); }

	void SetVariable(const MVariable* v) { m_pv = v; }

public:
	double value() const override { return m_pv->value(); }
	MItem* copy() const override { return new MVarRef(m_pv); }

private:
	const MVariable*	m_pv;
};

//-----------------------------------------------------------------------------
// Item for representing boolean values
class MBoolean : public MItem
{
public:
	MBoolean(bool b) : MItem(MBOOL), m_b(b) {}
	MItem* copy() const override { return new MBoolean(m_b); }

	bool value() { return m_b; }

private:
	bool	m_b;
};

//-----------------------------------------------------------------------------
// base class for unary operations
class MUnary : public MItem
{
public:
	MUnary(MItem* pi, Item_Type ntype) : MItem(ntype), m_pi(pi) {}
	~MUnary() { delete m_pi; }
	
	MItem* Item() { return m_pi; }
	const MItem* Item() const { return m_pi; }

protected:
	MItem*	m_pi;
};

//-----------------------------------------------------------------------------
// base class for binary operations
class MBinary : public MItem
{
public:
	MBinary(MItem* pl, MItem* pr, Item_Type ntype) : MItem(ntype), m_pleft(pl), m_pright(pr) {}
	~MBinary() { delete m_pleft; delete m_pright; }

	MItem* LeftItem() { return m_pleft; }
	MItem* RightItem() { return m_pright; }

	const MItem* LeftItem() const { return m_pleft; }
	const MItem* RightItem() const { return m_pright; }

protected:
	MItem *m_pleft, *m_pright;
};

//-----------------------------------------------------------------------------
// a list of items
class MSequence : public MItem
{
public:
	MSequence() : MItem(MSEQUENCE){}
	MSequence(const MSequence& s);
	~MSequence();

	MSequence* copy() const override { return new MSequence(*this); }
	MSequence& operator = (MSequence& s);

	MItem* operator [] (int i) { return m_item[i]; }
	const MItem* operator [] (int i) const { return m_item[i]; }

	int size() const { return (int) m_item.size(); }

	void add(MItem* pi) { m_item.push_back(pi); }

	void replace(int i, MItem* pi);

	void remove(int i);

protected:
	std::vector<MItem*>	m_item;
};

//-----------------------------------------------------------------------------
// base class for N-ary operators and functions
class MNary : public MItem
{
public:
	MNary(const MSequence& s, Item_Type ntype) : MItem(ntype), m_s(s) {}
	~MNary() {}
	int Params() const { return (int) m_s.size(); }
	MItem* Param(int n) { return m_s[n]; }
	const MItem* Param(int n) const { return m_s[n]; }
	MSequence& GetSequence() { return m_s; }
	const MSequence& GetSequence() const { return m_s; }

protected:
	MSequence	m_s;
};

//-----------------------------------------------------------------------------
// defines a negative operation
class MNeg : public MUnary
{
public: 
	MNeg(MItem* pi) : MUnary(pi, MNEG) {}
	MItem* copy() const override { return new MNeg(m_pi->copy()); }
};

//-----------------------------------------------------------------------------
// addition operator
class MAdd : public MBinary
{
public:
	MAdd(MItem* pl, MItem* pr) : MBinary(pl, pr, MADD){}
	MAdd(MItem* pl, double  r) : MBinary(pl, new MConstant(r), MADD){}
	MItem* copy() const override { return new MAdd(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// subtraction operator
class MSub : public MBinary
{
public:
	MSub(MItem* pl, MItem* pr) : MBinary(pl, pr, MSUB){}
	MSub(double  l, MItem* pr) : MBinary(new MConstant(l), pr, MSUB){}
	MSub(MItem* pl, double  r) : MBinary(pl, new MConstant(r), MSUB){}
	MItem* copy() const override { return new MSub(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// multiplication operator
class MMul : public MBinary
{
public:
	MMul(MItem* pl, MItem* pr) : MBinary(pl, pr, MMUL){}
	MMul(double  l, MItem* pr) : MBinary(new MConstant(l), pr, MMUL){}
	MItem* copy() const override { return new MMul(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// division operator
class MDiv : public MBinary
{
public:
	MDiv(MItem* pl, MItem* pr) : MBinary(pl, pr, MDIV){}
	MDiv(double  l, MItem* pr) : MBinary(new MConstant(l), pr, MDIV){}
	MDiv(MItem* pl, double  r) : MBinary(pl, new MConstant(r), MDIV){}
	MDiv(double  l, double  r) : MBinary(new MConstant(l), new MConstant(r), MDIV){}
	MItem* copy() const override { return new MDiv(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// power operator
class MPow : public MBinary
{
public:
	MPow(MItem* pl, MItem* pr) : MBinary(pl, pr, MPOW){}
	MPow(MItem* pl, double  r) : MBinary(pl, new MConstant(r), MPOW){}
	MItem* copy() const override { return new MPow(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// equal-than operator
class MEquation : public MBinary
{
public:
	MEquation(MItem* pl, MItem* pr) : MBinary(pl, pr, MEQUATION){}
	MItem* copy() const override { return new MEquation(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// equality test
class MEquality : public MBinary
{
public:
	MEquality(MItem* pl, MItem* pr) : MBinary(pl, pr, MEQUALITY){}
	MItem* copy() const override { return new MEquality(m_pleft->copy(), m_pright->copy()); }
};


//-----------------------------------------------------------------------------
// function of one variable
class MFunc1D : public MUnary
{
public: 
	MFunc1D(FUNCPTR pf, const std::string& s, MItem* pi) : MUnary(pi, MF1D), m_name(s), m_pf(pf) {}
	MItem* copy() const override { return new MFunc1D(m_pf, m_name, m_pi->copy()); }
	const std::string& Name() const { return m_name; }
	FUNCPTR funcptr() const { return m_pf; }

protected:
	std::string	m_name;
	FUNCPTR		m_pf;
};

//-----------------------------------------------------------------------------
// function of two variables
class MFunc2D : public MBinary
{
public: 
	MFunc2D(FUNC2PTR pf, const std::string& s, MItem* p1, MItem* p2) : MBinary(p1, p2, MF2D), m_name(s), m_pf(pf)  {}
	MItem* copy() const override { return new MFunc2D(m_pf, m_name, m_pleft->copy(), m_pright->copy()); }
	const std::string& Name() const { return m_name; }
	FUNC2PTR	funcptr() const { return m_pf; }

protected:
	std::string		m_name;
	FUNC2PTR		m_pf;
};

//-----------------------------------------------------------------------------
// function of N variables
class MFuncND : public MNary
{
public: 
	MFuncND(FUNCNPTR pf, const std::string& s, const MSequence& l) : MNary(l, MFND), m_name(s), m_pf(pf) {}
	MItem* copy() const override;
	const std::string& Name() const { return m_name; }

protected:
	std::string		m_name;
	FUNCNPTR		m_pf;
};

//-----------------------------------------------------------------------------
// class for a symbolic function definition of N variables
class MSFuncND : public MNary
{
public:
	MSFuncND(const std::string& s, MItem* pv, const MSequence* ps) : MNary(*ps, MSFNC), m_name(s), m_pval(pv) {}
	~MSFuncND(){ delete m_pval; }
	MItem* copy() const override { return new MSFuncND(m_name, m_pval->copy(), &GetSequence()); }

	const MItem* Value() const { return m_pval; }

	const std::string& Name() const { return m_name; }

private:
	std::string		m_name;
	MItem*			m_pval;	//!< result of function evaluation
};


//-----------------------------------------------------------------------------
// Symbolic operator
class MSymOp : public MItem
{
public:
	MSymOp(Item_Type ntype) : MItem(ntype) {}
};

//-----------------------------------------------------------------------------
class MOpIntegral : public MSymOp
{
public:
	MOpIntegral(MItem* pf, MVarRef* px) : MSymOp(MINTEGRAL), m_pf(pf), m_px(px){}
	~MOpIntegral() { delete m_pf; delete m_px; }

	MItem* copy() const override { return new MOpIntegral(m_pf->copy(), dynamic_cast<MVarRef*>(m_px->copy())); }

	MItem* Item() { return m_pf; }
	MVarRef* Var() { return m_px; }

protected:
	MItem*		m_pf;
	MVarRef*	m_px;
};

//-----------------------------------------------------------------------------
inline const MNumber*   mnumber  (const MItem* pi) { return static_cast<const MNumber*  >(pi); }
inline const MConstant* mconst   (const MItem* pi) { return static_cast<const MConstant*>(pi); }
inline const MFraction* mfrac    (const MItem* pi) { return static_cast<const MFraction*>(pi); }
inline const MNamedCt*  mnamed   (const MItem* pi) { return static_cast<const MNamedCt *>(pi); }
inline const MVarRef*   mvar     (const MItem* pi) { return static_cast<const MVarRef  *>(pi); }
inline const MBoolean*  mboolean (const MItem* pi) { return static_cast<const MBoolean* >(pi); }
inline const MNeg*      mneg     (const MItem* pi) { return static_cast<const MNeg     *>(pi); }
inline const MAdd*      madd     (const MItem* pi) { return static_cast<const MAdd     *>(pi); }
inline const MSub*      msub     (const MItem* pi) { return static_cast<const MSub     *>(pi); }
inline const MMul*      mmul     (const MItem* pi) { return static_cast<const MMul     *>(pi); }
inline const MDiv*      mdiv     (const MItem* pi) { return static_cast<const MDiv     *>(pi); }
inline const MPow*      mpow     (const MItem* pi) { return static_cast<const MPow     *>(pi); }
inline const MUnary*    munary   (const MItem* pi) { return static_cast<const MUnary   *>(pi); }
inline const MBinary*   mbinary  (const MItem* pi) { return static_cast<const MBinary  *>(pi); }
inline const MNary*     mnary    (const MItem* pi) { return static_cast<const MNary    *>(pi); }
inline const MFunc1D*   mfnc1d   (const MItem* pi) { return static_cast<const MFunc1D  *>(pi); }
inline const MFunc2D*   mfnc2d   (const MItem* pi) { return static_cast<const MFunc2D  *>(pi); }
inline const MFuncND*   mfncnd   (const MItem* pi) { return static_cast<const MFuncND  *>(pi); }
inline const MSFuncND*  msfncnd  (const MItem* pi) { return static_cast<const MSFuncND *>(pi); }
inline const MSequence* msequence(const MItem* pi) { return static_cast<const MSequence*>(pi); }
inline const MEquation* mequation(const MItem* pi) { return static_cast<const MEquation*>(pi); }
inline const MEquality* mequality(const MItem* pi) { return static_cast<const MEquality*>(pi); }

//-----------------------------------------------------------------------------
// This class is defined so that we can do arithmetic with the  math classes 
// without worrying about making copies and garbage collection. Each MITEM 
// stores its own expression. This will probably cause a proliferation of 
// copying, so perhaps I can come up with a better memory managment model.
// (e.g. like auto_ptr, where owner- ship of an pointer can be passed along)
class MITEM
{
public:
	MITEM() { m_pi = 0; }
	~MITEM() { delete m_pi; }
	//~MITEM() { }

	MITEM(MItem* pi) { m_pi = pi; } // Note that the item is not copied in this constructor
	MITEM(const MItem* pi) { m_pi = pi->copy(); }
	MITEM(const MITEM& it) { m_pi = (it.m_pi?it.m_pi->copy():0); }
	MITEM(const MVariable& v) { m_pi = new MVarRef(&v); }
	MITEM(double g) 
	{
		m_pi = new MConstant(fabs(g));
		if (g < 0) m_pi = new MNeg(m_pi);
	}
	MITEM(FRACTION& g)
	{
		m_pi = new MFraction(fabs(g.n), fabs(g.d));
		if (g < 0) m_pi = new MNeg(m_pi);
	}

	void operator = (const MITEM& it) { delete m_pi; m_pi = it.m_pi->copy(); }
	void operator = (MItem* pi) { delete m_pi; m_pi = pi; }

	MItem* ItemPtr() { return m_pi; }
	const MItem* ItemPtr() const { return m_pi; }

	MItem* copy() const { return m_pi->copy(); }

	MITEM Left () const { return (static_cast<MBinary*>(m_pi))->LeftItem()->copy(); }
	MITEM Right() const { return (static_cast<MBinary*>(m_pi))->RightItem()->copy(); }

	MITEM Item() const { return (static_cast<MUnary*>(m_pi))->Item()->copy(); }

	MITEM Param() const { return (static_cast<MFunc1D*>(m_pi))->Item()->copy(); }

	Item_Type Type() const { return m_pi->Type(); }

	double value() const { return (Type()==MNEG? -mnumber(mneg(m_pi)->Item())->value() : mnumber(m_pi)->value()); }

private:
	MItem*	m_pi;
};

//-----------------------------------------------------------------------------
inline const MNumber*   mnumber  (const MITEM& i) { return static_cast<const MNumber  *>(i.ItemPtr()); }
inline const MConstant* mconst   (const MITEM& i) { return static_cast<const MConstant*>(i.ItemPtr()); }
inline const MFraction* mfrac    (const MITEM& i) { return static_cast<const MFraction*>(i.ItemPtr()); }
inline const MNamedCt*  mnamed   (const MITEM& i) { return static_cast<const MNamedCt *>(i.ItemPtr()); }
inline const MVarRef*   mvar     (const MITEM& i) { return static_cast<const MVarRef  *>(i.ItemPtr()); }
inline const MBoolean*  mboolean (const MITEM& i) { return static_cast<const MBoolean* >(i.ItemPtr()); }
inline const MNeg*      mneg     (const MITEM& i) { return static_cast<const MNeg     *>(i.ItemPtr()); }
inline const MAdd*      madd     (const MITEM& i) { return static_cast<const MAdd     *>(i.ItemPtr()); }
inline const MSub*      msub     (const MITEM& i) { return static_cast<const MSub     *>(i.ItemPtr()); }
inline const MMul*      mmul     (const MITEM& i) { return static_cast<const MMul     *>(i.ItemPtr()); }
inline const MDiv*      mdiv     (const MITEM& i) { return static_cast<const MDiv     *>(i.ItemPtr()); }
inline const MPow*      mpow     (const MITEM& i) { return static_cast<const MPow     *>(i.ItemPtr()); }
inline const MUnary*    munary   (const MITEM& i) { return static_cast<const MUnary   *>(i.ItemPtr()); }
inline const MBinary*   mbinary  (const MITEM& i) { return static_cast<const MBinary  *>(i.ItemPtr()); }
inline const MNary*     mnary    (const MITEM& i) { return static_cast<const MNary    *>(i.ItemPtr()); }
inline const MFunc1D*   mfnc1d   (const MITEM& i) { return static_cast<const MFunc1D  *>(i.ItemPtr()); }
inline const MFunc2D*   mfnc2d   (const MITEM& i) { return static_cast<const MFunc2D  *>(i.ItemPtr()); }
inline const MSFuncND*  msfncnd  (const MITEM& i) { return static_cast<const MSFuncND *>(i.ItemPtr()); }
inline const MSequence* msequence(const MITEM& i) { return static_cast<const MSequence*>(i.ItemPtr()); }
inline const MEquation* mequation(const MITEM& i) { return static_cast<const MEquation*>(i.ItemPtr()); }
inline const MEquality* mequality(const MITEM& i) { return static_cast<const MEquality*>(i.ItemPtr()); }

//-----------------------------------------------------------------------------
// Some helper functions to determine the type of an MItem
inline bool is_number(const MItem* pi) { return (dynamic_cast<const MNumber  *>(pi) != 0); }
inline bool is_unary (const MItem* pi) { return (dynamic_cast<const MUnary   *>(pi) != 0); }
inline bool is_binary(const MItem* pi) { return (dynamic_cast<const MBinary  *>(pi) != 0); }
inline bool is_nary  (const MItem* pi) { return (dynamic_cast<const MNary    *>(pi) != 0); }

inline bool isConst    (const MItem* pi) { return (pi->Type() == MCONST   ); }
inline bool is_frac    (const MItem* pi) { return (pi->Type() == MFRAC    ); }
inline bool is_named   (const MItem* pi) { return (pi->Type() == MNAMED   ); }
inline bool is_var     (const MItem* pi) { return (pi->Type() == MVAR     ); }
inline bool is_boolean (const MItem* pi) { return (pi->Type() == MBOOL    ); }
inline bool is_neg     (const MItem* pi) { return (pi->Type() == MNEG     ); }
inline bool is_add     (const MItem* pi) { return (pi->Type() == MADD     ); }
inline bool is_sub     (const MItem* pi) { return (pi->Type() == MSUB     ); }
inline bool is_mul     (const MItem* pi) { return (pi->Type() == MMUL     ); }
inline bool is_div     (const MItem* pi) { return (pi->Type() == MDIV     ); }
inline bool is_pow     (const MItem* pi) { return (pi->Type() == MPOW     ); }
inline bool is_equation(const MItem* pi) { return (pi->Type() == MEQUATION); }
inline bool is_equality(const MItem* pi) { return (pi->Type() == MEQUALITY); }
inline bool is_func1d  (const MItem* pi) { return (pi->Type() == MF1D     ); }
inline bool is_matrix  (const MItem* pi) { return (pi->Type() == MMATRIX  ); }
inline bool is_sfunc   (const MItem* pi) { return (pi->Type() == MSFNC    ); }

inline bool isConst    (const MITEM& i) { return (i.Type() == MCONST   ); }
inline bool is_frac    (const MITEM& i) { return (i.Type() == MFRAC    ); }
inline bool is_named   (const MITEM& i) { return (i.Type() == MNAMED   ); }
inline bool is_var     (const MITEM& i) { return (i.Type() == MVAR     ); }
inline bool is_boolean (const MITEM& i) { return (i.Type() == MBOOL    ); }
inline bool is_neg     (const MITEM& i) { return (i.Type() == MNEG     ); }
inline bool is_add     (const MITEM& i) { return (i.Type() == MADD     ); }
inline bool is_sub     (const MITEM& i) { return (i.Type() == MSUB     ); }
inline bool is_mul     (const MITEM& i) { return (i.Type() == MMUL     ); }
inline bool is_div     (const MITEM& i) { return (i.Type() == MDIV     ); }
inline bool is_pow     (const MITEM& i) { return (i.Type() == MPOW     ); }
inline bool is_equation(const MITEM& i) { return (i.Type() == MEQUATION); }
inline bool is_equality(const MITEM& i) { return (i.Type() == MEQUALITY); }
inline bool is_func1d  (const MITEM& i) { return (i.Type() == MF1D     ); }
inline bool is_matrix  (const MITEM& i) { return (i.Type() == MMATRIX  ); }
inline bool is_sequence(const MITEM& i) { return (i.Type() == MSEQUENCE); }
inline bool is_sfunc   (const MITEM& i) { return (i.Type() == MSFNC    ); }

inline bool is_number(const MITEM& i) { return is_number(i.ItemPtr()); }

inline bool is_func(const MItem* pi) { Item_Type t = pi->Type(); return ((t==MF1D) || (t==MF2D) || (t==MFND)); }
inline bool is_int (double  r) { return (r - floor(r) == 0.0); }
inline bool is_abs (const MItem* pi) { return (is_func1d(pi) && (mfnc1d(pi)->Name().compare("abs") == 0)); }

inline bool is_func(const MITEM& i) { return is_func(i.ItemPtr()); }
inline bool is_abs (const MITEM& i) { return is_abs (i.ItemPtr()); }

inline bool is_int(const MItem* pi) { return isConst(pi) && is_int(mnumber(pi)->value()); }
inline bool is_int (const MITEM& i) { return is_int(i.ItemPtr()); }

inline bool is_fnc1d(const MITEM& i, const char* sz) { return (is_func1d(i) && (strcmp(mfnc1d(i)->Name().c_str(), sz) == 0)); }

bool is_rconst(const MItem* pi); // is recursive constant
bool is_dependent(const MItem* pi, const MVariable& x); // is i dependent on x
bool is_dependent(const MItem* pi, const MItem* px); // is i dependent on x
bool is_scalar(const MItem* pi);	// is the expression scalar

inline bool is_rconst(const MITEM& i) { return is_rconst(i.ItemPtr()); }

bool is_pi(const MItem* pi); // is pi the constant pi?
inline bool is_pi(const MITEM& i) { return is_pi(i.ItemPtr()); }

inline bool is_zero(const MItem* pi) { if (isConst(pi)) return (mconst(pi)->value() == 0.0); else return false; }

inline bool is_scalar(const MITEM& m) { return ::is_scalar(m.ItemPtr()); }

//-----------------------------------------------------------------------------
// Algebraic operators for MITEM's
MITEM operator - (const MITEM& l);
MITEM operator + (const MITEM& l, const MITEM& r);
MITEM operator + (const MITEM& l, double r);
MITEM operator + (double l, const MITEM& r);
MITEM operator - (const MITEM& l, const MITEM& r);
MITEM operator - (double l, const MITEM& r);
MITEM operator - (const MITEM& l, double r);
MITEM operator * (const MITEM& l, const MITEM& r);
MITEM operator * (double l, const MITEM& r);
MITEM operator / (const MITEM& l, const MITEM& r);
MITEM operator / (double l, const MITEM& r);
MITEM operator / (const MITEM& l, double r);
MITEM operator ^ (const MITEM& l, const MITEM& r);
MITEM operator ^ (const MITEM& l, double r);

//-----------------------------------------------------------------------------
// Fractions
MITEM Fraction(double n, double d);
inline MITEM Fraction(FRACTION f) { return Fraction(f.n, f.d); }

//-----------------------------------------------------------------------------
// Symbolic math functions
MITEM Abs(const MITEM& l);
MITEM Sgn(const MITEM& l);
MITEM Sin(const MITEM& l);
MITEM Cos(const MITEM& l);
MITEM Tan(const MITEM& l);
MITEM Cot(const MITEM& l);
MITEM Sec(const MITEM& l);
MITEM Csc(const MITEM& l);
MITEM Atan(const MITEM& l);
MITEM Cosh(const MITEM& l);
MITEM Sinh(const MITEM& l);
MITEM Exp(const MITEM& l);
MITEM Log(const MITEM& l);
MITEM Log10(const MITEM& l);
MITEM Float(const MITEM& l);
MITEM Sqrt(const MITEM& l);
MITEM Erf(const MITEM& l);
MITEM Erfc(const MITEM& l);

// Chebyshev polynomials
MITEM Tn(int n, MITEM& r);
MITEM Tn(const MITEM& l, const MITEM& r);

// Bessel functions (1st kind)
#ifdef WIN32
MITEM J0(const MITEM& l);
MITEM J1(const MITEM& l);
MITEM Jn(int n, const MITEM& r);
MITEM Jn(const MITEM& l, const MITEM& r);
#endif 

// Bessel functions (2e kind)
#ifdef WIN32
MITEM Y0(const MITEM& l);
MITEM Y1(const MITEM& l);
MITEM Yn(int n, const MITEM& r);
MITEM Yn(const MITEM& l, const MITEM& r);
#endif

// Miscellaneous functions
MITEM Fac(const MITEM& l);
MITEM Binomial(const MITEM& l, const MITEM& r);
MITEM Identity(const MITEM& n);

//-----------------------------------------------------------------------------
// check equality of two items
inline bool operator == (const MITEM& l, const std::string& r) { if (is_var(l)) { const std::string& s = mvar(l)->Name(); return (s.compare(r) == 0); } else return false; }
inline bool operator == (const MITEM& l, const MVariable& x) { if (is_var(l)) { const std::string& s = mvar(l)->Name(); return (s.compare(x.Name()) == 0); } else return false; }
inline bool operator == (const MITEM& l, double r) { if (isConst(l)) return (mnumber(l)->value() == r); else return false; }
inline bool operator != (const MITEM& l, double r) { return !(l==r); }
inline bool is_dependent(const MITEM& l, const MVariable& x) { return is_dependent(l.ItemPtr(), x); }
inline bool is_dependent(const MITEM& l, const MITEM& x) { return is_dependent(l.ItemPtr(), x.ItemPtr()); }

bool is_equal(const MItem* pl, const MItem* pr);
inline bool is_equal(const MItem* pl, const MVariable& x) { if (is_var(pl)) { const std::string& s = mvar(pl)->Name(); return (s.compare(x.Name()) == 0); } else return false; }
inline bool operator == (const MITEM& l, const MITEM& r) { return is_equal(l.ItemPtr(), r.ItemPtr()); }
inline bool operator != (const MITEM& l, const MITEM& r) { return (is_equal(l.ItemPtr(), r.ItemPtr()) == false); }

//-----------------------------------------------------------------------------
// functions for calculating and comparing heuristics
int op_count(MItem* pi);
inline int op_count(MITEM& i) { return op_count(i.ItemPtr()); }

//-----------------------------------------------------------------------------
bool is_format(MItem* pe, const char* sz);
inline bool is_format(MITEM& e, const char* sz) { return is_format(e.ItemPtr(), sz); }
