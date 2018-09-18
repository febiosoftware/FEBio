#pragma once
#include "MFunctions.h"
#include "MTypes.h"
#include <vector>
#include <assert.h>

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
	MVariable(const std::string& s) : m_v(0), m_name(s), m_index(-1) {}
	MVariable(const MVariable& v) : m_v(0), m_name(v.m_name), m_index(-1) {}
	void operator = (const MVariable& v) { m_name = v.m_name; m_index = v.m_index; }
	double value() { return m_v; }
	void value(double v) { m_v = v; }
	const std::string& Name() { return m_name; }
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
	virtual MItem* copy() = 0;

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
	virtual double value() = 0;
};

//-----------------------------------------------------------------------------
// defines a constant number
class MConstant : public MNumber
{
public:
	MConstant(double x) : MNumber(MCONST), m_v(x) { assert(m_v>=0); }
	void value(double v) { m_v = v; }

public:
	double value() { return m_v; }
	MItem* copy() { return new MConstant(m_v); }

private:
	double	m_v;
};

//-----------------------------------------------------------------------------
// defines a constant fraction
class MFraction : public MNumber
{
public:
	MFraction(double a, double b) : MNumber(MFRAC), m_v(a,b) {}
	MFraction(FRACTION& a) : MNumber(MFRAC), m_v(a) {}
	FRACTION fraction() { return m_v; }

public:
	double value() { return m_v.n / m_v.d; }
	MItem* copy() { return new MFraction(m_v); }

private:
	FRACTION	m_v;
};

//-----------------------------------------------------------------------------
// defines a named constant number
class MNamedCt : public MNumber
{
public:
	MNamedCt(double x, const std::string& sz) : MNumber(MNAMED), m_v(x), m_name(sz) {}
	const std::string& Name() { return m_name; }

public:
	double value() { return m_v; }
	MItem* copy() { return new MNamedCt(m_v, m_name); }

private:
	double			m_v;
	std::string		m_name;
};

//-----------------------------------------------------------------------------
// defines a reference to a variable
class MVarRef : public MNumber
{
public:
	MVarRef(MVariable* pv) : MNumber(MVAR), m_pv(pv) {}
	MVariable* GetVariable() { return m_pv; }
	const std::string& Name() { return m_pv->Name(); }
	int index() const { return m_pv->index(); }

public:
	double value() { return m_pv->value(); }
	MItem* copy() { return new MVarRef(m_pv); }

private:
	MVariable*	m_pv;
};

//-----------------------------------------------------------------------------
// Item for representing boolean values
class MBoolean : public MItem
{
public:
	MBoolean(bool b) : MItem(MBOOL), m_b(b) {}
	MItem* copy() { return new MBoolean(m_b); }

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

protected:
	MItem *m_pleft, *m_pright;
};

//-----------------------------------------------------------------------------
// a list of items
class MSequence : public MItem
{
public:
	MSequence() : MItem(MSEQUENCE){}
	MSequence(MSequence& s);
	~MSequence();

	MSequence* copy() { return new MSequence(*this); }
	MSequence& operator = (MSequence& s);

	MItem* operator [] (int i) { return m_item[i]; }

	int size() { return (int) m_item.size(); }

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
	MNary(MSequence& s, Item_Type ntype) : MItem(ntype), m_s(s) {}
	~MNary() {}
	int Params() { return (int) m_s.size(); }
	MItem* Param(int n) { return m_s[n]; }
	MSequence& GetSequence() { return m_s; }

protected:
	MSequence	m_s;
};

//-----------------------------------------------------------------------------
// defines a negative operation
class MNeg : public MUnary
{
public: 
	MNeg(MItem* pi) : MUnary(pi, MNEG) {}
	MItem* copy() { return new MNeg(m_pi->copy()); }
};

//-----------------------------------------------------------------------------
// addition operator
class MAdd : public MBinary
{
public:
	MAdd(MItem* pl, MItem* pr) : MBinary(pl, pr, MADD){}
	MAdd(MItem* pl, double  r) : MBinary(pl, new MConstant(r), MADD){}
	MItem* copy() { return new MAdd(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// subtraction operator
class MSub : public MBinary
{
public:
	MSub(MItem* pl, MItem* pr) : MBinary(pl, pr, MSUB){}
	MSub(double  l, MItem* pr) : MBinary(new MConstant(l), pr, MSUB){}
	MSub(MItem* pl, double  r) : MBinary(pl, new MConstant(r), MSUB){}
	MItem* copy() { return new MSub(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// multiplication operator
class MMul : public MBinary
{
public:
	MMul(MItem* pl, MItem* pr) : MBinary(pl, pr, MMUL){}
	MMul(double  l, MItem* pr) : MBinary(new MConstant(l), pr, MMUL){}
	MItem* copy() { return new MMul(m_pleft->copy(), m_pright->copy()); }
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
	MItem* copy() { return new MDiv(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// power operator
class MPow : public MBinary
{
public:
	MPow(MItem* pl, MItem* pr) : MBinary(pl, pr, MPOW){}
	MPow(MItem* pl, double  r) : MBinary(pl, new MConstant(r), MPOW){}
	MItem* copy() { return new MPow(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// equal-than operator
class MEquation : public MBinary
{
public:
	MEquation(MItem* pl, MItem* pr) : MBinary(pl, pr, MEQUATION){}
	MItem* copy() { return new MEquation(m_pleft->copy(), m_pright->copy()); }
};

//-----------------------------------------------------------------------------
// equality test
class MEquality : public MBinary
{
public:
	MEquality(MItem* pl, MItem* pr) : MBinary(pl, pr, MEQUALITY){}
	MItem* copy() { return new MEquality(m_pleft->copy(), m_pright->copy()); }
};


//-----------------------------------------------------------------------------
// function of one variable
class MFunc1D : public MUnary
{
public: 
	MFunc1D(FUNCPTR pf, const std::string& s, MItem* pi) : MUnary(pi, MF1D), m_name(s), m_pf(pf) {}
	MItem* copy() { return new MFunc1D(m_pf, m_name, m_pi->copy()); }
	const std::string& Name() { return m_name; }
	FUNCPTR	funcptr() { return m_pf; }

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
	MItem* copy() { return new MFunc2D(m_pf, m_name, m_pleft->copy(), m_pright->copy()); }
	const std::string& Name() { return m_name; }
	FUNC2PTR	funcptr() { return m_pf; }

protected:
	std::string		m_name;
	FUNC2PTR		m_pf;
};

//-----------------------------------------------------------------------------
// function of N variables
class MFuncND : public MNary
{
public: 
	MFuncND(FUNCNPTR pf, const std::string& s, MSequence& l) : MNary(l, MFND), m_name(s), m_pf(pf) {}
	MItem* copy();
	const std::string& Name() { return m_name; }

protected:
	std::string		m_name;
	FUNCNPTR		m_pf;
};

//-----------------------------------------------------------------------------
// class for a symbolic function definition of N variables
class MSFuncND : public MNary
{
public:
	MSFuncND(const std::string& s, MItem* pv, MSequence* ps) : MNary(*ps, MSFNC), m_name(s), m_pval(pv) {}
	~MSFuncND(){ delete m_pval; }
	MItem* copy() { return new MSFuncND(m_name, m_pval->copy(), &GetSequence()); }

	MItem* Value() { return m_pval; }

	const std::string& Name() { return m_name; }

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

	MItem* copy() { return new MOpIntegral(m_pf->copy(), dynamic_cast<MVarRef*>(m_px->copy())); }

	MItem* Item() { return m_pf; }
	MVarRef* Var() { return m_px; }

protected:
	MItem*		m_pf;
	MVarRef*	m_px;
};

//-----------------------------------------------------------------------------
inline MNumber*   mnumber  (MItem* pi) { return static_cast<MNumber*  >(pi); }
inline MConstant* mconst   (MItem* pi) { return static_cast<MConstant*>(pi); }
inline MFraction* mfrac    (MItem* pi) { return static_cast<MFraction*>(pi); }
inline MNamedCt*  mnamed   (MItem* pi) { return static_cast<MNamedCt *>(pi); }
inline MVarRef*   mvar     (MItem* pi) { return static_cast<MVarRef  *>(pi); }
inline MBoolean*  mboolean (MItem* pi) { return static_cast<MBoolean* >(pi); }
inline MNeg*      mneg     (MItem* pi) { return static_cast<MNeg     *>(pi); }
inline MAdd*      madd     (MItem* pi) { return static_cast<MAdd     *>(pi); }
inline MSub*      msub     (MItem* pi) { return static_cast<MSub     *>(pi); }
inline MMul*      mmul     (MItem* pi) { return static_cast<MMul     *>(pi); }
inline MDiv*      mdiv     (MItem* pi) { return static_cast<MDiv     *>(pi); }
inline MPow*      mpow     (MItem* pi) { return static_cast<MPow     *>(pi); }
inline MUnary*    munary   (MItem* pi) { return static_cast<MUnary   *>(pi); }
inline MBinary*   mbinary  (MItem* pi) { return static_cast<MBinary  *>(pi); }
inline MNary*     mnary    (MItem* pi) { return static_cast<MNary    *>(pi); }
inline MFunc1D*   mfnc1d   (MItem* pi) { return static_cast<MFunc1D  *>(pi); }
inline MFunc2D*   mfnc2d   (MItem* pi) { return static_cast<MFunc2D  *>(pi); }
inline MFuncND*   mfncnd   (MItem* pi) { return static_cast<MFuncND  *>(pi); }
inline MSFuncND*  msfncnd  (MItem* pi) { return static_cast<MSFuncND *>(pi); }
inline MSequence* msequence(MItem* pi) { return static_cast<MSequence*>(pi); }
inline MEquation* mequation(MItem* pi) { return static_cast<MEquation*>(pi); }
inline MEquality* mequality(MItem* pi) { return static_cast<MEquality*>(pi); }

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

	MITEM(MItem* pi) { m_pi = pi; } // Note that the item is not copied in this constructor
	MITEM(const MITEM& it) { m_pi = (it.m_pi?it.m_pi->copy():0); }
	MITEM(MVariable& v) { m_pi = new MVarRef(&v); }
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

	MItem* copy() { return m_pi->copy(); }

	MITEM Left () { return (static_cast<MBinary*>(m_pi))->LeftItem()->copy(); }
	MITEM Right() { return (static_cast<MBinary*>(m_pi))->RightItem()->copy(); }

	MITEM Item() { return (static_cast<MUnary*>(m_pi))->Item()->copy(); }

	MITEM Param() { return (static_cast<MFunc1D*>(m_pi))->Item()->copy(); }

	Item_Type Type() { return m_pi->Type(); }

	double value() { return (Type()==MNEG? -mnumber(mneg(m_pi)->Item())->value() : mnumber(m_pi)->value()); }

private:
	MItem*	m_pi;
};

//-----------------------------------------------------------------------------
inline MNumber*   mnumber  (MITEM& i) { return static_cast<MNumber  *>(i.ItemPtr()); }
inline MConstant* mconst   (MITEM& i) { return static_cast<MConstant*>(i.ItemPtr()); }
inline MFraction* mfrac    (MITEM& i) { return static_cast<MFraction*>(i.ItemPtr()); }
inline MNamedCt*  mnamed   (MITEM& i) { return static_cast<MNamedCt *>(i.ItemPtr()); }
inline MVarRef*   mvar     (MITEM& i) { return static_cast<MVarRef  *>(i.ItemPtr()); }
inline MBoolean*  mboolean (MITEM& i) { return static_cast<MBoolean* >(i.ItemPtr()); }
inline MNeg*      mneg     (MITEM& i) { return static_cast<MNeg     *>(i.ItemPtr()); }
inline MAdd*      madd     (MITEM& i) { return static_cast<MAdd     *>(i.ItemPtr()); }
inline MSub*      msub     (MITEM& i) { return static_cast<MSub     *>(i.ItemPtr()); }
inline MMul*      mmul     (MITEM& i) { return static_cast<MMul     *>(i.ItemPtr()); }
inline MDiv*      mdiv     (MITEM& i) { return static_cast<MDiv     *>(i.ItemPtr()); }
inline MPow*      mpow     (MITEM& i) { return static_cast<MPow     *>(i.ItemPtr()); }
inline MUnary*    munary   (MITEM& i) { return static_cast<MUnary   *>(i.ItemPtr()); }
inline MBinary*   mbinary  (MITEM& i) { return static_cast<MBinary  *>(i.ItemPtr()); }
inline MNary*     mnary    (MITEM& i) { return static_cast<MNary    *>(i.ItemPtr()); }
inline MFunc1D*   mfnc1d   (MITEM& i) { return static_cast<MFunc1D  *>(i.ItemPtr()); }
inline MFunc2D*   mfnc2d   (MITEM& i) { return static_cast<MFunc2D  *>(i.ItemPtr()); }
inline MSFuncND*  msfncnd  (MITEM& i) { return static_cast<MSFuncND *>(i.ItemPtr()); }
inline MSequence* msequence(MITEM& i) { return static_cast<MSequence*>(i.ItemPtr()); }
inline MEquation* mequation(MITEM& i) { return static_cast<MEquation*>(i.ItemPtr()); }
inline MEquality* mequality(MITEM& i) { return static_cast<MEquality*>(i.ItemPtr()); }

//-----------------------------------------------------------------------------
// Some helper functions to determine the type of an MItem
inline bool is_number(MItem* pi) { return (dynamic_cast<MNumber  *>(pi) != 0); }
inline bool is_unary (MItem* pi) { return (dynamic_cast<MUnary   *>(pi) != 0); }
inline bool is_binary(MItem* pi) { return (dynamic_cast<MBinary  *>(pi) != 0); }
inline bool is_nary  (MItem* pi) { return (dynamic_cast<MNary    *>(pi) != 0); }

inline bool isConst    (MItem* pi) { return (pi->Type() == MCONST   ); }
inline bool is_frac    (MItem* pi) { return (pi->Type() == MFRAC    ); }
inline bool is_named   (MItem* pi) { return (pi->Type() == MNAMED   ); }
inline bool is_var     (MItem* pi) { return (pi->Type() == MVAR     ); }
inline bool is_boolean (MItem* pi) { return (pi->Type() == MBOOL    ); }
inline bool is_neg     (MItem* pi) { return (pi->Type() == MNEG     ); }
inline bool is_add     (MItem* pi) { return (pi->Type() == MADD     ); }
inline bool is_sub     (MItem* pi) { return (pi->Type() == MSUB     ); }
inline bool is_mul     (MItem* pi) { return (pi->Type() == MMUL     ); }
inline bool is_div     (MItem* pi) { return (pi->Type() == MDIV     ); }
inline bool is_pow     (MItem* pi) { return (pi->Type() == MPOW     ); }
inline bool is_equation(MItem* pi) { return (pi->Type() == MEQUATION); }
inline bool is_equality(MItem* pi) { return (pi->Type() == MEQUALITY); }
inline bool is_func1d  (MItem* pi) { return (pi->Type() == MF1D     ); }
inline bool is_matrix  (MItem* pi) { return (pi->Type() == MMATRIX  ); }
inline bool is_sfunc   (MItem* pi) { return (pi->Type() == MSFNC    ); }

inline bool isConst    (MITEM& i) { return (i.Type() == MCONST   ); }
inline bool is_frac    (MITEM& i) { return (i.Type() == MFRAC    ); }
inline bool is_named   (MITEM& i) { return (i.Type() == MNAMED   ); }
inline bool is_var     (MITEM& i) { return (i.Type() == MVAR     ); }
inline bool is_boolean (MITEM& i) { return (i.Type() == MBOOL    ); }
inline bool is_neg     (MITEM& i) { return (i.Type() == MNEG     ); }
inline bool is_add     (MITEM& i) { return (i.Type() == MADD     ); }
inline bool is_sub     (MITEM& i) { return (i.Type() == MSUB     ); }
inline bool is_mul     (MITEM& i) { return (i.Type() == MMUL     ); }
inline bool is_div     (MITEM& i) { return (i.Type() == MDIV     ); }
inline bool is_pow     (MITEM& i) { return (i.Type() == MPOW     ); }
inline bool is_equation(MITEM& i) { return (i.Type() == MEQUATION); }
inline bool is_equality(MITEM& i) { return (i.Type() == MEQUALITY); }
inline bool is_func1d  (MITEM& i) { return (i.Type() == MF1D     ); }
inline bool is_matrix  (MITEM& i) { return (i.Type() == MMATRIX  ); }
inline bool is_sequence(MITEM& i) { return (i.Type() == MSEQUENCE); }
inline bool is_sfunc   (MITEM& i) { return (i.Type() == MSFNC    ); }

inline bool is_number(MITEM& i) { return is_number(i.ItemPtr()); }

inline bool is_func(MItem* pi) { Item_Type t = pi->Type(); return ((t==MF1D) || (t==MF2D) || (t==MFND)); }
inline bool is_int (double  r) { return (r - floor(r) == 0.0); }
inline bool is_abs (MItem* pi) { return (is_func1d(pi) && (mfnc1d(pi)->Name().compare("abs") == 0)); }

inline bool is_func(MITEM& i) { return is_func(i.ItemPtr()); }
inline bool is_abs (MITEM& i) { return is_abs (i.ItemPtr()); }

inline bool is_int(MItem* pi) { return isConst(pi) && is_int(mnumber(pi)->value()); }
inline bool is_int (MITEM& i) { return is_int(i.ItemPtr()); }

inline bool is_fnc1d(MITEM& i, const char* sz) { return (is_func1d(i) && (strcmp(mfnc1d(i)->Name().c_str(), sz) == 0)); }

bool is_rconst(MItem* pi); // is recursive constant
bool is_dependent(MItem* pi, MVariable& x); // is i dependent on x
bool is_dependent(MItem* pi, MItem* px); // is i dependent on x
bool is_scalar(MItem* pi);	// is the expression scalar

bool is_pi(MItem* pi); // is pi the constant pi?
inline bool is_pi(MITEM& i) { return is_pi(i.ItemPtr()); }

inline bool is_zero(MItem* pi) { if (isConst(pi)) return (mconst(pi)->value() == 0.0); else return false; }

inline bool is_scalar(MITEM& m) { return ::is_scalar(m.ItemPtr()); }

//-----------------------------------------------------------------------------
// Algebraic operators for MITEM's
MITEM operator - (MITEM& l);
MITEM operator + (MITEM& l, MITEM& r);
MITEM operator + (MITEM& l, double r);
MITEM operator + (double l, MITEM& r);
MITEM operator - (MITEM& l, MITEM& r);
MITEM operator - (double l, MITEM& r);
MITEM operator - (MITEM& l, double r);
MITEM operator * (MITEM& l, MITEM& r);
MITEM operator * (double l, MITEM& r);
MITEM operator / (MITEM& l, MITEM& r);
MITEM operator / (double l, MITEM& r);
MITEM operator / (MITEM& l, double r);
MITEM operator ^ (MITEM& l, MITEM& r);
MITEM operator ^ (MITEM& l, double r);

//-----------------------------------------------------------------------------
// Fractions
MITEM Fraction(double n, double d);
inline MITEM Fraction(FRACTION f) { return Fraction(f.n, f.d); }

//-----------------------------------------------------------------------------
// Symbolic math functions
MITEM Abs(MITEM& l);
MITEM Sgn(MITEM& l);
MITEM Sin(MITEM& l);
MITEM Cos(MITEM& l);
MITEM Tan(MITEM& l);
MITEM Cot(MITEM& l);
MITEM Sec(MITEM& l);
MITEM Csc(MITEM& l);
MITEM Atan(MITEM& l);
MITEM Cosh(MITEM& l);
MITEM Sinh(MITEM& l);
MITEM Exp(MITEM& l);
MITEM Log(MITEM& l);
MITEM Log10(MITEM& l);
MITEM Float(MITEM& l);
MITEM Sqrt(MITEM& l);
MITEM Erf(MITEM& l);
MITEM Erfc(MITEM& l);

// Chebyshev polynomials
MITEM Tn(int n, MITEM& r);
MITEM Tn(MITEM& l, MITEM& r);

// Bessel functions (1st kind)
MITEM J0(MITEM& l);
MITEM J1(MITEM& l);
MITEM Jn(int n, MITEM& r);
MITEM Jn(MITEM& l, MITEM& r);

// Bessel functions (2e kind)
MITEM Y0(MITEM& l);
MITEM Y1(MITEM& l);
MITEM Yn(int n, MITEM& r);
MITEM Yn(MITEM& l, MITEM& r);

// Miscellaneous functions
MITEM Fac(MITEM& l);
MITEM Binomial(MITEM& l, MITEM& r);
MITEM Identity(MITEM& n);

//-----------------------------------------------------------------------------
// check equality of two items
inline bool operator == (MITEM& l, const std::string& r) { if (is_var(l)) { const std::string& s = mvar(l)->Name(); return (s.compare(r) == 0); } else return false; }
inline bool operator == (MITEM& l, MVariable& x) { if (is_var(l)) { const std::string& s = mvar(l)->Name(); return (s.compare(x.Name()) == 0); } else return false; }
inline bool operator == (MITEM& l, double r) { if (isConst(l)) return (mnumber(l)->value() == r); else return false; }
inline bool operator != (MITEM& l, double r) { return !(l==r); }
inline bool is_dependent(MITEM& l, MVariable& x) { return is_dependent(l.ItemPtr(), x); }
inline bool is_dependent(MITEM& l, MITEM& x) { return is_dependent(l.ItemPtr(), x.ItemPtr()); }

bool is_equal(MItem* pl, MItem* pr);
inline bool is_equal(MItem* pl, MVariable& x) { if (is_var(pl)) { const std::string& s = mvar(pl)->Name(); return (s.compare(x.Name()) == 0); } else return false; }
inline bool operator == (MITEM& l, MITEM& r) { return is_equal(l.ItemPtr(), r.ItemPtr()); }
inline bool operator != (MITEM& l, MITEM& r) { return (is_equal(l.ItemPtr(), r.ItemPtr()) == false); }

//-----------------------------------------------------------------------------
// functions for calculating and comparing heuristics
int op_count(MItem* pi);
inline int op_count(MITEM& i) { return op_count(i.ItemPtr()); }

//-----------------------------------------------------------------------------
bool is_format(MItem* pe, const char* sz);
inline bool is_format(MITEM& e, const char* sz) { return is_format(e.ItemPtr(), sz); }
