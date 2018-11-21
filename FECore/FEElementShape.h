#pragma once
#include "fecore_enum.h"

//=============================================================================
// Base class for defining element shape classes.
// The element shape classes are used for evaluating shape functions and their derivatives.
class FEElementShape
{
public:
	FEElementShape(FE_Element_Shape eshape) : m_shape(eshape) {}
	virtual ~FEElementShape() {}

	FE_Element_Shape shape() const { return m_shape; }

private:
	FE_Element_Shape m_shape;
};

//=============================================================================
// Base class for defining element shape classes for (3D) solid elements
class FESolidElementShape : public FEElementShape
{
public:
	FESolidElementShape(FE_Element_Shape shape) : FEElementShape(shape) {}

	//! values of shape functions
	virtual void shape_fnc(double* H, double r, double s, double t) = 0;

	//! values of shape function derivatives
	virtual void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) = 0;

	//! values of shape function second derivatives
	virtual void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) = 0;
};

//=============================================================================
class FETet4 : public FESolidElementShape
{
public:
	FETet4() : FESolidElementShape(ET_TET4) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
class FEHex8 : public FESolidElementShape
{
public:
	FEHex8() : FESolidElementShape(ET_HEX8) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
class FEPenta6 : public FESolidElementShape
{
public:
	FEPenta6() : FESolidElementShape(ET_PENTA6) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
class FEPenta15 : public FESolidElementShape
{
public:
	FEPenta15(): FESolidElementShape(ET_PENTA15) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
class FETet10 : public FESolidElementShape
{
public:
	FETet10() : FESolidElementShape(ET_TET10) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
class FETet15 : public FESolidElementShape
{
public:
	FETet15() : FESolidElementShape(ET_TET15) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
class FETet20 : public FESolidElementShape
{
public:
	FETet20() : FESolidElementShape(ET_TET20) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
class FEHex20 : public FESolidElementShape
{
public:
	FEHex20() : FESolidElementShape(ET_HEX20) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
//! Base class for 27-node quadratic hexahedral element
class FEHex27 : public FESolidElementShape
{
public:
	FEHex27() : FESolidElementShape(ET_HEX27) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
class FEPyra5 : public FESolidElementShape
{
public:
	FEPyra5() : FESolidElementShape(ET_PYRA5) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};
