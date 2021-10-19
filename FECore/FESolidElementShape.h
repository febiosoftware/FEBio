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
#include "FEElementShape.h"

//=============================================================================
// Base class for defining element shape classes for (3D) solid elements
class FESolidElementShape : public FEElementShape
{
public:
	FESolidElementShape(FE_Element_Shape shape, int nodes) : FEElementShape(shape, nodes) {}

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
	FETet4() : FESolidElementShape(ET_TET4, 4) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
class FETet5 : public FESolidElementShape
{
public:
	FETet5() : FESolidElementShape(ET_TET5, 5) {}

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
	FEHex8() : FESolidElementShape(ET_HEX8, 8) {}

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
	FEPenta6() : FESolidElementShape(ET_PENTA6, 6) {}

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
	FEPenta15() : FESolidElementShape(ET_PENTA15, 15) {}

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
	FETet10() : FESolidElementShape(ET_TET10, 10) {}

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
	FETet15() : FESolidElementShape(ET_TET15, 15) {}

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
	FETet20() : FESolidElementShape(ET_TET20, 20) {}

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
	FEHex20() : FESolidElementShape(ET_HEX20, 20) {}

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
	FEHex27() : FESolidElementShape(ET_HEX27, 27) {}

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
	FEPyra5() : FESolidElementShape(ET_PYRA5, 5) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
class FEPyra13 : public FESolidElementShape
{
public:
    FEPyra13() : FESolidElementShape(ET_PYRA13, 13) {}
    
    //! values of shape functions
    void shape_fnc(double* H, double r, double s, double t);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);
    
    //! values of shape function second derivatives
    void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

