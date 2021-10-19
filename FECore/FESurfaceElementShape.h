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
class FESurfaceElementShape : public FEElementShape
{
public:
	FESurfaceElementShape(FE_Element_Shape shape, int nodes) : FEElementShape(shape, nodes) {}

	//! values of shape functions
	virtual void shape_fnc(double* H, double r, double s) = 0;

	//! values of shape function derivatives
	virtual void shape_deriv(double* Hr, double* Hs, double r, double s) = 0;

	//! values of shape function second derivatives
	virtual void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) = 0;
};

//=============================================================================
// Class for QUAD4 elements
class FEQuad4 : public FESurfaceElementShape
{
public:
	FEQuad4() : FESurfaceElementShape(ET_QUAD4, 4) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for QUAD8 elements
class FEQuad8 : public FESurfaceElementShape
{
public:
	FEQuad8() : FESurfaceElementShape(ET_QUAD8, 8) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for QUAD9 elements
class FEQuad9 : public FESurfaceElementShape
{
public:
	FEQuad9() : FESurfaceElementShape(ET_QUAD9, 9) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for TRI3 elements
class FETri3 : public FESurfaceElementShape
{
public:
	FETri3() : FESurfaceElementShape(ET_TRI3, 4) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for TRI6 elements
class FETri6 : public FESurfaceElementShape
{
public:
	FETri6() : FESurfaceElementShape(ET_TRI6, 6) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for TRI7 elements
class FETri7 : public FESurfaceElementShape
{
public:
	FETri7() : FESurfaceElementShape(ET_TRI7, 7) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for TRI10 elements
class FETri10 : public FESurfaceElementShape
{
public:
	FETri10() : FESurfaceElementShape(ET_TRI10, 10) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};
