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
#include "FEMembraneMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a wrinkle model for an Ogden material
//
class FEWrinkleOgdenMaterial : public FEMembraneMaterial
{
public:
	FEWrinkleOgdenMaterial(FEModel* pfem);

public: // material parameters
	double	m_u;
	double	m_a;
	double	m_l0[2];	// pre-stretch
	bool	m_bwrinkle;	// wrinkle flag

public: // material interface

	//! calculate stress at material point
	void Stress(FEMaterialPoint& pt, double s[3]) override;

	//! calculate tangent stiffness at material point
	void Tangent(FEMaterialPoint& pt, double D[3][3]) override;

protected:
	void principals(FEMaterialPoint& pt, double l[2], double v[4]);

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
