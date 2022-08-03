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
#include "FEElasticMaterial.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! Polynomial Hyperelastic
class FEPolynomialHyperElastic : public FEElasticMaterial
{
public:
	FEPolynomialHyperElastic(FEModel* pfem);

public:
	FEParamDouble	m_c[3][3];
	double			m_D1, m_D2;

public:
	//! calculate deviatoric stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp) override;

private:
	double U(double J);
	double UJ(double J);
	double UJJ(double J);

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
