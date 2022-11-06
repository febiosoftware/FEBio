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

//-----------------------------------------------------------------------------
// This class implements a version of orthotropic elasticity that remains
// valid for large deformations
class FEOrthoElastic :	public FEElasticMaterial
{
public:
	FEParamDouble	m_E1, m_E2, m_E3;		// Young's moduli
	FEParamDouble	m_v12, m_v23, m_v31;	// Poisson's ratio
	FEParamDouble	m_G12, m_G23, m_G31;	// Shear moduli

public:
	FEOrthoElastic(FEModel* pfem);

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization
	bool Validate() override;

private:
	// Helper function for evaluating Lame coefficients at material point
	bool EvaluateLameCoefficients(FEMaterialPoint& pt, double lam[3][3], double mu[3]);

	// declare parameter list
	DECLARE_FECORE_CLASS();
};
