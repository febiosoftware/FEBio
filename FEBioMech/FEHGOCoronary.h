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
#include "FEUncoupledMaterial.h"
#include "FEUncoupledFiberExpLinear.h"

//-----------------------------------------------------------------------------
// Constitutive formulation from:
// Holzapfel, et.a., "Determination of layer-specific mechanical properties of human
// coronary arteries with nonatherosclerotic intimal thickening
// and related constitutive modeling", Am J Physiol Heart Circ Physiol 289
class FEHGOCoronary : public FEUncoupledMaterial
{
public:
	FEHGOCoronary(FEModel* pfem);

public:
	double		m_rho;
	double		m_k1;
	double		m_k2;

	FEParamVec3 m_fiber;

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;

protected:

	// declare parameter list
	DECLARE_FECORE_CLASS();
};
