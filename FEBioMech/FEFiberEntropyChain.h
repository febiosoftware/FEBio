/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Exponential-power law
//! (Variation that includes a shear term)
class FEFiberEntropyChain : public FEFiberMaterial
{
public:
	FEFiberEntropyChain(FEModel* pfem);

	//! Initialization
	bool Validate() override;

	//! Cauchy stress
	mat3ds FiberStress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds FiberTangent(FEMaterialPoint& mp, const vec3d& a0) override;

	//! Strain energy density
	double FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

public:
	double          m_N;        // coefficient of micro-combination number
	FEParamDouble	m_ksi;		// measure of fiber modulus which equals to nkT
	int             m_term;     // how many Tayler approximation terms will be used
	FEParamDouble   m_mu;       // shear modulus

	double	m_epsf;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEElasticFiberEntropyChain : public FEElasticFiberMaterial_T<FEFiberEntropyChain>
{
public:
    FEElasticFiberEntropyChain(FEModel* fem) : FEElasticFiberMaterial_T<FEFiberEntropyChain>(fem) {}
    DECLARE_FECORE_CLASS();
};
