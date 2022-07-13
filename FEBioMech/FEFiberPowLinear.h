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
#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power law - linear

class FEFiberPowLinear : public FEFiberMaterial
{
public:
    FEFiberPowLinear(FEModel* pfem);
    
    //! Cauchy stress
    mat3ds FiberStress(FEMaterialPoint& mp, const vec3d& a0) override;
    
    // Spatial tangent
    tens4ds FiberTangent(FEMaterialPoint& mp, const vec3d& a0) override;
    
    //! Strain energy density
    double FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
    
public:
    FEParamDouble   m_E;		// fiber modulus
    FEParamDouble   m_lam0;     // stretch ratio at end of toe region
    FEParamDouble   m_beta;     // power law exponent in toe region
	double	m_epsf;
};

class FEElasticFiberPowLinear : public FEElasticFiberMaterial_T<FEFiberPowLinear>
{
public:
	FEElasticFiberPowLinear(FEModel* fem) : FEElasticFiberMaterial_T<FEFiberPowLinear>(fem) {}
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Exponential-Power law toe region - linear
class FEFiberExpPowLinear : public FEFiberMaterial
{
public:
	FEFiberExpPowLinear(FEModel* pfem);

	//! Initialization
	bool Validate() override;

	//! Cauchy stress
	mat3ds FiberStress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds FiberTangent(FEMaterialPoint& mp, const vec3d& a0) override;

	//! Strain energy density
	double FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

public:
    FEParamDouble   m_E;		// fiber modulus
    FEParamDouble   m_lam0;     // stretch ratio at end of toe region
    FEParamDouble   m_alpha;    // exponential coefficient
    FEParamDouble   m_beta;     // power law exponent in toe region
	double	m_epsf;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

class FEElasticFiberExpPowLinear : public FEElasticFiberMaterial_T<FEFiberExpPowLinear>
{
public:
	FEElasticFiberExpPowLinear(FEModel* fem) : FEElasticFiberMaterial_T<FEFiberExpPowLinear>(fem) {}
	DECLARE_FECORE_CLASS();
};
