/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power law - linear

class FEFiberPowLinear : public FEElasticFiberMaterial
{
public:
    FEFiberPowLinear(FEModel* pfem);
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp, const vec3d& a0) override;
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp, const vec3d& a0) override;
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
    double  m_beta;     // power law exponent in toe region
};

//-----------------------------------------------------------------------------
//! Power law toe region - linear
//! TODO: I want to delete one of these materials
class FEFiberPowerLinear : public FEElasticFiberMaterial
{
public:
	FEFiberPowerLinear(FEModel* pfem);

	//! Initialization
	bool Validate() override;

	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp, const vec3d& a0) override;

	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

public:
	double	m_E;		// fiber modulus
	double  m_lam0;     // stretch ratio at end of toe region
	double  m_beta;     // power law exponent in toe region

private:
	double  m_I0;       // m_lam0^2
	double  m_ksi;      // power law coefficient in toe region
	double  m_b;        // coefficient in linear region

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
