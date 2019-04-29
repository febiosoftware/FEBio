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
#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power-law linear response (uncoupled)

class FEFiberPowLinearUncoupled : public FEElasticFiberMaterialUC
{
public:
    FEFiberPowLinearUncoupled(FEModel* pfem);

	// validation
	bool Validate() override;
    
    //! Cauchy stress
    virtual mat3ds DevStress(FEMaterialPoint& mp, const vec3d& n0) override;
    
    // Spatial tangent
    virtual tens4ds DevTangent(FEMaterialPoint& mp, const vec3d& n0) override;
    
    //! Strain energy density
    virtual double DevStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0) override;
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
	double  m_beta;     // power law exponent in toe region

private:
    double  m_ksi;      // power law coefficient in toe region
	double  m_I0;       // m_lam0^2
    double  m_b;        // coefficient in linear region

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
