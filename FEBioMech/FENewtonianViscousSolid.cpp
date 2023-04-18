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



#include "stdafx.h"
#include "FENewtonianViscousSolid.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianViscousSolid, FEElasticMaterial)
	ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(      0.0), "kappa")->setUnits("P.t")->setLongName("bulk viscosity");
	ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(      0.0), "mu"   )->setUnits("P.t")->setLongName("shear viscosity");
    ADD_PARAMETER(m_secant_tangent, "secant_tangent");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENewtonianViscousSolid::FENewtonianViscousSolid(FEModel* pfem) : FEElasticMaterial(pfem) 
{
    m_kappa = 0.0;
    m_mu = 0.0;
    m_secant_tangent = false;
}

//-----------------------------------------------------------------------------
mat3ds FENewtonianViscousSolid::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    mat3ds D = pt.RateOfDeformation();
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds s = I*(D.tr()*(m_kappa - 2*m_mu/3)) + D*(2*m_mu);
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FENewtonianViscousSolid::Tangent(FEMaterialPoint& mp)
{
    const FETimeInfo& tp = GetTimeInfo();
    tens4ds Cv;
    
    if (tp.timeIncrement > 0) {
        mat3dd I(1);
        double tmp = tp.alphaf*tp.gamma/(tp.beta*tp.timeIncrement);
        Cv = (dyad1s(I)*(m_kappa - 2 * m_mu / 3) + dyad4s(I)*(2 * m_mu))*tmp;
    }
    else Cv.zero();
    
    return Cv;
}

//-----------------------------------------------------------------------------
double FENewtonianViscousSolid::StrainEnergyDensity(FEMaterialPoint& mp)
{
    return 0;
}

