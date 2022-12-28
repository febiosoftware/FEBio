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
#include "FENewtonianViscousSolidUC.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianViscousSolidUC, FEUncoupledMaterial)
	ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(      0.0), "kappa")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(      0.0), "mu"   )->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_secant_tangent, "secant_tangent");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENewtonianViscousSolidUC::FENewtonianViscousSolidUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
    m_kappa = 0.0;
    m_mu = 0.0;
    m_secant_tangent = false;
}

//-----------------------------------------------------------------------------
mat3ds FENewtonianViscousSolidUC::DevStress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    
    mat3ds D = pe.RateOfDeformation();
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds s = I*(D.tr()*(m_kappa - 2*m_mu/3)) + D*(2*m_mu);
    
    // determinant of deformation gradient
    double J = pe.m_J;
    
    return s.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FENewtonianViscousSolidUC::DevTangent(FEMaterialPoint& mp)
{
    mat3dd I(1);
    tens4ds Cv;

    double dt = CurrentTimeIncrement();
    if (dt > 0)
        Cv = (dyad1s(I, I)*(m_kappa - 2 * m_mu / 3) + dyad4s(I, I)*(2 * m_mu)) / (2 * dt);
    else
        Cv.zero();
    
    return Cv;
}

//-----------------------------------------------------------------------------
double FENewtonianViscousSolidUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    return 0;
}

