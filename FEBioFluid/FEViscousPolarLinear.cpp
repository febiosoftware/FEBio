/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include "FEViscousPolarLinear.h"
#include "FEFluidMaterialPoint.h"
#include "FEPolarFluidMaterialPoint.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEViscousPolarLinear, FEViscousPolarFluid)
    ADD_PARAMETER(m_tau, FE_RANGE_GREATER_OR_EQUAL(0.0), "tau");
    ADD_PARAMETER(m_alpha, "alpha");
    ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(0.0), "beta");
    ADD_PARAMETER(m_gamma, FE_RANGE_GREATER_OR_EQUAL(0.0), "gamma");
END_FECORE_CLASS();


//! constructor
FEViscousPolarLinear::FEViscousPolarLinear(FEModel* pfem) : FEViscousPolarFluid(pfem)
{
    m_tau = m_alpha = m_beta = m_gamma = 0;
}

//! dual vector of non-symmetric part of viscous stress in polar fluid
vec3d FEViscousPolarLinear::SkewStressDualVector(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
    FEPolarFluidMaterialPoint* pf = (mp.ExtractData<FEPolarFluidMaterialPoint>());
    
    vec3d h = pf ? pf->m_gf - pt->Vorticity()/2 : - pt->Vorticity()/2;
    
    return h*(-2*m_tau);
}

//! non-symmetric part of viscous stress in polar fluid
mat3da FEViscousPolarLinear::SkewStress(FEMaterialPoint& mp)
{
    return mat3da(SkewStressDualVector(mp));
}

//! tangent of stress with respect to strain J
mat3da FEViscousPolarLinear::SkewTangent_Strain(FEMaterialPoint& mp)
{
    return mat3da(0,0,0);
}

//! tangent of stress with respect to relative rate of rotation tensor H
mat3d FEViscousPolarLinear::SkewTangent_RateOfRotation(FEMaterialPoint& mp)
{
    return mat3dd(-2*m_tau);
}

//! viscous couple stress in polar fluid
mat3d FEViscousPolarLinear::CoupleStress(FEMaterialPoint& mp)
{
    FEPolarFluidMaterialPoint* pf = (mp.ExtractData<FEPolarFluidMaterialPoint>());
    mat3d Psi = pf->m_Psi;
    mat3d M = mat3dd(Psi.trace()*m_alpha) + Psi*(m_beta+m_gamma) + Psi.transpose()*(m_beta - m_gamma);
    
    return M;
}

//! tangent of viscous couple stress with respect to strain J
mat3d FEViscousPolarLinear::CoupleTangent_Strain(FEMaterialPoint& mp)
{
    return mat3d(0.);
}

//! tangent of viscous couple stress with respect to rate of rotation tensor g (left-minor transpose of...)
tens4d FEViscousPolarLinear::CoupleTangent_RateOfRotation(FEMaterialPoint& mp)
{
    mat3dd I(1);
    tens4d m = dyad1(I,I)*m_alpha + dyad2(I,I)*(m_beta-m_gamma) + dyad3(I,I)*(m_beta+m_gamma);
    
    return m;
}

//! dynamic viscosity
double FEViscousPolarLinear::RelativeRotationalViscosity(FEMaterialPoint& mp)
{
    return m_tau;
}

