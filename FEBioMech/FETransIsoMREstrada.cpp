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
#include "FETransIsoMREstrada.h"

//////////////////////////////////////////////////////////////////////
// FETransIsoMREstrada
//////////////////////////////////////////////////////////////////////

// define the material parameters
BEGIN_FECORE_CLASS(FETransIsoMREstrada, FEUncoupledMaterial)
    ADD_PARAMETER(c1, "c1");
    ADD_PARAMETER(c2, "c2");
    ADD_PARAMETER(m_fib.m_c3, "c3");
    ADD_PARAMETER(m_fib.m_c4, "c4");
    ADD_PARAMETER(m_fib.m_c5, "c5");
    ADD_PARAMETER(m_fib.m_lam1, "lam_max");
    
    ADD_PROPERTY(m_fiber, "fiber");
    ADD_PROPERTY(m_ac, "active_contraction", FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETransIsoMREstrada::FETransIsoMREstrada(FEModel* pfem) : FEUncoupledMaterial(pfem), m_fib(pfem)
{
    m_ac = nullptr;
    m_fib.SetParent(this);
    m_fiber = nullptr;
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPointData* FETransIsoMREstrada::CreateMaterialPointData() 
{
    // create the elastic solid material point
    FEMaterialPointData* ep = new FEElasticMaterialPoint;
    
    // create the material point from the active contraction material
    if (m_ac) ep->Append(m_ac->CreateMaterialPointData());

	return ep;
}

//-----------------------------------------------------------------------------
mat3ds FETransIsoMREstrada::DevStress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    
    // calculate deviatoric left Cauchy-Green tensor
    mat3ds B = pt.DevLeftCauchyGreen();
    
    // calculate square of B
    mat3ds B2 = B.sqr();

    // material axes
    mat3d Q = GetLocalCS(mp);

    // get the fiber vector in local coordinates
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d a0 = Q * fiber;
    
    // Invariants of B (= invariants of C)
    // Note that these are the invariants of Btilde, not of B!
    double I1 = B.tr();
    
    // --- TODO: put strain energy derivatives here ---
    // Wi = dW/dIi
    double W1 = c1(mp);
    double W2 = c2(mp);
    // ------------------------------------------------
    
    // calculate T = F*dW/dC*Ft
    mat3ds T = B*(W1 + W2*I1) - B2*W2;
    
    // calculate stress s = pI + 2/J * dev(T)
    mat3ds s = T.dev()*(2.0/J);
    
    // calculate the passive fiber stress
    mat3ds fs = m_fib.DevFiberStress(mp, a0);
    
    // calculate the active fiber stress (if provided)
    if (m_ac) fs += m_ac->ActiveStress(mp, a0);
    
    return s + fs;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FETransIsoMREstrada::DevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    double Ji = 1.0/J;
    
    // calculate deviatoric left Cauchy-Green tensor: B = F*Ft
    mat3ds B = pt.DevLeftCauchyGreen();
    
    // calculate square of B
    mat3ds B2 = B.sqr();
    
    // Invariants of B (= invariants of C)
    double I1 = B.tr();
    double I2 = 0.5*(I1*I1 - B2.tr());
    
    // --- TODO: put strain energy derivatives here ---
    // Wi = dW/dIi
    double W1, W2;
    W1 = c1(mp);
    W2 = c2(mp);
    // ------------------------------------
    
    // calculate dWdC:C
    double WC = W1*I1 + 2*W2*I2;
    
    // calculate C:d2WdCdC:C
    double CWWC = 2*I2*W2;
    
    mat3dd I(1);    // Identity
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    tens4ds BxB = dyad1s(B);
    tens4ds B4  = dyad4s(B);

    // material axes
    mat3d Q = GetLocalCS(mp);

    // get the fiber vector in local coordinates
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d a0 = Q * fiber;
    
    // deviatoric cauchy-stress, trs = trace[s]/3
    mat3ds devs = pt.m_s.dev();
    
    // d2W/dCdC:C
    mat3ds WCCxC = B*(W2*I1) - B2*W2;
    
    tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);
    tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;
    
    // add the passive fiber stiffness
    c += m_fib.DevFiberTangent(mp, a0);
    
    // add the active fiber stiffness
    if (m_ac) c += m_ac->ActiveStiffness(mp, a0);
    
    return c;
}

//-----------------------------------------------------------------------------
double FETransIsoMREstrada::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // calculate deviatoric left Cauchy-Green tensor
    mat3ds B = pt.DevLeftCauchyGreen();
    
    // calculate square of B
    mat3ds B2 = B.sqr();
    
    // material axes
    mat3d Q = GetLocalCS(mp);

    // get the fiber vector in local coordinates
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d a0 = Q * fiber;

    // Invariants of B (= invariants of C)
    // Note that these are the invariants of Btilde, not of B!
    double I1 = B.tr();
    double I2 = 0.5*(I1*I1 - B2.tr());
    
    // calculate sed
    double sed = c1(mp)*(I1-3) + c2(mp)*(I2-3);
    
    // add the fiber sed
    sed += m_fib.DevFiberStrainEnergyDensity(mp, a0);
    
    return sed;
}

//-----------------------------------------------------------------------------
// update force-velocity material point
void FETransIsoMREstrada::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& timeInfo)
{
    // material axes
    mat3d Q = GetLocalCS(mp);

    // get the fiber vector in local coordinates
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d a0 = Q * fiber;

    if (m_ac) m_ac->UpdateSpecializedMaterialPoints(mp, timeInfo, a0);
}
