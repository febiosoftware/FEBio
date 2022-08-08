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
#include "FEHolmesMowUC.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEHolmesMowUC, FEUncoupledMaterial)
    ADD_PARAMETER(m_mu, FE_RANGE_GREATER(0.0), "mu")->setLongName("shear modulus");
    ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "beta")->setLongName("power exponent");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FEHolmesMowUC::DevStress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // calculate left Cauchy-Green tensor
    mat3ds bt = pt.DevLeftCauchyGreen();
    
    // calculate invariants of B
    double I1 = bt.tr();
    
    // Exponential term
    double eQ = exp(m_b*(I1-3));
    
    // calculate stress
    mat3ds st = bt*(m_mu*eQ/pt.m_J);
    
    return st.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEHolmesMowUC::DevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // calculate left Cauchy-Green tensor
    mat3ds bt = pt.DevLeftCauchyGreen();
    
    // calculate invariants of B
    double I1 = bt.tr();
    
    // Exponential term
    double eQ = exp(m_b*(I1-3));
    
    // calculate stress
    mat3ds st = bt*(m_mu*eQ/pt.m_J);
    
    // calculate identity tensor
    mat3dd I(1);

    // calculate elasticity tensor
    tens4ds ct = dyad1s(bt)*(2*m_b*m_mu*eQ/pt.m_J);
    
    // This is the final value of the elasticity tensor
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    
    ct += - 1./3.*(ddots(ct,IxI) - IxI*(ct.tr()/3.))
    + 2./3.*((I4-IxI/3.)*st.tr()-dyad1s(st.dev(),I));
    
    return ct;
}

//-----------------------------------------------------------------------------
double FEHolmesMowUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // calculate left Cauchy-Green tensor
    mat3ds bt = pt.DevLeftCauchyGreen();
    
    // calculate invariants of B
    double I1 = bt.tr();
    
    // Exponential term
    double eQ = exp(m_b*(I1-3));
    
    // calculate strain energy density
    double sed = m_mu/(2*m_b)*(eQ-1);
    
    return sed;
}
