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
#include "FEBinghamFluid.h"
#include "FEFluid.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEBinghamFluid, FEViscousFluid)
    ADD_PARAMETER(m_mu  , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"  )->setLongName("shear viscosity");
    ADD_PARAMETER(m_tauy, FE_RANGE_GREATER_OR_EQUAL(0.0), "tauy")->setLongName("yield stress");
    ADD_PARAMETER(m_n   , FE_RANGE_GREATER_OR_EQUAL(0.0), "n"   )->setLongName("exponent");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBinghamFluid::FEBinghamFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_mu = 0;
    m_tauy = 0;
    m_n = 1;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FEBinghamFluid::Stress(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double mu = ShearViscosity(pt);
    
    mat3ds s = D*(2*mu);
    
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FEBinghamFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FEBinghamFluid::Tangent_RateOfDeformation(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D.sqr()).tr());
    tens4ds c;
    mat3dd I(1.0);

    if (gdot > 0) {
        double mu = m_mu + m_tauy/gdot*(1-exp(-m_n*gdot));
        double dmu = m_tauy/pow(gdot,2)*((1+m_n*gdot)*exp(-m_n*gdot)-1);
        c = dyad1s(D)*(4/gdot*dmu) + dyad4s(I)*(2*mu);
    }
    else {
        double mu = m_mu + m_tauy*m_n;
        c = dyad4s(I)*(2*mu);
    }
    return c;
}

//-----------------------------------------------------------------------------
//! dynamic viscosity
double FEBinghamFluid::ShearViscosity(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D.sqr()).tr());
    double mu = (gdot > 0) ? m_mu + m_tauy/gdot*(1-exp(-m_n*gdot)) : m_mu + m_tauy*m_n;
    return mu;
}

//-----------------------------------------------------------------------------
//! bulk viscosity
double FEBinghamFluid::BulkViscosity(FEMaterialPoint& pt)
{
    return 2*ShearViscosity(pt)/3.;
}
