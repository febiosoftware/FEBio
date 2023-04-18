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
#include "FECarreauYasudaFluid.h"
#include "FEFluid.h"

// define the material parameters
BEGIN_FECORE_CLASS(FECarreauYasudaFluid, FEViscousFluid)
	ADD_PARAMETER(m_mu0, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0")->setUnits("P.t")->setLongName("zero shear rate viscosity");
	ADD_PARAMETER(m_mui, FE_RANGE_GREATER_OR_EQUAL(0.0), "mui")->setUnits("P.t")->setLongName("infinite shear rate viscosity");
	ADD_PARAMETER(m_lam, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda")->setUnits(UNIT_TIME)->setLongName("relaxation time");
	ADD_PARAMETER(m_n  , FE_RANGE_GREATER_OR_EQUAL(0.0), "n")->setLongName("power index");
	ADD_PARAMETER(m_a  , FE_RANGE_GREATER_OR_EQUAL(2.0), "a")->setLongName("power denominator");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FECarreauYasudaFluid::FECarreauYasudaFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_mu0 = 0;
    m_mui = 0;
    m_lam = 0;
    m_n = 1;
    m_a = 2;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FECarreauYasudaFluid::Stress(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double mu = ShearViscosity(pt);
    
    mat3ds s = D*(2*mu);
    
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FECarreauYasudaFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FECarreauYasudaFluid::Tangent_RateOfDeformation(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D.sqr()).tr());
    double lamga = pow(m_lam*gdot,m_a);
    
    double mu = m_mui + (m_mu0 - m_mui)*pow(1+lamga, (m_n-1)/m_a);
    double dmu = 2*(m_mu0 - m_mui)*(m_n-1)*pow(m_lam,m_a)*pow(gdot,m_a-2)
    *pow(1+lamga, (m_n-m_a-1)/m_a);
    mat3dd I(1.0);
    tens4ds c = dyad1s(D)*(2*dmu) + dyad4s(I)*(2*mu);
    return c;
}

//-----------------------------------------------------------------------------
//! dynamic viscosity
double FECarreauYasudaFluid::ShearViscosity(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D.sqr()).tr());
    double lamga = pow(m_lam*gdot,m_a);
    double mu = m_mui + (m_mu0 - m_mui)*pow(1+lamga, (m_n-1)/m_a);
    return mu;
}

//-----------------------------------------------------------------------------
//! bulk viscosity
double FECarreauYasudaFluid::BulkViscosity(FEMaterialPoint& pt)
{
    return 2*ShearViscosity(pt)/3.;
}
