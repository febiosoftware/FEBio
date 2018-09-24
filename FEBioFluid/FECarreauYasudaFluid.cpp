//
//  FECarreauYasudaFluid.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 4/1/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FECarreauYasudaFluid.h"
#include "FEFluid.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FECarreauYasudaFluid, FEViscousFluid)
	ADD_PARAMETER(m_mu0, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0");
	ADD_PARAMETER(m_mui, FE_RANGE_GREATER_OR_EQUAL(0.0), "mui");
	ADD_PARAMETER(m_lam, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda");
	ADD_PARAMETER(m_n  , FE_RANGE_GREATER_OR_EQUAL(0.0), "n");
	ADD_PARAMETER(m_a  , FE_RANGE_GREATER_OR_EQUAL(2.0), "a");
END_PARAMETER_LIST();

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
    double gdot = sqrt(2*(D*D).tr());
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
    double gdot = sqrt(2*(D*D).tr());
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
