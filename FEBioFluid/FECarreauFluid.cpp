//
//  FECarreauFluid.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/29/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#include "FECarreauFluid.h"
#include "FEFluid.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FECarreauFluid, FEViscousFluid)
	ADD_PARAMETER(m_mu0, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0");
	ADD_PARAMETER(m_mui, FE_RANGE_GREATER_OR_EQUAL(0.0), "mui");
	ADD_PARAMETER(m_lam, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda");
	ADD_PARAMETER(m_n  , FE_RANGE_GREATER_OR_EQUAL(0.0), "n");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FECarreauFluid::FECarreauFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_mu0 = 0;
    m_mui = 0;
    m_lam = 0;
    m_n = 1;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FECarreauFluid::Stress(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double mu = ShearViscosity(pt);
    
    mat3ds s = D*(2*mu);
    
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FECarreauFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FECarreauFluid::Tangent_RateOfDeformation(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D*D).tr());
    double lamg2 = m_lam*m_lam*gdot*gdot;
    
    double mu = m_mui + (m_mu0 - m_mui)*pow(1+lamg2, (m_n-1)*0.5);
    double dmu = 2*(m_mu0 - m_mui)*(m_n-1)*m_lam*m_lam*pow(1+lamg2, (m_n-3)*0.5);
    mat3dd I(1.0);
    tens4ds c = dyad1s(D)*(2*dmu) + dyad4s(I)*(2*mu);
    return c;
}

//-----------------------------------------------------------------------------
//! dynamic viscosity
double FECarreauFluid::ShearViscosity(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D*D).tr());
    double mu = m_mui + (m_mu0 - m_mui)*pow(1+m_lam*m_lam*gdot*gdot, (m_n-1)*0.5);
    return mu;
}

//-----------------------------------------------------------------------------
//! bulk viscosity
double FECarreauFluid::BulkViscosity(FEMaterialPoint& pt)
{
    return 2*ShearViscosity(pt)/3.;
}
