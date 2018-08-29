//
//  FECrossFluid.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 4/1/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FECrossFluid.h"
#include "FEFluid.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FECrossFluid, FEViscousFluid)
ADD_PARAMETER2(m_mu0, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0");
ADD_PARAMETER2(m_mui, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mui");
ADD_PARAMETER2(m_lam, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda");
ADD_PARAMETER2(m_m  , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(2.0), "m");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FECrossFluid::FECrossFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_mu0 = 0;
    m_mui = 0;
    m_lam = 0;
    m_m = 2;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FECrossFluid::Stress(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double mu = ShearViscosity(pt);
    
    mat3ds s = D*(2*mu);
    
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FECrossFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FECrossFluid::Tangent_RateOfDeformation(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D*D).tr());
    double lamg = m_lam*gdot;
    
    double mu = ShearViscosity(pt);
    double dmu = -2*(m_mu0 - m_mui)*m_m*pow(m_lam,m_m)*pow(gdot,m_m-2)/pow(1+pow(lamg,m_m),2);
    mat3dd I(1.0);
    tens4ds c = dyad1s(D)*(2*dmu) + dyad4s(I)*(2*mu);
    return c;
}

//-----------------------------------------------------------------------------
//! dynamic viscosity
double FECrossFluid::ShearViscosity(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D*D).tr());
    double lamg = m_lam*gdot;
    double mu = m_mui + (m_mu0 - m_mui)/(1+pow(lamg, m_m));
    return mu;
}

//-----------------------------------------------------------------------------
//! bulk viscosity
double FECrossFluid::BulkViscosity(FEMaterialPoint& pt)
{
    return 2*ShearViscosity(pt)/3.;
}
