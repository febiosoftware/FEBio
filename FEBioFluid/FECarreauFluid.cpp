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
ADD_PARAMETER2(m_mu0, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0");
ADD_PARAMETER2(m_mui, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mui");
ADD_PARAMETER2(m_lam, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda");
ADD_PARAMETER2(m_n  , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "n");
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
    double I2 = (D*D).tr();
    
    double mu = m_mui + (m_mu0 - m_mui)*pow(1+m_lam*m_lam*I2*I2, (m_n-1)*0.5);
    
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
    double I2 = (D*D).tr();
    double lam2 = m_lam*m_lam;
    
    double mu = m_mui + (m_mu0 - m_mui)*pow(1+lam2*I2*I2, (m_n-1)*0.5);
    double dmu = (m_mu0 - m_mui)*(m_n-1)*lam2*I2*pow(1+lam2*I2*I2, (m_n-3)*0.5);
    mat3dd I(1.0);
    tens4ds c = dyad1s(D)*(2*dmu) + dyad4s(I)*(2*mu);
    return c;
}
