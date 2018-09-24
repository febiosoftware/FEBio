//
//  FEPowellEyringFluid.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 4/1/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEPowellEyringFluid.h"
#include "FEFluid.h"
#include <math.h>

// define the material parameters
BEGIN_PARAMETER_LIST(FEPowellEyringFluid, FEViscousFluid)
	ADD_PARAMETER(m_mu0, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0");
	ADD_PARAMETER(m_mui, FE_RANGE_GREATER_OR_EQUAL(0.0), "mui");
	ADD_PARAMETER(m_lam, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEPowellEyringFluid::FEPowellEyringFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_mu0 = 0;
    m_mui = 0;
    m_lam = 0;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FEPowellEyringFluid::Stress(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double mu = ShearViscosity(pt);
    
    mat3ds s = D*(2*mu);
    
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FEPowellEyringFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FEPowellEyringFluid::Tangent_RateOfDeformation(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D*D).tr());
    double lamg = m_lam*gdot;
    
    double mu = ShearViscosity(pt);
    double dmu = (lamg < 1e-3) ? -2*(m_mu0 - m_mui)*m_lam*m_lam/3. :
    (2*(m_mu0 - m_mui)*(gdot/sqrt(1 + pow(lamg,2)) - asinh(lamg)/m_lam))/pow(gdot,3);
    mat3dd I(1.0);
    tens4ds c = dyad1s(D)*(2*dmu) + dyad4s(I)*(2*mu);
    return c;
}

//-----------------------------------------------------------------------------
//! dynamic viscosity
double FEPowellEyringFluid::ShearViscosity(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = vt.RateOfDeformation();
    double gdot = sqrt(2*(D*D).tr());
    double lamg = m_lam*gdot;
    double mu = (lamg < 1e-3) ? m_mu0 : m_mui + (m_mu0 - m_mui)*asinh(lamg)/lamg;
    return mu;
}

//-----------------------------------------------------------------------------
//! bulk viscosity
double FEPowellEyringFluid::BulkViscosity(FEMaterialPoint& pt)
{
    return 2*ShearViscosity(pt)/3.;
}
