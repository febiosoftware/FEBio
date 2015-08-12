//
//  FEOsmoticVirialExpansion.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/30/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//
#include "stdafx.h"
#include "FEOsmoticVirialExpansion.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEOsmoticVirialExpansion, FEElasticMaterial)
    ADD_PARAMETER(m_phiwr, FE_PARAM_DOUBLE, "phiw0");
    ADD_PARAMETER(m_cr, FE_PARAM_DOUBLE, "cr");
    ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
    ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
    ADD_PARAMETER(m_c3, FE_PARAM_DOUBLE, "c3");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEDonnanEquilibrium
void FEOsmoticVirialExpansion::Init()
{
	FEElasticMaterial::Init();
    if (m_phiwr < 0 || m_phiwr > 1) throw MaterialError("phiw0 must be between 0. and 1.");
    if (m_cr < 0) throw MaterialError("cr must be positive.");
}

//-----------------------------------------------------------------------------
mat3ds FEOsmoticVirialExpansion::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // jacobian
    double J = pt.m_J;
    
    // calculate concentration in current configuration
    double c = m_phiwr*m_cr/(J-1+m_phiwr);
    
    // calculate osmotic pressure
    double p = m_c1*c + m_c2*c*c + m_c3*pow(c,3);
    
    // calculate T = -p*I
    mat3dd I(1.0);	// identity tensor
    mat3ds s = -p*I;
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEOsmoticVirialExpansion::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // jacobian
    double J = pt.m_J;
    
    // calculate concentration in current configuration
    double c = m_phiwr*m_cr/(J-1+m_phiwr);
    
    // calculate osmotic pressure
    double p = m_c1*c + m_c2*c*c + m_c3*pow(c,3);
    
    // calculate derivative of osmotic pressure w.r.t. J
    double dpdJ = -(m_c1 + 2*m_c2*c + 3*m_c3*c*c)*c/(J -1 + m_phiwr);
    
    mat3dd I(1.0);	// Identity
    
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    
    // calculate tangent osmotic modulus
    tens4ds C = IxI*(-J*dpdJ) + (2.0*I4 - IxI)*p;
    
    return C;
}
