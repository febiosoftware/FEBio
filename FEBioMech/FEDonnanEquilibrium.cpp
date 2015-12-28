/*
 *  FEDonnanEquilibrium.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *
 */
#include "stdafx.h"
#include "FEDonnanEquilibrium.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEDonnanEquilibrium, FEElasticMaterial)
	ADD_PARAMETER2(m_phiwr, FE_PARAM_DOUBLE, FE_RANGE_LEFT_OPEN(0.0, 1.0), "phiw0");
    ADD_PARAMETER(m_phisr, FE_PARAM_DOUBLE, "phis0");
	ADD_PARAMETER(m_cFr, FE_PARAM_DOUBLE, "cF0");
	ADD_PARAMETER2(m_bosm, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "bosm");
    ADD_PARAMETER2(m_Phi, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "Phi");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEDonnanEquilibrium
bool FEDonnanEquilibrium::Init()
{
    if (!m_binit) {
        if (m_phisr >= 0) {
            m_bnew = true;
            m_phiwr = 1 - m_phisr;  // use value at t=0 to initialize
        }
        m_binit = true;
    }

	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");
	
	if (m_Rgas <= 0) return MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) return MaterialError("A positive absolute temperature T must be defined in Globals section");
	
	return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
mat3ds FEDonnanEquilibrium::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// calculate fixed charge density in current configuration
    double cF;
    if (m_bnew)
        cF = m_phiwr*m_cFr/(J-m_phisr);
    else
        cF = m_phiwr*m_cFr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*m_Phi*(sqrt(cF*cF+m_bosm*m_bosm) - m_bosm);
	
	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEDonnanEquilibrium::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;

	// calculate fixed charge density in current configuration
    double cF;
    if (m_bnew)
        cF = m_phiwr*m_cFr/(J-m_phisr);
    else
        cF = m_phiwr*m_cFr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double tosm = sqrt(cF*cF+m_bosm*m_bosm);	// tissue osmolarity
	double p = m_Rgas*m_Tabs*m_Phi*(tosm - m_bosm);	// osmotic pressure
	
	// calculate derivative of osmotic pressure w.r.t. J
    double bpi;
    if (m_bnew)
        bpi = m_Rgas*m_Tabs*m_Phi*J*cF*cF/(J-m_phisr)/tosm;
    else
        bpi = m_Rgas*m_Tabs*m_Phi*J*cF*cF/(J-1+m_phiwr)/tosm;
	
	mat3dd I(1.0);	// Identity
	
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds c = bpi*IxI + p*(2.0*I4 - IxI);
	return c;
}
