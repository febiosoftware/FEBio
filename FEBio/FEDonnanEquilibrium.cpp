/*
 *  FEDonnanEquilibrium.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *
 */
#include "stdafx.h"
#include "FEDonnanEquilibrium.h"

// register the material with the framework
REGISTER_MATERIAL(FEDonnanEquilibrium, "Donnan equilibrium");

// define the material parameters
BEGIN_PARAMETER_LIST(FEDonnanEquilibrium, FEElasticMaterial)
	ADD_PARAMETER(m_phiwr, FE_PARAM_DOUBLE, "phiw0");
	ADD_PARAMETER(m_cFr, FE_PARAM_DOUBLE, "cF0");
	ADD_PARAMETER(m_bosm, FE_PARAM_DOUBLE, "bosm");
	ADD_PARAMETER(m_Rgas, FE_PARAM_DOUBLE, "R");
	ADD_PARAMETER(m_Tabs, FE_PARAM_DOUBLE, "T");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEDonnanEquilibrium
//-----------------------------------------------------------------------------

void FEDonnanEquilibrium::Init()
{
	if (m_unstable) throw MaterialError("This material is unstable (produces infinite swelling) when used alone.  Combine it in a solid mixture with a material that can resist tension.");
	if (m_phiwr < 0 || m_phiwr > 1) throw MaterialError("phiw0 must be between 0. and 1.");
	if (m_Rgas < 0) throw MaterialError("R must be positive.");
	if (m_Tabs < 0) throw MaterialError("T must be positive.");
	if (m_bosm < 0) throw MaterialError("bosm must be positive.");
}

//-----------------------------------------------------------------------------
mat3ds FEDonnanEquilibrium::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.J;
	
	// calculate fixed charge density in current configuration
	double cF = m_phiwr*m_cFr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(sqrt(cF*cF+m_bosm*m_bosm) - m_bosm);
	
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
	double J = pt.J;

	// calculate fixed charge density in current configuration
	double cF = m_phiwr*m_cFr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double tosm = sqrt(cF*cF+m_bosm*m_bosm);	// tissue osmolarity
	double p = m_Rgas*m_Tabs*(tosm - m_bosm);	// osmotic pressure
	
	// calculate derivative of osmotic pressure w.r.t. J
	double bpi = m_Rgas*m_Tabs*J*cF*cF/(J-1+m_phiwr)/tosm;
	
	mat3dd I(1.0);	// Identity
	
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds c = bpi*IxI + p*(2.0*I4 - IxI);
	return c;
}

//-----------------------------------------------------------------------------
double FEDonnanEquilibrium::BulkModulus()
{
	// calculate osmotic pressure (assume J=1)
	double tosm = sqrt(m_cFr*m_cFr+m_bosm*m_bosm);	// tissue osmolarity
	double pi = m_Rgas*m_Tabs*(tosm - m_bosm);	// osmotic pressure
	
	// calculate derivative of osmotic pressure w.r.t. J
	double bpi = m_Rgas*m_Tabs*m_cFr*m_cFr/m_phiwr/tosm;
	
	return -(pi/3.0 + bpi);
}
