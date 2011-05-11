/*
 *  FEPerfectOsmometer.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 7/14/10.
 *
 */
#include "stdafx.h"
#include "FEPerfectOsmometer.h"
#include "fem.h"

// register the material with the framework
REGISTER_MATERIAL(FEPerfectOsmometer, "perfect osmometer");

// define the material parameters
BEGIN_PARAMETER_LIST(FEPerfectOsmometer, FEElasticMaterial)
	ADD_PARAMETER(m_phiwr, FE_PARAM_DOUBLE, "phiw0");
	ADD_PARAMETER(m_iosm, FE_PARAM_DOUBLE, "iosm");
	ADD_PARAMETER(m_bosm, FE_PARAM_DOUBLE, "bosm");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEPerfectOsmometer
//-----------------------------------------------------------------------------

void FEPerfectOsmometer::Init()
{
	if (m_unstable) throw MaterialError("This material is unstable (cannot sustain shear) when used alone.  Combine it in a solid mixture with a material that can resist tension.");
	if (m_phiwr < 0 || m_phiwr > 1) throw MaterialError("phiw0 must be between 0. and 1.");
	if (m_iosm < 0) throw MaterialError("iosm must be positive.");
	if (m_bosm < 0) throw MaterialError("bosm must be positive.");
	
	m_Rgas = FEM::GetGlobalConstant("R");
	m_Tabs = FEM::GetGlobalConstant("T");
	
	if (m_Rgas <= 0) throw MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) throw MaterialError("A positive absolute temperature T must be defined in Globals section");
	
}

//-----------------------------------------------------------------------------
mat3ds FEPerfectOsmometer::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.J;
	
	// calculate internal concentration in current configuration
	double iosm = m_iosm*m_phiwr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(iosm - m_bosm);
	
	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPerfectOsmometer::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.J;

	// calculate internal osmolarity in current configuration
	double iosm = m_iosm*m_phiwr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(iosm - m_bosm);
	
	// calculate derivative of osmotic pressure w.r.t. J
	double dp = -m_Rgas*m_Tabs*iosm/(J-1+m_phiwr);
	
	mat3dd I(1.0);	// Identity
	
	tens4ds I1 = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds c = -J*dp*I1 + p*(2.0*I4 - I1);
	return c;
}

