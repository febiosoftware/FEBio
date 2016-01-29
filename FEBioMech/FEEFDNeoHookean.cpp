#include "stdafx.h"
#include "FEEFDNeoHookean.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_NH.m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_NH.m_v, FE_PARAM_DOUBLE, "v");
	ADD_PARAMETERV(m_EFD.m_beta, FE_PARAM_DOUBLE, 3, "beta");
	ADD_PARAMETERV(m_EFD.m_ksi , FE_PARAM_DOUBLE, 3, "ksi" );
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FEEFDNeoHookean
//////////////////////////////////////////////////////////////////////

bool FEEFDNeoHookean::Init()
{
	if (FEElasticMaterial::Init() == false) return false;
	if (m_NH.Init()  == false) return false;
	if (m_EFD.Init() == false) return false;
	return true;
}

void FEEFDNeoHookean::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	m_NH.Serialize(ar);
	m_EFD.Serialize(ar);
}

mat3ds FEEFDNeoHookean::Stress(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	mat3ds s = m_NH.Stress(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	s += m_EFD.Stress(mp);
	
	return s;
}

tens4ds FEEFDNeoHookean::Tangent(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	tens4ds c = m_NH.Tangent(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	c += m_EFD.Tangent(mp);
	
	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEEFDNeoHookean::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	double sed = m_NH.StrainEnergyDensity(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	sed += m_EFD.StrainEnergyDensity(mp);
    
    return sed;
}

//////////////////////////////////////////////////////////////////////
// FEEFDNeoHookeanOld
//////////////////////////////////////////////////////////////////////

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDNeoHookeanOld, FEElasticMaterial)
	ADD_PARAMETER(m_NH.m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_NH.m_v, FE_PARAM_DOUBLE, "v");
	ADD_PARAMETERV(m_EFD.m_beta, FE_PARAM_DOUBLE, 3, "beta");
	ADD_PARAMETERV(m_EFD.m_ksi , FE_PARAM_DOUBLE, 3, "ksi" );
END_PARAMETER_LIST();

bool FEEFDNeoHookeanOld::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	if (m_NH.Init() == false) return false;
	if (m_EFD.Init() == false) return false;
	return true;
}

void FEEFDNeoHookeanOld::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	m_NH.Serialize(ar);
	m_EFD.Serialize(ar);
}

mat3ds FEEFDNeoHookeanOld::Stress(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	mat3ds s = m_NH.Stress(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	s += m_EFD.Stress(mp);
	
	return s;
}

tens4ds FEEFDNeoHookeanOld::Tangent(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	tens4ds c = m_NH.Tangent(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	c += m_EFD.Tangent(mp);
	
	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEEFDNeoHookeanOld::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	double sed = m_NH.StrainEnergyDensity(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	sed += m_EFD.StrainEnergyDensity(mp);
    
    return sed;
}
