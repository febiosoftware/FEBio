#include "stdafx.h"
#include "FERandomFiberVerondaWestmann.h"

// register the material with the framework
REGISTER_MATERIAL(FERandomFiberVerondaWestmann, "random fiber Veronda-Westmann");

// define the material parameters
BEGIN_PARAMETER_LIST(FERandomFiberVerondaWestmann, FEIncompressibleMaterial)
	ADD_PARAMETER(m_VW.m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_VW.m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETERV(m_EFD.m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_EFD.m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FERandomFiberVerondaWestmann
//////////////////////////////////////////////////////////////////////

void FERandomFiberVerondaWestmann::Init()
{
	FEIncompressibleMaterial::Init();
	
	m_VW.Init();
	// ellipsoidal fiber distribution is stable when combined with a ground matrix
	m_EFD.m_unstable = false;
	m_EFD.Init();
}

mat3ds FERandomFiberVerondaWestmann::Stress(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	mat3ds s = m_VW.Stress(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	s += m_EFD.Stress(mp);
	
	return s;
}

tens4ds FERandomFiberVerondaWestmann::Tangent(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	tens4ds c = m_VW.Tangent(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	c += m_EFD.Tangent(mp);
	
	return c;
}

double FERandomFiberVerondaWestmann::BulkModulus()
{
	// Evaluate bulk modulus in reference configuration
	return m_VW.BulkModulus() + m_EFD.BulkModulus();
}
