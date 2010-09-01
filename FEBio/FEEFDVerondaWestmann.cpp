#include "stdafx.h"
#include "FEEFDVerondaWestmann.h"

// register the material with the framework
REGISTER_MATERIAL(FEEFDVerondaWestmann, "EFD Veronda-Westmann");

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDVerondaWestmann, FEUncoupledMaterial)
	ADD_PARAMETER(m_VW.m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_VW.m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETERV(m_EFD.m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_EFD.m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FEEFDVerondaWestmann
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void FEEFDVerondaWestmann::Init()
{
	FEUncoupledMaterial::Init();
	
	m_VW.Init();
	// ellipsoidal fiber distribution is stable when combined with a ground matrix
	m_EFD.m_unstable = false;
	m_EFD.Init();
}

//-----------------------------------------------------------------------------
mat3ds FEEFDVerondaWestmann::DevStress(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	mat3ds s = m_VW.DevStress(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	// TODO: is this correct? Should we not add a deviatoric fiber stress?
	s += m_EFD.Stress(mp);
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEEFDVerondaWestmann::DevTangent(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	tens4ds c = m_VW.DevTangent(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	// TODO: is this correct? Should we not add a deviatoric fiber tangent
	c += m_EFD.Tangent(mp);
	
	return c;
}

//-----------------------------------------------------------------------------
double FEEFDVerondaWestmann::BulkModulus()
{
	// Evaluate bulk modulus in reference configuration
	return m_VW.BulkModulus() + m_EFD.BulkModulus();
}
