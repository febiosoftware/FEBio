#include "stdafx.h"
#include "FEEFDDonnanEquilibrium.h"

// register the material with the framework
REGISTER_MATERIAL(FEEFDDonnanEquilibrium, "EFD Donnan equilibrium");

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDDonnanEquilibrium, FEElasticMaterial)
	ADD_PARAMETER(m_DEQ.m_phiwr, FE_PARAM_DOUBLE, "phiw0");
	ADD_PARAMETER(m_DEQ.m_cFr, FE_PARAM_DOUBLE, "cF0");
	ADD_PARAMETER(m_DEQ.m_bosm, FE_PARAM_DOUBLE, "bosm");
	ADD_PARAMETER(m_DEQ.m_Rgas, FE_PARAM_DOUBLE, "R");
	ADD_PARAMETER(m_DEQ.m_Tabs, FE_PARAM_DOUBLE, "T");
	ADD_PARAMETERV(m_Fib.m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_Fib.m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEEFDDonnanEquilibrium
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void FEEFDDonnanEquilibrium::Init()
{
	FEElasticMaterial::Init();
	// Donnan equilibrium is stable when combined with fibers that can resist swelling
	m_DEQ.m_unstable = false;
	m_DEQ.Init();
	// ellipsoidal fiber distribution is stable when combined with a ground matrix
	m_Fib.m_unstable = false;
	m_Fib.Init();
}

//-----------------------------------------------------------------------------
mat3ds FEEFDDonnanEquilibrium::Stress(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	mat3ds s = m_DEQ.Stress(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	s += m_Fib.Stress(mp);
		
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEEFDDonnanEquilibrium::Tangent(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	tens4ds c = m_DEQ.Tangent(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	c += m_Fib.Tangent(mp);

	return c;
}

//-----------------------------------------------------------------------------
double FEEFDDonnanEquilibrium::BulkModulus()
{
	// Evaluate bulk modulus in reference configuration
	// Only the matrix contributes to the bulk modulus
	return m_DEQ.BulkModulus();
}
