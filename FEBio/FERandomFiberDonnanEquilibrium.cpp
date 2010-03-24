#include "stdafx.h"
#include "FERandomFiberDonnanEquilibrium.h"

// register the material with the framework
REGISTER_MATERIAL(FERandomFiberDonnanEquilibrium, "EFD Donnan equilibrium");

// define the material parameters
BEGIN_PARAMETER_LIST(FERandomFiberDonnanEquilibrium, FEElasticMaterial)
	ADD_PARAMETER(m_DEQ.m_phiwr, FE_PARAM_DOUBLE, "phiw0");
	ADD_PARAMETER(m_DEQ.m_cFr, FE_PARAM_DOUBLE, "cF0");
	ADD_PARAMETER(m_DEQ.m_bosm, FE_PARAM_DOUBLE, "bosm");
	ADD_PARAMETER(m_DEQ.m_Rgas, FE_PARAM_DOUBLE, "R");
	ADD_PARAMETER(m_DEQ.m_Tabs, FE_PARAM_DOUBLE, "T");
	ADD_PARAMETERV(m_Fib.m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_Fib.m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FERandomFiberDonnanEquilibrium
//-----------------------------------------------------------------------------

FERandomFiberDonnanEquilibrium::FERandomFiberDonnanEquilibrium()
{
}

//-----------------------------------------------------------------------------
void FERandomFiberDonnanEquilibrium::Init()
{
	FEElasticMaterial::Init();
	
	m_DEQ.Init();
	m_Fib.Init();
}

//-----------------------------------------------------------------------------
mat3ds FERandomFiberDonnanEquilibrium::Stress(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	mat3ds s = m_DEQ.Stress(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	s += m_Fib.Stress(mp);
		
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FERandomFiberDonnanEquilibrium::Tangent(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	tens4ds c = m_DEQ.Tangent(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	c += m_Fib.Tangent(mp);

	return c;
}

//-----------------------------------------------------------------------------
double FERandomFiberDonnanEquilibrium::BulkModulus()
{
	// Evaluate bulk modulus in reference configuration
	// Only the matrix contributes to the bulk modulus
	return m_DEQ.BulkModulus();
}
