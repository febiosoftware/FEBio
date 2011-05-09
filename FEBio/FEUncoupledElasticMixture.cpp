#include "stdafx.h"
#include "FEBioLib/tens4d.h"
#include "FEUncoupledElasticMixture.h"

// register the material with the framework
REGISTER_MATERIAL(FEUncoupledElasticMixture, "uncoupled solid mixture");

// define the material parameters
// BEGIN_PARAMETER_LIST(FEUncoupledElasticMixture, FEUncoupledMaterial)
// END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// Mixture of uncoupled elastic solids
//////////////////////////////////////////////////////////////////////

void FEUncoupledElasticMixture::Init()
{
	FEUncoupledMaterial::Init();
	for (int i=0; i < m_pMat.size(); ++i) {
		m_pMat[i]->Init();
		m_K += m_pMat[i]->m_K;	// Sum up all the values of the bulk moduli
	}
}

mat3ds FEUncoupledElasticMixture::DevStress(FEMaterialPoint& mp)
{
	mat3ds s;
	
	// calculate stress
	s.zero();
	for (int i=0; i < m_pMat.size(); ++i)
		s += m_pMat[i]->DevStress(mp);
	
	return s;
}

tens4ds FEUncoupledElasticMixture::DevTangent(FEMaterialPoint& mp)
{
	tens4ds c(0.);
	
	// calculate elasticity tensor
	for (int i=0; i < m_pMat.size(); ++i)
		c += m_pMat[i]->DevTangent(mp);
	
	return c;
}
