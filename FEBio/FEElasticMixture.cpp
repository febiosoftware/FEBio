#include "stdafx.h"
#include "FECore/tens4d.h"
#include "FEElasticMixture.h"

// register the material with the framework
REGISTER_MATERIAL(FEElasticMixture, "solid mixture");

// define the material parameters
// BEGIN_PARAMETER_LIST(FEElasticMixture, FEElasticMaterial)
// END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// Mixture of elastic solids
//////////////////////////////////////////////////////////////////////

void FEElasticMixture::Init()
{
	FEElasticMaterial::Init();
	for (int i=0; i < m_pMat.size(); ++i)
		m_pMat[i]->Init();
}

mat3ds FEElasticMixture::Stress(FEMaterialPoint& mp)
{
	mat3ds s;
	
	// calculate stress
	s.zero();
	for (int i=0; i < m_pMat.size(); ++i)
		s += m_pMat[i]->Stress(mp);

	return s;
}

tens4ds FEElasticMixture::Tangent(FEMaterialPoint& mp)
{
	tens4ds c(0.);

	// calculate elasticity tensor
	for (int i=0; i < m_pMat.size(); ++i)
		c += m_pMat[i]->Tangent(mp);

	return c;
}

double FEElasticMixture::BulkModulus()
{
	double k = 0;
	
	// calculate bulk modulus
	for (int i=0; i < m_pMat.size(); ++i)
		k += m_pMat[i]->BulkModulus();
	
	return k;
}

