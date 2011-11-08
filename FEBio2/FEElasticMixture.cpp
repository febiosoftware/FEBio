#include "stdafx.h"
#include "FECore/tens4d.h"
#include "FEElasticMixture.h"

//-----------------------------------------------------------------------------
// register the material with the framework
REGISTER_MATERIAL(FEElasticMixture, "solid mixture");

//-----------------------------------------------------------------------------
void FEElasticMixture::Init()
{
	FEElasticMaterial::Init();
	for (int i=0; i < (int)m_pMat.size(); ++i)
		m_pMat[i]->Init();
}

//-----------------------------------------------------------------------------
mat3ds FEElasticMixture::Stress(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	mat3ds s;
	
	// calculate stress
	s.zero();
	for (int i=0; i < (int) m_pMat.size(); ++i)
		s += m_pMat[i]->Stress(mp)*w[i];

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEElasticMixture::Tangent(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	tens4ds c(0.);

	// calculate elasticity tensor
	for (int i=0; i < (int) m_pMat.size(); ++i)
		c += m_pMat[i]->Tangent(mp)*w[i];

	return c;
}
