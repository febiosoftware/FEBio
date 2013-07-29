#include "stdafx.h"
#include "FEElasticMixture.h"

//-----------------------------------------------------------------------------
FEElasticMixture::FEElasticMixture()
{

}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticMixture::CreateMaterialPointData() 
{ 
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint();
	pt->m_w.resize(m_pMat.size());
	return pt;
}

//-----------------------------------------------------------------------------
void FEElasticMixture::Init()
{
	FEElasticMaterial::Init();
	for (int i=0; i < (int)m_pMat.size(); ++i)
	{
		m_pMat[i]->m_pParent = m_pParent;
		m_pMat[i]->Init();
	}
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

//-----------------------------------------------------------------------------
//! For elastic mixtures, the parameter name is defined as follows:
//!		material.param
//! where material refers to the name of one of the mixture components and
//! param is the parameter name.
//!
FEParam* FEElasticMixture::GetParameter(const ParamString& s)
{
	// see if this is a composite name
	if (s.count() == 1) return FEElasticMaterial::GetParameter(s);

	// else, find the variable name and search the mixture components
	int NMAT = Materials();
	for (int i=0; i<NMAT; ++i) 
	{
		FEElasticMaterial* pmi = GetMaterial(i);
		if (s == pmi->GetName()) return pmi->GetParameter(s.next());
	}

	// no match found
	return 0;
}
