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

//-----------------------------------------------------------------------------
//! For elastic mixtures, the parameter name is defined as follows:
//!		material.param
//! where material refers to the name of one of the mixture components and
//! param is the parameter name.
//!
FEParam* FEElasticMixture::GetParameter(const char* sz)
{
	// see if this is a composite name
	char* ch = strchr((char*)sz, '.');

	// if not, find the parameter in the base class
	if (ch == 0) return FEElasticMaterial::GetParameter(sz);

	// else, find the variable name and search the mixture components
	*ch = 0;
	const char* szvar2 = ch+1;
	int NMAT = Materials();
	for (int i=0; i<NMAT; ++i) 
	{
		FEElasticMaterial* pmi = GetMaterial(i);
		if (strcmp(sz, pmi->GetName()) == 0) return pmi->GetParameter(szvar2);
	}

	// no match found
	return 0;
}
