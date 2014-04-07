#include "stdafx.h"
#include "FEElasticMixture.h"

//-----------------------------------------------------------------------------
FEElasticMixtureMaterialPoint::FEElasticMixtureMaterialPoint()
{ 
	m_pt = new FEElasticMaterialPoint; 
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticMixtureMaterialPoint::Copy()
{
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint;
	pt->m_w = m_w;
	pt->m_mp = m_mp;
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEElasticMixtureMaterialPoint::Init(bool bflag)
{
	if (bflag)
	{
		for (int i=0; i<(int) m_w.size(); ++i) m_w[i] = 1.0;
	}

	m_pt->Init(bflag);
	for (int i=0; i<(int)m_mp.size(); ++i) m_mp[i]->Init(bflag);
}

//-----------------------------------------------------------------------------
void FEElasticMixtureMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_w;
	}
	else
	{
		dmp >> m_w;
	}
	for (int i=0; i<(int)m_mp.size(); ++i) m_mp[i]->ShallowCopy(dmp, bsave);

	if (m_pt) m_pt->ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FEElasticMixtureMaterialPoint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_w;
	}
	else
	{
		ar >> m_w;
	}
}

//=============================================================================
//								FEElasticMixture
//=============================================================================

//-----------------------------------------------------------------------------
FEElasticMixture::FEElasticMixture(FEModel* pfem) : FEElasticMaterial(pfem)
{

}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticMixture::CreateMaterialPointData() 
{ 
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint();
	int NMAT = Materials();
	pt->m_w.resize(NMAT);
	pt->m_mp.resize(NMAT);
	for (int i=0; i<NMAT; ++i) pt->m_mp[i] = m_pMat[i]->CreateMaterialPointData();
	return pt;
}

//-----------------------------------------------------------------------------
//! \todo why does the base gets this material's parent?
void FEElasticMixture::Init()
{
	FEElasticMaterial::Init();
	for (int i=0; i < (int)m_pMat.size(); ++i)
	{
		m_pMat[i]->SetParent(GetParent());
		m_pMat[i]->Init();
	}
}

//-----------------------------------------------------------------------------
void FEElasticMixture::AddMaterial(FEElasticMaterial* pm) 
{ 
	m_pMat.push_back(pm); 
}

//-----------------------------------------------------------------------------
//! Find the index of a material property
int FEElasticMixture::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "solid") == 0) return (int) m_pMat.size();
	return -1;
}

//-----------------------------------------------------------------------------
//! Set a material property
bool FEElasticMixture::SetProperty(int n, FECoreBase* pm)
{
	assert(n <= (int) m_pMat.size());
	FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
	if (pme == 0) return false;
	AddMaterial(pme);
	return true;
}

//-----------------------------------------------------------------------------
//! This function evaluates the stress at the material point by evaluating the
//! individual stress components. 

//! \todo This function copies some material point data from the mixture material
//!       to the component materials, but not all. I need to check if more data needs
//!       to be copied.
mat3ds FEElasticMixture::Stress(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate stress
	mat3ds s(0.0);
	for (int i=0; i < (int) m_pMat.size(); ++i)
	{
		// copy the elastic material point data to the components
		FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
		epi.m_Q = ep.m_Q;
		s += epi.m_s = m_pMat[i]->Stress(*pt.m_mp[i])*w[i];
	}

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEElasticMixture::Tangent(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate elasticity tensor
	tens4ds c(0.);
	for (int i=0; i < (int) m_pMat.size(); ++i)
	{
		// copy the elastic material point data to the components
		FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
		epi.m_Q = ep.m_Q;
		c += m_pMat[i]->Tangent(*pt.m_mp[i])*w[i];
	}

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
