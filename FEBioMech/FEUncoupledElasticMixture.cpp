#include "stdafx.h"
#include "FEUncoupledElasticMixture.h"

// define the material parameters
// BEGIN_PARAMETER_LIST(FEUncoupledElasticMixture, FEUncoupledMaterial)
// END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// Mixture of uncoupled elastic solids
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
FEMaterialPoint* FEUncoupledElasticMixture::CreateMaterialPointData() 
{ 
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint();
	int NMAT = Materials();
	pt->m_w.resize(NMAT);
	pt->m_mp.resize(NMAT);
	for (int i=0; i<NMAT; ++i) pt->m_mp[i] = m_pMat[i]->CreateMaterialPointData();
	return pt;
}

//-----------------------------------------------------------------------------
void FEUncoupledElasticMixture::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
    FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
    
	// check the local coordinate systems for each component
	for (int j=0; j<Materials(); ++j)
	{
		FEElasticMaterial* pmj = GetMaterial(j)->GetElasticMaterial();
		FEMaterialPoint& mpj = *mp.GetPointData(j);
        FEElasticMaterialPoint& pj = *(mpj.ExtractData<FEElasticMaterialPoint>());
        pj.m_Q = pt.m_Q;    // copy mixture material's coordinate system into component
		pmj->SetLocalCoordinateSystem(el, n, mpj);
	}
}

//-----------------------------------------------------------------------------
void FEUncoupledElasticMixture::Init()
{
	FEUncoupledMaterial::Init();
	m_K = 0.0;
	for (int i=0; i < (int)m_pMat.size(); ++i) {
		m_pMat[i]->Init();
		m_K += m_pMat[i]->m_K;	// Sum up all the values of the bulk moduli
	}
}

//-----------------------------------------------------------------------------
void FEUncoupledElasticMixture::AddMaterial(FEUncoupledMaterial* pm) 
{ 
	m_pMat.push_back(pm); 
}

//-----------------------------------------------------------------------------
//! Find the index of a material property
int FEUncoupledElasticMixture::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "solid") == 0) return (int) m_pMat.size();
	return -1;
}

//-----------------------------------------------------------------------------
//! Set a material property
bool FEUncoupledElasticMixture::SetProperty(int n, FECoreBase* pm)
{
	assert(n <= (int) m_pMat.size());
	FEUncoupledMaterial* pme = dynamic_cast<FEUncoupledMaterial*>(pm);
	if (pme == 0) return false;
	AddMaterial(pme);
	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEUncoupledElasticMixture::DevStress(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate stress
	mat3ds s(0.0);
	for (int i=0; i < (int)m_pMat.size(); ++i)
	{
		// copy the elastic material point data to the components
        // but don't copy m_Q since correct value was set in SetLocalCoordinateSystem
		FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
		epi.m_s = m_pMat[i]->DevStress(*pt.m_mp[i])*w[i];
		s += epi.m_s;
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEUncoupledElasticMixture::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate elasticity tensor
	tens4ds c(0.);
	for (int i=0; i < (int)m_pMat.size(); ++i)
	{
		// copy the elastic material point data to the components
        // but don't copy m_Q since correct value was set in SetLocalCoordinateSystem
		FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
		c += m_pMat[i]->DevTangent(*pt.m_mp[i])*w[i];
	}
	
	return c;
}

//-----------------------------------------------------------------------------
//! For elastic mixtures, the parameter name is defined as follows:
//!		material.param
//! where material refers to the name of one of the mixture components and
//! param is the parameter name.
//!
FEParam* FEUncoupledElasticMixture::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEUncoupledMaterial::GetParameter(s);

	int NMAT = Materials();
	for (int i=0; i<NMAT; ++i) 
	{
		FEUncoupledMaterial* pmi = GetMaterial(i);
		if (s == pmi->GetName()) return pmi->GetParameter(s.next());
	}

	// no match found
	return 0;
}
