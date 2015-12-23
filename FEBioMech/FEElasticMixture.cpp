#include "stdafx.h"
#include "FEElasticMixture.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
FEElasticMixtureMaterialPoint::FEElasticMixtureMaterialPoint() : FEMaterialPoint(new FEElasticMaterialPoint)
{ 
}

//-----------------------------------------------------------------------------
void FEElasticMixtureMaterialPoint::AddMaterialPoint(FEMaterialPoint* pt)
{
	m_mp.push_back(pt);
	pt->SetPrev(this);
	m_w.push_back(1.0);
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticMixtureMaterialPoint::Copy()
{
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint;
	pt->m_w = m_w;
	pt->m_mp = m_mp;
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEElasticMixtureMaterialPoint::Init(bool bflag)
{
	if (bflag)
	{
		for (int i=0; i<(int) m_w.size(); ++i) m_w[i] = 1.0;
	}

	// don't forget to initialize the base class
    FEMaterialPoint::Init(bflag);

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

	if (m_pNext) m_pNext->ShallowCopy(dmp, bsave);
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
	for (int i=0; i<(int)m_mp.size(); ++i) m_mp[i]->Serialize(ar);
    
	if (m_pNext) m_pNext->Serialize(ar);
}

//=============================================================================
//								FEElasticMixture
//=============================================================================

//-----------------------------------------------------------------------------
FEElasticMixture::FEElasticMixture(FEModel* pfem) : FEElasticMaterial(pfem)
{
	AddProperty(&m_pMat, "solid");
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticMixture::CreateMaterialPointData() 
{ 
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint();
	int NMAT = Materials();
	for (int i=0; i<NMAT; ++i) pt->AddMaterialPoint(m_pMat[i]->CreateMaterialPointData());
	return pt;
}

//-----------------------------------------------------------------------------
void FEElasticMixture::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
    FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEElasticMixtureMaterialPoint& mmp = *(mp.ExtractData<FEElasticMixtureMaterialPoint>());
    
	// check the local coordinate systems for each component
	for (int j=0; j<Materials(); ++j)
	{
		FEElasticMaterial* pmj = GetMaterial(j)->GetElasticMaterial();
		FEMaterialPoint& mpj = *mmp.GetPointData(j);
        FEElasticMaterialPoint& pj = *(mpj.ExtractData<FEElasticMaterialPoint>());
		pj.m_Q = pt.m_Q;
		pmj->SetLocalCoordinateSystem(el, n, mpj);
	}
}

//-----------------------------------------------------------------------------
void FEElasticMixture::AddMaterial(FEElasticMaterial* pm) 
{ 
	m_pMat.SetProperty(pm); 
}

//-----------------------------------------------------------------------------
void FEElasticMixture::Init()
{
    for (int i=0; i < (int)m_pMat.size(); ++i) m_pMat[i]->Init();
    
    FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! This function evaluates the stress at the material point by evaluating the
//! individual stress components. 
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
        // but don't copy m_Q since correct value was set in SetLocalCoordinateSystem
		FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;

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
        // but don't copy m_Q since correct value was set in SetLocalCoordinateSystem
		FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
        
		c += m_pMat[i]->Tangent(*pt.m_mp[i])*w[i];
	}

	return c;
}

//-----------------------------------------------------------------------------
//! This function evaluates the stress at the material point by evaluating the
//! individual stress components.
double FEElasticMixture::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());
    
	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate strain energy density
	double sed = 0.0;
	for (int i=0; i < (int) m_pMat.size(); ++i)
	{
		// copy the elastic material point data to the components
        // but don't copy m_Q since correct value was set in SetLocalCoordinateSystem
		FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
        
		sed += m_pMat[i]->StrainEnergyDensity(*pt.m_mp[i])*w[i];
	}
    
	return sed;
}
