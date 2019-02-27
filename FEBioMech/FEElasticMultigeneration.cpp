#include "stdafx.h"
#include "FECore/tens4d.h"
#include "FECore/log.h"
#include "FECore/FECoreKernel.h"
#include "FEElasticMultigeneration.h"

//=============================================================================
// define the material parameters
BEGIN_FECORE_CLASS(FEGenerationMaterial, FEElasticMaterial)
	
	// material parameters
	ADD_PARAMETER(btime, "start_time");

	// the solid property
	ADD_PROPERTY(m_pMat, "solid");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEGenerationMaterial::FEGenerationMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_pMat = 0;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEGenerationMaterial::Stress(FEMaterialPoint& pt)
{
	return m_pMat->Stress(pt);
}
		
//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEGenerationMaterial::Tangent(FEMaterialPoint& pt)
{
	return m_pMat->Tangent(pt);
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEGenerationMaterial::StrainEnergyDensity(FEMaterialPoint& pt)
{
	return m_pMat->StrainEnergyDensity(pt);
}

//=============================================================================
FEMultigenerationMaterialPoint::FEMultigenerationMaterialPoint() : FEMaterialPoint(new FEElasticMaterialPoint)
{
    m_tgen = 0.0;
    m_ngen = 1;     // the first generation is always active
}

//-----------------------------------------------------------------------------
void FEMultigenerationMaterialPoint::AddMaterialPoint(FEMaterialPoint* pt)
{
	m_mp.push_back(pt);
	pt->SetPrev(this);
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEMultigenerationMaterialPoint::Copy()
{
	FEMultigenerationMaterialPoint* pt = new FEMultigenerationMaterialPoint(*this);
    pt->m_mp = m_mp;
	pt->m_pmat = m_pmat;
    pt->m_tgen = m_tgen;
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEMultigenerationMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	if (ar.IsShallow())
	{
		ar & m_tgen & m_ngen;
		for (int i=0; i < (int)m_mp.size(); i++) m_mp[i]->Serialize(ar);
		// TODO: shallow copy m_pmat
	}
	else
	{
		if (ar.IsSaving())
		{
			ar << m_tgen << m_ngen;
			ar << (int)m_mp.size();
			for (int i=0; i < (int)m_mp.size(); i++) m_mp[i]->Serialize(ar);
		}
		else
		{
			ar >> m_tgen >> m_ngen;
			int mp_size;
			ar >> mp_size;
			m_mp.resize(mp_size);
			for (int i=0; i < mp_size; i++)
			{
				m_mp[i] = new FEElasticMaterialPoint;
				m_mp[i]->Serialize(ar);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEMultigenerationMaterialPoint::Init()
{
    FEMaterialPoint::Init();
    for (int i=0; i<(int)m_mp.size(); ++i) m_mp[i]->Init();

	m_tgen = 0.0;
	m_ngen = 1;
}

//-----------------------------------------------------------------------------
void FEMultigenerationMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    FEMaterialPoint::Update(timeInfo);
    for (int i=0; i<(int)m_mp.size(); ++i) m_mp[i]->Update(timeInfo);

	// get the time
	double t = timeInfo.currentTime;

	// Check if this constitutes a new generation
	int igen = m_pmat->CheckGeneration(t);
	t = m_pmat->m_MG[igen]->btime;
	if (t>m_tgen)
	{
		FEElasticMaterialPoint& pt = *((*this).ExtractData<FEElasticMaterialPoint>());
					
		// push back F and J to define relative deformation gradient of this generation
		mat3d F = pt.m_F;
		double J = pt.m_J;
        FEElasticMaterialPoint& pe = *(m_mp[m_ngen]->ExtractData<FEElasticMaterialPoint>());
        pe.m_F = F.inverse();
        pe.m_J = 1.0/J;
		m_tgen = t;
        ++m_ngen;
	}
}

//=============================================================================

BEGIN_FECORE_CLASS(FEElasticMultigeneration, FEElasticMaterial)
	ADD_PROPERTY(m_MG, "generation");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEElasticMultigeneration::FEElasticMultigeneration(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEElasticMultigeneration::CreateMaterialPointData()
{
    // use the zero-th generation material point as the base elastic material point
    FEMultigenerationMaterialPoint* pt = new FEMultigenerationMaterialPoint();
    pt->m_pmat = this;
    int NMAT = Materials();
    for (int i=0; i<NMAT; ++i) pt->AddMaterialPoint(m_MG[i]->CreateMaterialPointData());
    return pt;
}

//--------------------------------------------------------------------------------
// Check if time t constitutes a new generation and return that generation
int FEElasticMultigeneration::CheckGeneration(const double t)
{
	int ngen = (int)m_MG.size();
	for (int igen=1; igen<ngen; ++igen)
	{
		if (t < m_MG[igen]->btime) return igen - 1;
	}
	return ngen - 1;
}

//-----------------------------------------------------------------------------
mat3ds FEElasticMultigeneration::Stress(FEMaterialPoint& mp)
{
	FEMultigenerationMaterialPoint& pt = *mp.ExtractData<FEMultigenerationMaterialPoint>();
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds s;
	
	// calculate stress
	s.zero();
	
	// extract deformation gradient
	mat3d Fs = ep.m_F;
	double Js = ep.m_J;
	
	for (int i=0; i < pt.m_ngen; ++i)
	{
        // evaluate deformation gradient for this generation
        FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
        
        // store safe copies of Fi and Ji for this generation
        mat3d Fi = epi.m_F;
        double Ji = epi.m_J;
        
        // copy the elastic material point data to the components
        epi.m_rt = ep.m_rt;
        epi.m_r0 = ep.m_r0;
        
        // evaluate relative deformation gradient
       	epi.m_F = Fs*Fi;
        epi.m_J = Js*Ji;
        
        // evaluate stress for this generation
        s += epi.m_s = Ji*m_MG[i]->Stress(*pt.m_mp[i]);
        
        // restore the material point deformation gradient
        epi.m_F = Fi;
        epi.m_J = Ji;
	}
    
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEElasticMultigeneration::Tangent(FEMaterialPoint& mp)
{
	FEMultigenerationMaterialPoint& pt = *mp.ExtractData<FEMultigenerationMaterialPoint>();
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
	
	tens4ds c(0.);
	
	// extract deformation gradient
	mat3d Fs = ep.m_F;
	double Js = ep.m_J;
	
	for (int i=0; i < pt.m_ngen; ++i)
	{
		// evaluate deformation gradient for this generation
        FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
        
        // store safe copies of Fi and Ji for this generation
        mat3d Fi = epi.m_F;
        double Ji = epi.m_J;
        
        // copy the elastic material point data to the components
        epi.m_rt = ep.m_rt;
        epi.m_r0 = ep.m_r0;
        
        // evaluate relative deformation gradient
       	epi.m_F = Fs*Fi;
        epi.m_J = Js*Ji;
        
		// evaluate tangent for this generation
		c += Ji*m_MG[i]->Tangent(*pt.m_mp[i]);

        // restore the material point deformation gradient
        epi.m_F = Fi;
        epi.m_J = Ji;
	}
	
	return c;
}

//-----------------------------------------------------------------------------
double FEElasticMultigeneration::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEMultigenerationMaterialPoint& pt = *mp.ExtractData<FEMultigenerationMaterialPoint>();
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
	double sed = 0.0;
	
	// extract deformation gradient
	mat3d Fs = ep.m_F;
	double Js = ep.m_J;
	
	for (int i=0; i < pt.m_ngen; ++i)
	{
        // evaluate deformation gradient for this generation
        FEElasticMaterialPoint& epi = *pt.m_mp[i]->ExtractData<FEElasticMaterialPoint>();
        
        // store safe copies of Fi and Ji for this generation
        mat3d Fi = epi.m_F;
        double Ji = epi.m_J;
        
        // copy the elastic material point data to the components
        epi.m_rt = ep.m_rt;
        epi.m_r0 = ep.m_r0;
        
        // evaluate relative deformation gradient
       	epi.m_F = Fs*Fi;
        epi.m_J = Js*Ji;
        
        // evaluate strain energy density for this generation
//        sed += Ji*m_MG[i]->StrainEnergyDensity(*pt.m_mp[i]);
        double dsed = Ji*m_MG[i]->StrainEnergyDensity(*pt.m_mp[i]);
        sed += dsed;
        
        // restore the material point deformation gradient
        epi.m_F = Fi;
        epi.m_J = Ji;
	}
     
	return sed;
}
