#include "stdafx.h"
#include "FECore/tens4d.h"
#include "FEElasticMultigeneration.h"

//=============================================================================
// define the material parameters
BEGIN_PARAMETER_LIST(FEGenerationMaterial, FEElasticMaterial)
	ADD_PARAMETER(btime, FE_PARAM_DOUBLE, "start_time");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEGenerationMaterial::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "solid") == 0) return 0;
	else return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEGenerationMaterial::SetProperty(int i, FECoreBase* pm)
{
	if (i==0)
	{
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
		if (pme == 0) return false;
		m_pMat = pme;
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEGenerationMaterial::Init()
{
	assert(m_pMat);
	m_pMat->Init();
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

//=============================================================================
FEMaterialPoint* FEMultigenerationMaterialPoint::Copy()
{
	FEMultigenerationMaterialPoint* pt = new FEMultigenerationMaterialPoint(*this);
	pt->m_pmat = m_pmat;
	pt->Fi = Fi;
	pt->Ji = Ji;
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEMultigenerationMaterialPoint::Serialize(DumpFile& ar)
{
	if (m_pt) m_pt->Serialize(ar);
			
	if (ar.IsSaving())
	{
		for (int i=0; i < (int)Fi.size(); i++)
			ar << Fi[i];
		for (int i=0; i < (int)Ji.size(); i++)
			ar << Ji[i];
	}
	else
	{
		for (int i=0; i < (int)Fi.size(); i++)
			ar >> Fi[i];
		for (int i=0; i < (int)Ji.size(); i++)
			ar >> Ji[i];
	}
			
	if (m_pt) m_pt->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEMultigenerationMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (m_pt) m_pt->ShallowCopy(dmp, bsave);
			
	if (bsave)
	{
		for (int i=0; i < (int)Fi.size(); i++)
			dmp << Fi[i];
		for (int i=0; i < (int)Ji.size(); i++)
			dmp << Ji[i];
	}
	else
	{
		for (int i=0; i < (int)Fi.size(); i++)
			dmp >> Fi[i];
		for (int i=0; i < (int)Ji.size(); i++)
			dmp >> Ji[i];
	}
}

//-----------------------------------------------------------------------------
void FEMultigenerationMaterialPoint::Init(bool bflag)
{
	if (m_pt) m_pt->Init(bflag);

	if (bflag)
	{
		Fi.clear();
		Ji.clear();
		m_tgen = 0.0;
	}

	// get the time
	double t = FEMaterialPoint::time;

	// Check if this constitutes a new generation
	int igen = m_pmat->CheckGeneration(t);
	t = m_pmat->m_MG[igen]->btime;
	if ((bflag == false) && (t>m_tgen))
	{
		FEElasticMaterialPoint& pt = *((*this).ExtractData<FEElasticMaterialPoint>());
					
		// push back F and J to define relative deformation gradient of this generation
		mat3d F = pt.m_F;
		double J = pt.m_J; 
		Fi.push_back(F.inverse());
		Ji.push_back(1.0/J);

		m_tgen = t;
	}
}

//=============================================================================

//-----------------------------------------------------------------------------
//! Find the index of a material property
int FEElasticMultigeneration::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "generation") == 0) return (int) m_MG.size();
	return -1;
}

//-----------------------------------------------------------------------------
//! Set a material property
bool FEElasticMultigeneration::SetProperty(int n, FECoreBase* pm)
{
	assert(n == (int)m_MG.size());
	FEGenerationMaterial* pmg = dynamic_cast<FEGenerationMaterial*>(pm);
	if (pmg)
	{
		m_MG.push_back(pmg);
		return true;
	}
	return false;
}

//--------------------------------------------------------------------------------
// Check if time t constitutes a new generation and return that generation
int FEElasticMultigeneration::CheckGeneration(const double t)
{
	int ngen = m_MG.size();
	for (int igen=1; igen<ngen; ++igen)
	{
		if (t < m_MG[igen]->btime) return igen - 1;
	}
	return ngen - 1;
}

//-----------------------------------------------------------------------------
void FEElasticMultigeneration::Init()
{
	FEElasticMaterial::Init();
	for (int i=0; i<(int)m_MG.size(); i++) m_MG[i]->Init();
}

//-----------------------------------------------------------------------------
mat3ds FEElasticMultigeneration::Stress(FEMaterialPoint& mp)
{
	FEMultigenerationMaterialPoint& pt = *mp.ExtractData<FEMultigenerationMaterialPoint>();
	FEElasticMaterialPoint& mpt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds s;
	
	// calculate stress
	s.zero();
	
	s = m_MG[0]->Stress(mp);
	
	// extract deformation gradient
	mat3d Fs = mpt.m_F;
	double Js = mpt.m_J;
	
	for (int i=0; i < (int)pt.Fi.size(); ++i) 
	{
		// evaluate deformation gradient for this generation
		mat3d Fi = pt.Fi[i];
       	mpt.m_F = Fs*Fi;
		double Ji = pt.Ji[i];
		mpt.m_J = Js*Ji;
		// evaluate stress for this generation
		s += Ji*m_MG[i+1]->Stress(mp);
	}

	// restore the material point deformation gradient
	mpt.m_F = Fs;
	mpt.m_J = Js;
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEElasticMultigeneration::Tangent(FEMaterialPoint& mp)
{
	FEMultigenerationMaterialPoint& pt = *mp.ExtractData<FEMultigenerationMaterialPoint>();
	FEElasticMaterialPoint& mpt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	tens4ds c(0.);
	
	c = m_MG[0]->Tangent(mp);

	// extract deformation gradient
	mat3d Fs = mpt.m_F;
	double Js = mpt.m_J;
	
	for (int i=0; i < (int)pt.Fi.size(); ++i) 
	{
		// evaluate deformation gradient for this generation
		mat3d Fi = pt.Fi[i];
       	mpt.m_F = Fs*Fi;
		double Ji = pt.Ji[i];
		mpt.m_J = Js*Ji;
		// evaluate stress for this generation
		c += Ji*m_MG[i+1]->Tangent(mp);
	}
	
	// restore the material point deformation gradient
	mpt.m_F = Fs;
	mpt.m_J = Js;
	
	return c;
}
