#include "stdafx.h"
#include "FECore/tens4d.h"
#include "FEElasticMultigeneration.h"

// define the material parameters
// BEGIN_PARAMETER_LIST(FEElasticMultigeneration, FEElasticMaterial)
// END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// Mixture of elastic solids
//////////////////////////////////////////////////////////////////////

void FEElasticMultigeneration::Init()
{
	FEElasticMaterial::Init();
	for (int i=0; i<(int)m_pMat.size(); i++)
		m_pMat[i]->Init();
}

mat3ds FEElasticMultigeneration::Stress(FEMaterialPoint& mp)
{
	FEMultigenerationMaterialPoint& pt = *mp.ExtractData<FEMultigenerationMaterialPoint>();
	FEElasticMaterialPoint& mpt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds s;
	
	// calculate stress
	s.zero();
	
	s = m_pMat[0]->Stress(mp);
	
	// extract deformation gradient
	mat3d Fs = mpt.F;
	double Js = mpt.J;
	
	for (int i=0; i < (int)pt.Fi.size(); ++i) 
	{
		// evaluate deformation gradient for this generation
		mat3d Fi = pt.Fi[i];
       	mpt.F = Fs*Fi;
		double Ji = pt.Ji[i];
		mpt.J = Js*Ji;
		// evaluate stress for this generation
		s += Ji*m_pMat[i+1]->Stress(mp);
	}

	// restore the material point deformation gradient
	mpt.F = Fs;
	mpt.J = Js;
	
	return s;
}

tens4ds FEElasticMultigeneration::Tangent(FEMaterialPoint& mp)
{
	FEMultigenerationMaterialPoint& pt = *mp.ExtractData<FEMultigenerationMaterialPoint>();
	FEElasticMaterialPoint& mpt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	tens4ds c(0.);
	
	c = m_pMat[0]->Tangent(mp);

	// extract deformation gradient
	mat3d Fs = mpt.F;
	double Js = mpt.J;
	
	for (int i=0; i < (int)pt.Fi.size(); ++i) 
	{
		// evaluate deformation gradient for this generation
		mat3d Fi = pt.Fi[i];
       	mpt.F = Fs*Fi;
		double Ji = pt.Ji[i];
		mpt.J = Js*Ji;
		// evaluate stress for this generation
		c += Ji*m_pMat[i+1]->Tangent(mp);
	}
	
	// restore the material point deformation gradient
	mpt.F = Fs;
	mpt.J = Js;
	
	return c;
}

bool FEElasticMultigeneration::HasGeneration(const int igen)
{
	for (int i=0; i<(int)m_pMat.size(); ++i)
		if (m_pMat[i]->GetID() == igen) return true;

	return false;
}

