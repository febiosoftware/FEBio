#include "stdafx.h"
#include "FEPreStrainElastic.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEPreStrainElastic, FEElasticMaterial)
	ADD_PARAMETER(m_Fp, FE_PARAM_MAT3D, "F0");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPreStrainElastic::FEPreStrainElastic(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_pmat = 0;
}

//-----------------------------------------------------------------------------
void FEPreStrainElastic::Init()
{
	if (m_pmat == 0) throw MaterialError("Base material undefined");
	m_pmat->Init();
}

//-----------------------------------------------------------------------------
//! Set the base material
void FEPreStrainElastic::SetBaseMaterial(FEElasticMaterial* pbase) { m_pmat = pbase; }

//-----------------------------------------------------------------------------
//! This material only has one property
int FEPreStrainElastic::Properties()
{
	return 1;
}

//-----------------------------------------------------------------------------
FECoreBase* FEPreStrainElastic::GetProperty(int i)
{
	if (i == 0) return m_pmat;
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEPreStrainElastic::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "elastic") == 0) return 0; else return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEPreStrainElastic::SetProperty(int i, FECoreBase* pm)
{
	if (i==0)
	{
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
		if (pme)
		{ 
			SetBaseMaterial(pme);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEPreStrainElastic::CreateMaterialPointData()
{ 
	return m_pmat->CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
mat3ds FEPreStrainElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());

	// store the original deformation gradient
	mat3d F0 = ep.m_F;
	double J0 = ep.m_J;

	// pre-multiply the pre-strain
	ep.m_F = ep.m_F*m_Fp;
	ep.m_J = ep.m_F.det();

	// evaluate the stress
	mat3ds s = m_pmat->Stress(mp);

	// restore original deformation gradient
	ep.m_F = F0;
	ep.m_J = J0;

	// return stress
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPreStrainElastic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());

	// store the original deformation gradient
	mat3d F0 = ep.m_F;
	double J0 = ep.m_J;

	// pre-multiply the pre-strain
	ep.m_F = ep.m_F*m_Fp;
	ep.m_J = ep.m_F.det();

	// evaluate the tangent
	tens4ds c = m_pmat->Tangent(mp);

	// restore original deformation gradient
	ep.m_F = F0;
	ep.m_J = J0;

	// return spatial tangent
	return c;
}
