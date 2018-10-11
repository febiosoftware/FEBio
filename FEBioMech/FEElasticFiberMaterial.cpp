#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterialPoint.h"

BEGIN_FECORE_CLASS(FEElasticFiberMaterial, FEElasticMaterial)
	ADD_PROPERTY(m_fiberGenerator, "fiber");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEElasticFiberMaterial::FEElasticFiberMaterial(FEModel* pfem) : FEElasticMaterial(pfem), m_fiberGenerator(nullptr)
{
}

//-----------------------------------------------------------------------------
bool FEElasticFiberMaterial::Init()
{
	if (m_fiberGenerator == nullptr) return false;
	return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
vec3d FEElasticFiberMaterial::GetFiberVector(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	return pt.m_Q*m_fiberGenerator->GetVector(mp);
}
