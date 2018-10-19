#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterialPoint.h"

//-----------------------------------------------------------------------------
FEElasticFiberMaterial::FEElasticFiberMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

mat3ds FEElasticFiberMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	vec3d a0 = pt.m_Q.col(0);
	return Stress(mp, a0);
}

tens4ds FEElasticFiberMaterial::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	vec3d a0 = pt.m_Q.col(0);
	return Tangent(mp, a0);
}

double FEElasticFiberMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	vec3d a0 = pt.m_Q.col(0);
	return StrainEnergyDensity(mp, a0);
}
