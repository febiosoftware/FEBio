#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
FEElasticFiberMaterialUC::FEElasticFiberMaterialUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
}

mat3ds FEElasticFiberMaterialUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	vec3d a0 = pt.m_Q.col(0);
	return DevStress(mp, a0);
}

tens4ds FEElasticFiberMaterialUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	vec3d a0 = pt.m_Q.col(0);
	return DevTangent(mp, a0);
}

double FEElasticFiberMaterialUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	vec3d a0 = pt.m_Q.col(0);
	return DevStrainEnergyDensity(mp, a0);
}
