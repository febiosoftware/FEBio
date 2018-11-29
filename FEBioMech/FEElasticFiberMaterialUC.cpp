#include "stdafx.h"
#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
FEElasticFiberMaterialUC::FEElasticFiberMaterialUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
}

mat3ds FEElasticFiberMaterialUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q.col(0);
	return DevStress(mp, a0);
}

tens4ds FEElasticFiberMaterialUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q.col(0);
	return DevTangent(mp, a0);
}

double FEElasticFiberMaterialUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q.col(0);
	return DevStrainEnergyDensity(mp, a0);
}
