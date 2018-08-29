#include "stdafx.h"
#include "FECentrifugalBodyForce.h"
#include "FEElasticMaterial.h"

BEGIN_PARAMETER_LIST(FECentrifugalBodyForce, FEBodyForce);
	ADD_PARAMETER(w, FE_PARAM_DOUBLE, "angular_speed");
	ADD_PARAMETER(n, FE_PARAM_VEC3D, "rotation_axis");
	ADD_PARAMETER(c, FE_PARAM_VEC3D, "rotation_center");
END_PARAMETER_LIST();

FECentrifugalBodyForce::FECentrifugalBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
	w = 0.0;
	n = vec3d(0,0,1);
	c = vec3d(0,0,0);
}

vec3d FECentrifugalBodyForce::force(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3ds K = stiffness(mp);
	return K*(pt.m_rt - c);
}

mat3ds FECentrifugalBodyForce::stiffness(FEMaterialPoint& mp) 
{ 
	return (mat3dd(1) - dyad(n))*(-w*w); 
}
