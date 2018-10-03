#include "stdafx.h"
#include "FECentrifugalBodyForce.h"
#include "FEElasticMaterial.h"

BEGIN_FECORE_CLASS(FECentrifugalBodyForce, FEBodyForce);
	ADD_PARAMETER(w, "angular_speed");
	ADD_PARAMETER(n, "rotation_axis");
	ADD_PARAMETER(c, "rotation_center");
END_FECORE_CLASS();

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
