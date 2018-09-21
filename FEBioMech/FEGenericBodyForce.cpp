#include "stdafx.h"
#include "FEGenericBodyForce.h"
#include "FEElasticMaterial.h"

BEGIN_PARAMETER_LIST(FEGenericBodyForce, FEBodyForce);
	ADD_PARAMETER(m_force, FE_PARAM_VEC3D_MAPPED, "force");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEGenericBodyForce::FEGenericBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
}

//-----------------------------------------------------------------------------
vec3d FEGenericBodyForce::force(FEMaterialPoint &mp)
{
	return m_force.eval(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEGenericBodyForce::stiffness(FEMaterialPoint& pt)
{
	return mat3ds(0, 0, 0, 0, 0, 0);
}
