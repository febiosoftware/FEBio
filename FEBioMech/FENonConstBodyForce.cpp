#include "stdafx.h"
#include "FENonConstBodyForce.h"
#include "FEElasticMaterial.h"

BEGIN_PARAMETER_LIST(FENonConstBodyForce, FEBodyForce);
	ADD_PARAMETER(m_val[0], FE_PARAM_DOUBLE_MAPPED, "x");
	ADD_PARAMETER(m_val[1], FE_PARAM_DOUBLE_MAPPED, "y");
	ADD_PARAMETER(m_val[2], FE_PARAM_DOUBLE_MAPPED, "z");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FENonConstBodyForce::FENonConstBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
}

//-----------------------------------------------------------------------------
vec3d FENonConstBodyForce::force(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate the force
	double f[3] = { 0 };
	for (int i = 0; i<3; ++i)
	{
		f[i] = m_val[i].eval(mp);;
	}

	return vec3d(f[0], f[1], f[2]);
}

//-----------------------------------------------------------------------------
mat3ds FENonConstBodyForce::stiffness(FEMaterialPoint& pt)
{
	return mat3ds(0, 0, 0, 0, 0, 0);
}
