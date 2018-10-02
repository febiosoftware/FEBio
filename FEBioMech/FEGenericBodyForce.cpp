#include "stdafx.h"
#include "FEGenericBodyForce.h"
#include "FEElasticMaterial.h"

BEGIN_PARAMETER_LIST(FEGenericBodyForce, FEBodyForce);
	ADD_PARAMETER(m_force, "force");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEGenericBodyForce::FEGenericBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
}

//-----------------------------------------------------------------------------
vec3d FEGenericBodyForce::force(FEMaterialPoint &mp)
{
	return m_force(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEGenericBodyForce::stiffness(FEMaterialPoint& pt)
{
	return mat3ds(0, 0, 0, 0, 0, 0);
}

//=============================================================================

BEGIN_PARAMETER_LIST(FENonConstBodyForceOld, FEGenericBodyForce);
	ADD_PARAMETER(m_force[0], "x");
	ADD_PARAMETER(m_force[1], "y");
	ADD_PARAMETER(m_force[2], "z");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FENonConstBodyForceOld::FENonConstBodyForceOld(FEModel* pfem) : FEGenericBodyForce(pfem)
{
}

//-----------------------------------------------------------------------------
bool FENonConstBodyForceOld::Init()
{
	FEParameterList& PL = GetParameterList();

	FEParam* px = PL.FindFromName("x");
	FEParam* py = PL.FindFromName("y");
	FEParam* pz = PL.FindFromName("z");

	FEParam* paramForce = PL.FindFromName("force");
	FEParamVec3& v = paramForce->value<FEParamVec3>();

	v.setValuator(new FEMathExpressionVec3(m_force[0], m_force[1], m_force[2]));

	paramForce->SetLoadCurve(px->GetLoadCurve());
	px->SetLoadCurve(-1);
	py->SetLoadCurve(-1);
	pz->SetLoadCurve(-1);

	return FEGenericBodyForce::Init();
}
