#include "stdafx.h"
#include "FEGenericBodyForce.h"
#include "FEElasticMaterial.h"
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FEGenericBodyForce, FEBodyForce);
	ADD_PARAMETER(m_force, "force");
END_FECORE_CLASS();

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
BEGIN_FECORE_CLASS(FEConstBodyForceOld, FEBodyForce);
	ADD_PARAMETER(m_f.x, "x");
	ADD_PARAMETER(m_f.y, "y");
	ADD_PARAMETER(m_f.z, "z");
END_FECORE_CLASS();


//=============================================================================
BEGIN_FECORE_CLASS(FENonConstBodyForceOld, FEGenericBodyForce);
	ADD_PARAMETER(m_force[0], "x");
	ADD_PARAMETER(m_force[1], "y");
	ADD_PARAMETER(m_force[2], "z");
END_FECORE_CLASS();

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

	FEModel* fem = GetFEModel();
	FEMathValueVec3* val = new FEMathValueVec3(fem);
	val->create(m_force[0], m_force[1], m_force[2]);
	v.setValuator(val);

	FELoadController* plc = fem->GetLoadController(px);
	fem->AttachLoadController(paramForce, plc);
	fem->DetachLoadController(px);
	fem->DetachLoadController(py);
	fem->DetachLoadController(pz);

	return FEGenericBodyForce::Init();
}
