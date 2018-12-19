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
BEGIN_FECORE_CLASS(FENonConstBodyForceOld, FEBodyForce);
	ADD_PARAMETER(m_f[0], "x");
	ADD_PARAMETER(m_f[1], "y");
	ADD_PARAMETER(m_f[2], "z");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENonConstBodyForceOld::FENonConstBodyForceOld(FEModel* pfem) : FEBodyForce(pfem)
{
}

//-----------------------------------------------------------------------------
vec3d FENonConstBodyForceOld::force(FEMaterialPoint& pt)
{
	vec3d F;
	F.x = m_f[0](pt);
	F.y = m_f[0](pt);
	F.z = m_f[0](pt);
	return F;
}
