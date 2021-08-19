#include "pch.h"
#include "MassDamping.h"
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FEMassDamping, FEBodyForce)
	ADD_PARAMETER(m_C, "C");
END_FECORE_CLASS();

FEMassDamping::FEMassDamping(FEModel* fem) : FEBodyForce(fem)
{
	m_C = 0.0;
}

//! calculate the body force at a material point
vec3d FEMassDamping::force(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = dynamic_cast<FEElasticMaterialPoint&>(mp);
	return ep.m_v*m_C;
}

//! calculate constribution to stiffness matrix
mat3ds FEMassDamping::stiffness(FEMaterialPoint& pt)
{
	FETimeInfo& ti = GetFEModel()->GetTime();
	double dt = ti.timeIncrement;
	return mat3dd(-m_C / dt);
}
