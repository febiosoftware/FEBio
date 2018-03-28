#pragma once
#include "FEBodyForce.h"

//-----------------------------------------------------------------------------
//! This class defines a deformation-independent constant force (e.g. gravity)
class FEConstBodyForce : public FEBodyForce
{
public:
	FEConstBodyForce(FEModel* pfem) : FEBodyForce(pfem) { m_f = vec3d(0,0,0); }
	vec3d force(FEMaterialPoint& pt) override { return m_f; }
	mat3ds stiffness(FEMaterialPoint& pt) override { return mat3ds(0,0,0,0,0,0); }

protected:
	vec3d	m_f;

	DECLARE_PARAMETER_LIST();
};
