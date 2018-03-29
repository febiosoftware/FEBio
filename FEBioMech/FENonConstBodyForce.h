#pragma once
#include "FEBodyForce.h"
#include <FECore/FEMathValue.h>

//-----------------------------------------------------------------------------
//! This class defines a non-homogeneous force, i.e. the force depends
//! on the spatial position
class FENonConstBodyForce : public FEBodyForce
{
public:
	FENonConstBodyForce(FEModel* pfem);
	vec3d force(FEMaterialPoint& pt) override;
	mat3ds stiffness(FEMaterialPoint& pt) override;

public:
	FEMathDouble m_val[3];

	DECLARE_PARAMETER_LIST();
};
