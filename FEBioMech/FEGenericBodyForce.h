#pragma once
#include "FEBodyForce.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! This class defines a non-homogeneous body force, i.e. the force can depend
//! on the reference position
class FEGenericBodyForce : public FEBodyForce
{
public:
	//! constructor
	FEGenericBodyForce(FEModel* pfem);

	//! evaluate the body force
	vec3d force(FEMaterialPoint& pt) override;

	//! stiffness 
	mat3ds stiffness(FEMaterialPoint& pt) override;

public:
	FEModelParam m_val[3];

	DECLARE_PARAMETER_LIST();
};
