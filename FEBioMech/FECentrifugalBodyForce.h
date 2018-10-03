#pragma once
#include "FEBodyForce.h"

//-----------------------------------------------------------------------------
//! This class defines a centrigufal force

class FECentrifugalBodyForce : public FEBodyForce
{
public:
	FECentrifugalBodyForce(FEModel* pfem);

	vec3d force(FEMaterialPoint& mp) override;

	mat3ds stiffness(FEMaterialPoint& mp) override;

public:
	vec3d	n;	// rotation axis
	vec3d	c;	// point on axis of rotation (e.g., center of rotation)
	double	w;	// angular speed

	DECLARE_FECORE_CLASS();
};
