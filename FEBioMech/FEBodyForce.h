#pragma once
#include "FECore/FEMaterialPoint.h"
#include "FECore/FEBodyLoad.h"

//-----------------------------------------------------------------------------
//! This class is the base class for body forces
//! Derived classes need to implement the force and stiffness functions.
//
class FEBodyForce : public FEBodyLoad
{
public:
	//! constructor
	FEBodyForce(FEModel* pfem);

	//! update
	virtual void Update(){}

public:
	//! calculate the body force at a material point
	virtual vec3d force(FEMaterialPoint& pt) = 0;

	//! calculate constribution to stiffness matrix
	virtual mat3ds stiffness(FEMaterialPoint& pt) = 0;
};
