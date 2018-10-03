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
	FEParamVec3 m_force;

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
// Class that implements old body force
class FENonConstBodyForceOld : public FEGenericBodyForce
{
public:
	FENonConstBodyForceOld(FEModel* fem);

	bool Init() override;

private:
	std::string	m_force[3];

	DECLARE_FECORE_CLASS();
};
