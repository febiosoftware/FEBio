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
//! This class defines a deformation-independent constant force (e.g. gravity)
// NOTE: This class is obsolete!
class FEConstBodyForceOld : public FEBodyForce
{
public:
	FEConstBodyForceOld(FEModel* pfem) : FEBodyForce(pfem) { m_f = vec3d(0, 0, 0); }
	vec3d force(FEMaterialPoint& pt) override { return m_f; }
	mat3ds stiffness(FEMaterialPoint& pt) override { return mat3ds(0, 0, 0, 0, 0, 0); }

protected:
	vec3d	m_f;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Class that implements old body force
// NOTE: This class is obsolete!
class FENonConstBodyForceOld : public FEBodyForce
{
public:
	FENonConstBodyForceOld(FEModel* fem);
	vec3d force(FEMaterialPoint& pt) override;
	mat3ds stiffness(FEMaterialPoint& pt) override { return mat3ds(0, 0, 0, 0, 0, 0); }

private:
	FEParamDouble	m_f[3];

	DECLARE_FECORE_CLASS();
};
