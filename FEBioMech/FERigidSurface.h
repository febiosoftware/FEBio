#pragma once
#include <FECore/FECoreBase.h>

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This class is the base class for rigid surfaces

//! Rigid surfaces are used in the rigid wall contact interface, where the
//! master surface is defined by an implicit surface

//! \todo Introduce parameter lists so that we can remove references to load curves.

class FERigidSurface : public FECoreBase
{
public: // interface
	FERigidSurface(FEModel* pfem) : FECoreBase(FERIGIDOBJECT_ID) {}

	//! intialize surface
	virtual void Init() = 0;

	//! returns the normal at point r, where r is assumed on the surface
	virtual vec3d Normal(const vec3d& r) = 0;

	//! projects the point on the surface
	virtual vec3d Project(const vec3d& r) = 0;
};

//-----------------------------------------------------------------------------
//! This class implements a rigid plane

//! The FEPlane is used to describe the (moving) rigid wall in a FERigidWallInterface
class FEPlane : public FERigidSurface
{
public:
	//! constructor
	FEPlane(FEModel* pfem);

	//! initialization
	void Init();

	//! return plane normal
	vec3d Normal(const vec3d& r);

	//! project node onto plane
	vec3d Project(const vec3d& r);

	//! get the initial plane equation
	double* GetEquation() { return a; }

public:
	double	a[4];	//!< plane equation

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Rigid Sphere class

class FERigidSphere : public FERigidSurface
{
public:
	//! constructor
	FERigidSphere(FEModel* pfem);

	//! initialization
	void Init();

	//! return the normal
	vec3d Normal(const vec3d& r);

	//! project on surface
	vec3d Project(const vec3d& r);

protected:
	vec3d Center();

public:
	vec3d	m_rc;		//!< center of sphere
	double	m_R;		//!< radius

	DECLARE_PARAMETER_LIST();
};
