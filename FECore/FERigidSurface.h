#pragma once
#include <FECore/FECoreBase.h>

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This class is the base class for rigid surfaces

//! Rigid surfaces are used in the rigid sliding contact interface, where the
//! master surface is defined by an implicit surface

class FERigidSurface : public FECoreBase
{
public: // interface
	FERigidSurface(FEModel* pfem) : FECoreBase(FERIGIDOBJECT_ID) {}

	//! intialize surface
	virtual bool Init() { return true; }

	//! returns the normal at point r, where r is assumed on the surface
	virtual vec3d Normal(const vec3d& r) = 0;

	//! projects the point on the surface
	virtual vec3d Project(const vec3d& r) = 0;
};

//-----------------------------------------------------------------------------
//! This class implements a rigid plane

//! The FEPlane is used to describe the (moving) rigid wall in a FERigidWallInterface
class FERigidPlane : public FERigidSurface
{
public:
	//! constructor
	FERigidPlane(FEModel* pfem);

	//! initialization
	bool Init();

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
	bool Init();

	//! return the normal
	vec3d Normal(const vec3d& r);

	//! project on surface
	vec3d Project(const vec3d& r);

protected:
	vec3d Center();

public:
	vec3d	m_rc;		//!< center of sphere
	vec3d	m_uc;		//!< displacement of center
	double	m_R;		//!< radius

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Rigid cylinder class
class FERigidCylinder : public FERigidSurface
{
public:
	//! constructor
	FERigidCylinder(FEModel* fem);

	//! Initialization
	bool Init();

	//! return closest point projection
	vec3d Project(const vec3d& r);

	//! return the normal (assumes point lies on cylinder)
	vec3d Normal(const vec3d& r);

private:
	double	m_R;	//!< radius
	vec3d	m_rc;	//!< "center" of cylinder
	vec3d	m_n;	//!< axis of cylinder

	vec3d	m_uc;	//!< center displacement

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Rigid ellipsoid class
class FERigidEllipsoid : public FERigidSurface
{
public:
	//! constructor
	FERigidEllipsoid(FEModel* fem);

	//! Initialization
	bool Init();

	//! return closest point projection
	vec3d Project(const vec3d& r);

	//! return the normal (assumes point lies on cylinder)
	vec3d Normal(const vec3d& r);

private:
	double	m_R[3];	//!< principle axes radii
	vec3d	m_rc;	//!< "center" of ellipsoid
	mat3d	m_Q, m_Qt;	//!< principle axes of cylinder (and transpose)

	vec3d	m_a, m_d;	//!< vectors used to generate principle axes

	vec3d	m_uc;	//!< center displacement

	DECLARE_PARAMETER_LIST();
};
