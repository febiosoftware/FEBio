#include "stdafx.h"
#include "FERigidSurface.h"

///////////////////////////////////////////////////////////////////////////////
// FEPlane
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEPlane, FERigidSurface)
	ADD_PARAMETERV(a, FE_PARAM_DOUBLE, 4, "plane");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPlane::FEPlane(FEModel* pfem) : FERigidSurface(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initializes data for FEPlane

void FEPlane::Init()
{
}

vec3d FEPlane::Normal(const vec3d& r)
{
	vec3d n(a[0], a[1], a[2]);
	n.unit();
	return n;
}

vec3d FEPlane::Project(const vec3d& r)
{
	double d = a[3];

	double l = a[0]*r.x + a[1]*r.y + a[2]*r.z - d;
	return vec3d(r.x-l*a[0], r.y-l*a[1], r.z-l*a[2]);
}

///////////////////////////////////////////////////////////////////////////////
// FERigidSphere
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidSphere, FERigidSurface)
	ADD_PARAMETER(m_R, FE_PARAM_DOUBLE, "radius");
	ADD_PARAMETER(m_rc, FE_PARAM_VEC3D, "center");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor

FERigidSphere::FERigidSphere(FEModel *pfem) : FERigidSurface(pfem)
{
	m_rc = vec3d(0,0,0);
	m_uc = vec3d(0,0,0);
	m_R = 1.0;
}

//-----------------------------------------------------------------------------
//! initialize data for rigid sphere

void FERigidSphere::Init()
{
}

//-----------------------------------------------------------------------------
//! returns the center of the sphere

vec3d FERigidSphere::Center()
{
	vec3d rc = m_rc + m_uc;
	return rc;
}

//-----------------------------------------------------------------------------
//! project node on sphere

vec3d FERigidSphere::Project(const vec3d& r)
{
	vec3d rc = Center();
	vec3d q = r - rc;
	q.unit();
	return rc + q*m_R;
}

//-----------------------------------------------------------------------------
//! return the local normal. This function assumes that r is on the surface

vec3d FERigidSphere::Normal(const vec3d& r)
{
	vec3d q = r - Center();
	q.unit();
	return q;
}
