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
//! Initializes data for FEPlane

void FEPlane::Init()
{
	// set the plane displacement curve
	if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);
}

bool FEPlane::SetParameterAttribute(FEParam& p, const char* szatt, const char* szval)
{
	if (strcmp(szatt, "lc") == 0)
	{
		m_nplc = atoi(szval) - 1;
		return true;
	}
	else return false;
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
	vec3d rc = m_rc;
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
