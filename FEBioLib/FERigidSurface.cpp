#include "stdafx.h"
#include "FERigidSurface.h"

///////////////////////////////////////////////////////////////////////////////
// FEPlane
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! Initializes data for FEPlane

void FEPlane::Init()
{
	// set the plane displacement curve
	if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);
}

///////////////////////////////////////////////////////////////////////////////
// FERigidSphere
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor

FERigidSphere::FERigidSphere(FEModel *pfem) : FERigidSurface(pfem)
{
	m_nplc[0] = m_nplc[1] = m_nplc[2] = -1;
	m_pplc[0] = m_pplc[1] = m_pplc[2] = 0;
}

//-----------------------------------------------------------------------------
//! initialize data for rigid sphere

void FERigidSphere::Init()
{
	if (m_nplc[0] >= 0) m_pplc[0] = m_pfem->GetLoadCurve(m_nplc[0]);
	if (m_nplc[1] >= 0) m_pplc[1] = m_pfem->GetLoadCurve(m_nplc[1]);
	if (m_nplc[2] >= 0) m_pplc[2] = m_pfem->GetLoadCurve(m_nplc[2]);
}

//-----------------------------------------------------------------------------
//! returns the center of the sphere

vec3d FERigidSphere::Center()
{
	vec3d rc = m_rc;
	if (m_pplc[0]) rc.x += m_pplc[0]->Value();
	if (m_pplc[1]) rc.y += m_pplc[1]->Value();
	if (m_pplc[2]) rc.z += m_pplc[2]->Value();

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
