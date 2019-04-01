/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "stdafx.h"
#include "FERigidSurface.h"

REGISTER_SUPER_CLASS(FERigidSurface, FERIGIDOBJECT_ID);

///////////////////////////////////////////////////////////////////////////////
// FEPlane
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidPlane, FERigidSurface)
	ADD_PARAMETER(a, 4, "plane");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidPlane::FERigidPlane(FEModel* pfem) : FERigidSurface(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initializes data for FEPlane

bool FERigidPlane::Init()
{
	return FERigidSurface::Init();
}

vec3d FERigidPlane::Normal(const vec3d& r)
{
	vec3d n(a[0], a[1], a[2]);
	n.unit();
	return n;
}

vec3d FERigidPlane::Project(const vec3d& r)
{
	double d = a[3];

	double l = a[0]*r.x + a[1]*r.y + a[2]*r.z - d;
	return vec3d(r.x-l*a[0], r.y-l*a[1], r.z-l*a[2]);
}

///////////////////////////////////////////////////////////////////////////////
// FERigidSphere
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidSphere, FERigidSurface)
	ADD_PARAMETER(m_R, "radius");
	ADD_PARAMETER(m_rc, "center");
	ADD_PARAMETER(m_uc.x, "ux");
	ADD_PARAMETER(m_uc.y, "uy");
	ADD_PARAMETER(m_uc.z, "uz");
END_FECORE_CLASS();

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

bool FERigidSphere::Init()
{
	return FERigidSurface::Init();
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

///////////////////////////////////////////////////////////////////////////////
// FERigidCylinder
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidCylinder, FERigidSurface)
	ADD_PARAMETER(m_R, "radius");
	ADD_PARAMETER(m_rc, "center");
	ADD_PARAMETER(m_n , "axis"  );
	ADD_PARAMETER(m_uc.x, "ux");
	ADD_PARAMETER(m_uc.y, "uy");
	ADD_PARAMETER(m_uc.z, "uz");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor

FERigidCylinder::FERigidCylinder(FEModel *pfem) : FERigidSurface(pfem)
{
	m_rc = vec3d(0,0,0);
	m_uc = vec3d(0,0,0);
	m_n  = vec3d(0,0,1);
	m_R = 1.0;
}

//-----------------------------------------------------------------------------
//! initialize data for rigid sphere

bool FERigidCylinder::Init()
{
	return FERigidSurface::Init();
}

//-----------------------------------------------------------------------------
//! project node on sphere

vec3d FERigidCylinder::Project(const vec3d& r)
{
	vec3d c = m_rc + m_uc;

	double NN = m_n*m_n;
	vec3d p = c + m_n*(((r - c)*m_n)/NN);
	vec3d d = r - p;
	d.unit();
	return p + d*m_R;
}

//-----------------------------------------------------------------------------
//! return the local normal. This function assumes that r is on the surface

vec3d FERigidCylinder::Normal(const vec3d& r)
{
	vec3d c = m_rc + m_uc;

	double NN = m_n*m_n;
	vec3d p = c + m_n*(((r - c)*m_n) / NN);
	vec3d d = r - p;

	d.unit();

	return d;
}

///////////////////////////////////////////////////////////////////////////////
// FERigidEllipsoid
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidEllipsoid, FERigidSurface)
	ADD_PARAMETER(m_R[0], "x_radius");
	ADD_PARAMETER(m_R[1], "y_radius");
	ADD_PARAMETER(m_R[2], "z_radius");
	ADD_PARAMETER(m_rc, "center");
	ADD_PARAMETER(m_a, "a"  );
	ADD_PARAMETER(m_d, "d");
	ADD_PARAMETER(m_uc.x, "ux");
	ADD_PARAMETER(m_uc.y, "uy");
	ADD_PARAMETER(m_uc.z, "uz");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor

FERigidEllipsoid::FERigidEllipsoid(FEModel *pfem) : FERigidSurface(pfem)
{
	m_rc = vec3d(0,0,0);
	m_uc = vec3d(0,0,0);

	m_R[0] = m_R[1] = m_R[2] = 0.0;
	m_a = vec3d(1,0,0);
	m_d = vec3d(0,1,0);

	m_Q.unit();
	m_Qt.unit();
}

//-----------------------------------------------------------------------------
//! initialize data for rigid sphere

bool FERigidEllipsoid::Init()
{
	vec3d e1 = m_a; e1.unit();
	vec3d e3 = m_a ^ m_d; e3.unit();
	vec3d e2 = e3 ^ e1; e2.unit();

	m_Q[0][0] = e1.x; m_Q[0][1] = e2.x; m_Q[0][2] = e3.x;
	m_Q[1][0] = e1.y; m_Q[1][1] = e2.y; m_Q[1][2] = e3.y;
	m_Q[2][0] = e1.z; m_Q[2][1] = e2.z; m_Q[2][2] = e3.z;

	m_Qt = m_Q.transpose();

	return FERigidSurface::Init();
}

//-----------------------------------------------------------------------------
//! project node on sphere

vec3d FERigidEllipsoid::Project(const vec3d& r)
{
	// get center
	vec3d c = m_rc + m_uc;

	// calculate relative position in local coordinates
	vec3d p = m_Qt*(r - c);
	p.x /= m_R[0];
	p.y /= m_R[1];
	p.z /= m_R[2];
	p.unit();

	// scale back by radii
	p.x *= m_R[0];
	p.y *= m_R[1];
	p.z *= m_R[2];

	// rotate back to global coordinates
	p = m_Q*p;

	// don't forget to add translation
	return c + p;
}

//-----------------------------------------------------------------------------
//! return the local normal. This function assumes that r is on the surface

vec3d FERigidEllipsoid::Normal(const vec3d& r)
{
	vec3d c = m_rc + m_uc;

	// calculate relative position in local coordinates
	vec3d p = m_Qt*(r - c);

	// calculate normal
	vec3d n;
	n.x = p.x / (m_R[0] * m_R[0]);
	n.y = p.y / (m_R[1] * m_R[1]);
	n.z = p.z / (m_R[2] * m_R[2]);
	n.unit();

	// rotate back to global coordinates
	n = m_Q*n;

	return n;
}
