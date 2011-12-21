#pragma once

#include "FECore/FESurface.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// TODO: introduce parameter lists so that we can remove references to load curves.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//! This class is the base class for rigid surfaces

//! Rigid surfaces are used in the rigid wall contact interface, where the
//! master surface is defined by an implicit surface

class FERigidSurface
{
public: // interface
	FERigidSurface(FEModel* pfem) : m_pfem(pfem) {}

	//! intialize surface
	virtual void Init() = 0;

	//! returns the normal at point r, where r is assumed on the surface
	virtual vec3d Normal(const vec3d& r) = 0;

	//! projects the point on the surface
	virtual vec3d Project(const vec3d& r) = 0;

protected:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! This class implements a rigid plane

//! The FEPlane is used to describe the (moving) rigid wall in a FERigidWallInterface
class FEPlane : public FERigidSurface
{
public:
	//! constructor
	FEPlane(FEModel* pfem) : FERigidSurface(pfem)
	{
		m_nplc = -1;
		m_pplc = 0;
	}

	//! initialization
	void Init();

	//! return plane normal
	vec3d Normal(const vec3d& r)
	{
		vec3d n(a[0], a[1], a[2]);
		n.unit();
		return n;
	}

	//! project node onto plane
	vec3d Project(const vec3d& r)
	{
		double d = a[3];
		if (m_pplc) d += m_pplc->Value();

		double l = a[0]*r.x + a[1]*r.y + a[2]*r.z - d;
		return vec3d(r.x-l*a[0], r.y-l*a[1], r.z-l*a[2]);
	}

	//! get the initial plane equation
	double* GetEquation() { return a; }

protected:
	double	a[4];	//!< plane equation

public:
	int				m_nplc;		//!< plane loadcurve number
	FELoadCurve*	m_pplc;		//!< plane load curve
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

	int				m_nplc[3];	//!< displacement load curve numbers
	FELoadCurve*	m_pplc[3];	//!< displacement load curves
};
