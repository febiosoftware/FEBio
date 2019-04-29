/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

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



#pragma once
#include <FECore/FECoreBase.h>
#include "febiomech_api.h"
//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This class is the base class for rigid surfaces

//! Rigid surfaces are used in the rigid sliding contact interface, where the
//! master surface is defined by an implicit surface

class FEBIOMECH_API FERigidSurface : public FECoreBase
{
	FECORE_SUPER_CLASS

public: // interface
	FERigidSurface(FEModel* fem) : FECoreBase(fem) {}

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
	bool Init() override;

	//! return plane normal
	vec3d Normal(const vec3d& r) override;

	//! project node onto plane
	vec3d Project(const vec3d& r) override;

	//! get the initial plane equation
	double* GetEquation() { return a; }

public:
	double	a[4];	//!< plane equation

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Rigid Sphere class

class FERigidSphere : public FERigidSurface
{
public:
	//! constructor
	FERigidSphere(FEModel* pfem);

	//! initialization
	bool Init() override;

	//! return the normal
	vec3d Normal(const vec3d& r) override;

	//! project on surface
	vec3d Project(const vec3d& r) override;

protected:
	vec3d Center();

public:
	vec3d	m_rc;		//!< center of sphere
	vec3d	m_uc;		//!< displacement of center
	double	m_R;		//!< radius

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Rigid cylinder class
class FERigidCylinder : public FERigidSurface
{
public:
	//! constructor
	FERigidCylinder(FEModel* fem);

	//! Initialization
	bool Init() override;

	//! return closest point projection
	vec3d Project(const vec3d& r) override;

	//! return the normal (assumes point lies on cylinder)
	vec3d Normal(const vec3d& r) override;

private:
	double	m_R;	//!< radius
	vec3d	m_rc;	//!< "center" of cylinder
	vec3d	m_n;	//!< axis of cylinder

	vec3d	m_uc;	//!< center displacement

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Rigid ellipsoid class
class FERigidEllipsoid : public FERigidSurface
{
public:
	//! constructor
	FERigidEllipsoid(FEModel* fem);

	//! Initialization
	bool Init() override;

	//! return closest point projection
	vec3d Project(const vec3d& r) override;

	//! return the normal (assumes point lies on cylinder)
	vec3d Normal(const vec3d& r) override;

private:
	double	m_R[3];	//!< principle axes radii
	vec3d	m_rc;	//!< "center" of ellipsoid
	mat3d	m_Q, m_Qt;	//!< principle axes of cylinder (and transpose)

	vec3d	m_a, m_d;	//!< vectors used to generate principle axes

	vec3d	m_uc;	//!< center displacement

	DECLARE_FECORE_CLASS();
};
