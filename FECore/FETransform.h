/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "vec3d.h"
#include "quatd.h"

//-----------------------------------------------------------------------------
// Class that defines an affine transformation (scale, rotate, translate).
// This currently applies the transformation as follows: 
// 1. scale : the scale is applied in the local coordinate system
// 2. rotate: rotation from local to global coordinates
// 3. translate: translate to a global position
// 
class Transform
{
public:
	FECORE_API Transform();

	// Reset the transform
	FECORE_API void Reset();

	// set the scale factors
	FECORE_API void SetScale(double sx, double sy, double sz);

	// set the scale of the object
	void SetScale(const vec3d& s) { m_scl = s; }

	//! get scale of the object
	const vec3d& GetScale() const { return m_scl; }

	// set the position (or translation)
	FECORE_API void SetPosition(const vec3d& t);

	// get position of object
	const vec3d& GetPosition() const { return m_pos; }

	// set the rotation quaternion
	FECORE_API void SetRotation(const quatd& q);

	// set the rotation vector (uses degrees)
	FECORE_API void SetRotation(const vec3d& r);

	// set rotation via Euler angles Tait-Bryan (Z,Y,X) convention (in degrees)
	FECORE_API void SetRotation(double X, double Y, double Z);

	//! get orientation
	const quatd& GetRotation() const { return m_rot; }

	// get inverse of rotation
	quatd GetRotationInverse() const { return m_roti; }

	// apply transformation
	FECORE_API vec3d Apply(const vec3d& r) const;

	// translate the transform
	void Translate(const vec3d& dr);

	// scale an object
	void Scale(double s, vec3d r, vec3d rc);

	// rotate around the center rc
	void Rotate(quatd q, vec3d rc);

	// Rotate angle w around an axis defined by the position vectors a, b.
	void Rotate(const vec3d& a, const vec3d& b, double w);

	// comparison
	bool operator == (const Transform& T) const;

public:
	// convert from local to global coordinates
	vec3d LocalToGlobal(const vec3d& r) const;

	// convert from global to local coordinates
	vec3d GlobalToLocal(const vec3d& r) const;

	//! get a normal-like vector from global to local
	vec3d LocalToGlobalNormal(const vec3d& n) const;

	//! get a normal-like vector from global to local
	vec3d GlobalToLocalNormal(const vec3d& n) const;

private:
	vec3d	m_scl;		// scale factors
	vec3d	m_pos;		// translation (global space)
	quatd	m_rot;		// rotation
	quatd	m_roti;		// inverse rotation
};

inline bool Transform::operator == (const Transform& T) const
{
	return ((m_pos == T.m_pos) && (m_scl == T.m_scl) && (m_rot == T.m_rot));
}

inline void Transform::Translate(const vec3d& dr) { m_pos += dr; }

// convert from local to global coordinates
inline vec3d Transform::LocalToGlobal(const vec3d& r) const
{
	return m_pos + m_rot * vec3d(r.x * m_scl.x, r.y * m_scl.y, r.z * m_scl.z);
}

// convert from global to local coordinates
inline vec3d Transform::GlobalToLocal(const vec3d& r) const
{
	vec3d p = m_roti * (r - m_pos);
	return vec3d(p.x / m_scl.x, p.y / m_scl.y, p.z / m_scl.z);
}

//! get a normal-like vector from global to local
inline vec3d Transform::LocalToGlobalNormal(const vec3d& n) const
{
	// NOTE: scaling is turned off because this is used in the generation of material axes.
	//       If I use scaling the axes may no longer be orthogonal. Maybe I should create another
	//       function for this since this is now inconsistent with the reverse operation.
//		return m_rot*vec3d(n.x / m_scl.x, n.y / m_scl.y, n.z / m_scl.z);
	return m_rot * vec3d(n.x, n.y, n.z);
}

//! get a normal-like vector from global to local
inline vec3d Transform::GlobalToLocalNormal(const vec3d& n) const
{
	vec3d m = m_roti * n;
	m.x /= m_scl.x; m.y /= m_scl.y; m.z /= m_scl.z;
	m.Normalize();
	return m;
}

// scale 
inline void Transform::Scale(double s, vec3d r, vec3d rc)
{
	vec3d r0 = GlobalToLocal(rc);

	double a = s - 1;
	m_roti.RotateVector(r);
	r.Normalize();

	r.x = 1 + a * fabs(r.x);
	r.y = 1 + a * fabs(r.y);
	r.z = 1 + a * fabs(r.z);

	m_scl.x *= r.x;
	m_scl.y *= r.y;
	m_scl.z *= r.z;

	m_pos -= LocalToGlobal(r0) - rc;
}

// rotate around the center rc
inline void Transform::Rotate(quatd q, vec3d rc)
{
	m_rot = q * m_rot;
	m_roti = m_rot.Inverse();

	m_rot.MakeUnit();
	m_roti.MakeUnit();

	m_pos = rc + q * (m_pos - rc);
}

// Rotate angle w around an axis defined by the position vectors a, b.
inline void Transform::Rotate(const vec3d& a, const vec3d& b, double w)
{
	double wr = PI * w / 180.0;
	vec3d n = (b - a); n.Normalize();
	quatd q(wr, n);
	Rotate(q, a);
}
