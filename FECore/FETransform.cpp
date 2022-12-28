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



#include "stdafx.h"
#include "FETransform.h"

Transform::Transform()
{
	Reset();
}

void Transform::Reset()
{
	m_scl = vec3d(1, 1, 1);
	m_pos = vec3d(0, 0, 0);
	m_rot = quatd(0, 0, 0, 1);
	m_roti = quatd(0, 0, 0, 1);
}

void Transform::SetPosition(const vec3d& t)
{
	m_pos = t;
}

void Transform::SetScale(double sx, double sy, double sz)
{
	m_scl.x = sx;
	m_scl.y = sy;
	m_scl.z = sz;
}

void Transform::SetRotation(const quatd& q)
{
	m_rot = q;
	m_roti = m_rot.Inverse();
}

void Transform::SetRotation(const vec3d& r)
{
	m_rot = quatd(r*DEG2RAD);
	m_roti = m_rot.Inverse();
}

void Transform::SetRotation(double X, double Y, double Z)
{
	X *= DEG2RAD;
	Y *= DEG2RAD;
	Z *= DEG2RAD;
	m_rot.SetEuler(X, Y, Z);
	m_roti = m_rot.Inverse();
}

vec3d Transform::Apply(const vec3d& r) const
{
	vec3d p(m_scl.x * r.x, m_scl.y * r.y, m_scl.z * r.z);
	m_rot.RotateVector(p);
	p += m_pos;
	return p;
}
