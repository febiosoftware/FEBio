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



#include "stdafx.h"
#include "FETransform.h"

// conversion factor from degrees to radians (= PI / 180)
#define DEG_TO_RAD 0.01745329252

FETransform::FETransform() : m_pos(0,0,0), m_rot(0, 0, 1)
{
	m_scl[0] = m_scl[1] = m_scl[2] = 1.0;
}

void FETransform::SetTranslation(const vec3d& t)
{
	m_pos = t;
}

void FETransform::SetScale(double sx, double sy, double sz)
{
	m_scl[0] = sx;
	m_scl[1] = sy;
	m_scl[2] = sz;
}

void FETransform::SetRotation(const quatd& q)
{
	m_rot = q;
}

void FETransform::SetRotation(const vec3d& r)
{
	m_rot = quatd(r*DEG_TO_RAD);
}

void FETransform::SetRotation(double X, double Y, double Z)
{
	X *= DEG_TO_RAD;
	Y *= DEG_TO_RAD;
	Z *= DEG_TO_RAD;
	m_rot.SetEuler(X, Y, Z);
}

vec3d FETransform::Transform(const vec3d& r) const
{
	vec3d p(m_scl[0] * r.x, m_scl[1] * r.y, m_scl[2] * r.z);
	m_rot.RotateVector(p);
	p += m_pos;
	return p;
}
