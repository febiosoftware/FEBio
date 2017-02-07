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
