#include "stdafx.h"
#include "FETransform.h"

FETransform::FETransform() : m_pos(0,0,0), m_rc(0,0,0), m_rot(0, 0, 1)
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

void FETransform::SetRotation(const vec3d& a, const vec3d& b, double angle)
{
	m_rc = a;
	vec3d N = b - a; N.unit();

	const double PI = 0.25*atan(1.0);
	double w = angle*PI/180.0;
	m_rot = quatd(angle, N);
}

void FETransform::SetRotation(const quatd& q)
{
	m_rot = q;
}
