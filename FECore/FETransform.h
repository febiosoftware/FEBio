#pragma once
#include "vec3d.h"
#include "quatd.h"

//-----------------------------------------------------------------------------
// Class that defines an affine transformation (scale, rotate, translate).
// This currently applies the transformation as follows: first scale, then translate, then rotate.
// 
class FETransform
{
public:
	FETransform();

	// set the scale factors
	void SetScale(double sx, double sy, double sz);

	// set the translation
	void SetTranslation(const vec3d& t);

	// set the rotation
	void SetRotation(const quatd& q);

	// set the rotation using axis (a,b) and rotation angle (degrees)
	void SetRotation(const vec3d& a, const vec3d& b, double angle);

	// apply transformation
	vec3d Transform(const vec3d& r) const;

private:
	double	m_scl[3];	// scale factors
	vec3d	m_pos;		// translation
	quatd	m_rot;		// rotation
	vec3d	m_rc;		// rotation center
};

inline vec3d FETransform::Transform(const vec3d& r) const
{
	vec3d p(m_scl[0]*r.x, m_scl[1]*r.y, m_scl[2]*r.z);

	p += m_pos;
	p -= m_rc;
	m_rot.RotateVector(p);
	p += m_rc;

	return p;
}
