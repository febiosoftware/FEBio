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
class FETransform
{
public:
	FETransform();

	// set the scale factors
	void SetScale(double sx, double sy, double sz);

	// set the translation
	void SetTranslation(const vec3d& t);

	// set the rotation quaternion
	void SetRotation(const quatd& q);

	// set the rotation vector (uses degrees)
	void SetRotation(const vec3d& r);

	// set rotation via Euler angles Tait-Bryan (Z,Y,X) convention (in degrees)
	void SetRotation(double X, double Y, double Z);

	// apply transformation
	vec3d Transform(const vec3d& r) const;

private:
	double	m_scl[3];	// scale factors
	vec3d	m_pos;		// translation (global space)
	quatd	m_rot;		// rotation
};
