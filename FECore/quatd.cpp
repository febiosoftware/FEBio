#include "stdafx.h"
#include "quatd.h"
#include <math.h>

#ifndef PI
#define PI 3.14159265358979
#endif

//-----------------------------------------------------------------------------
//! Spherical linear interpolation between two quaternions.
quatd quatd::slerp(quatd &q1, quatd &q2, const double t) 
{
	quatd q3;
	double dot = quatd::dot(q1, q2);

	/*	dot = cos(theta)
		if (dot < 0), q1 and q2 are more than 90 degrees apart,
		so we can invert one to reduce spinning	*/
	if (dot < 0)
	{
		dot = -dot;
		q3 = -q2;
	} else q3 = q2;
		
	if (dot < 0.95f)
	{
		double angle = acos(dot);
		return (q1*sin(angle*(1-t)) + q3*sin(angle*t))/sin(angle);
	} else // if the angle is small, use linear interpolation								
		return quatd::lerp(q1,q3,t);
}


//-----------------------------------------------------------------------------
void quatd::SetEuler(double X, double Y, double Z)
{
	// calculate cos and sin of angles
	double cz = cos(Z*0.5);
	double sz = sin(Z*0.5);
	double cx = cos(X*0.5);
	double sx = sin(X*0.5);
	double cy = cos(Y*0.5);
	double sy = sin(Y*0.5);

	// define quaternion
	w = cz * cx * cy + sz * sx * sy;
	x = cz * sx * cy - sz * cx * sy;
	y = cz * cx * sy + sz * sx * cy;
	z = sz * cx * cy - cz * sx * sy;
}

//-----------------------------------------------------------------------------
void quatd::GetEuler(double& X, double& Y, double& Z) const
{
	// roll (x-axis rotation)
	double t0 = +2.0 * (w * x + y * z);
	double t1 = +1.0 - 2.0 * (x*x + y*y);
	X = atan2(t0, t1);

	// pitch (y-axis rotation)
	double t2 = +2.0 * (w*y - z*x);
	t2 = t2 > 1.0 ? 1.0 : t2;
	t2 = t2 < -1.0 ? -1.0 : t2;
	Y = asin(t2);

	// yaw (z-axis rotation)
	double t3 = +2.0 * (w * z + x * y);
	double t4 = +1.0 - 2.0 * (y*y + z*z);
	Z = atan2(t3, t4);
}

//-----------------------------------------------------------------------------
// convert euler angles to a rotation matrix
// l[0] = psi   (x-rotation)
// l[1] = theta (y-rotation)
// l[2] = phi   (z-rotation)
mat3d euler2rot(double l[3])
{
	double c0 = cos(l[0]), s0 = sin(l[0]);
	double c1 = cos(l[1]), s1 = sin(l[1]);
	double c2 = cos(l[2]), s2 = sin(l[2]);
	mat3d Rx(1.0, 0.0, 0.0, 0.0, c0, -s0, 0.0, s0, c0);
	mat3d Ry(c1, 0.0, s1, 0.0, 1.0, 0.0, -s1, 0.0, c1);
	mat3d Rz(c2, -s2, 0.0, s2, c2, 0.0, 0.0, 0.0, 1.0);
	return Rz*Ry*Rx;
}

//-----------------------------------------------------------------------------
// convert a rotation matrix to euler angles
// l[0] = psi   (x-rotation)
// l[1] = theta (y-rotation)
// l[2] = phi   (z-rotation)
void rot2euler(mat3d& m, double l[3])
{
	const double e = 1e-12;
	if (fabs(m(2,0) - 1.0) < e)
	{
		if (m(2,0) < 0.0)
		{
			l[2] = 0.0;
			l[1] = PI/2.0;
			l[0] = atan2(m(0,1),m(0,2));
		}
		else
		{
			l[2] = 0.0;
			l[1] = -PI/2.0;
			l[0] = atan2(-m(0,1),-m(0,2));
		}
	}
	else
	{
		l[1] = -asin(m(2,0));
		double c1 = cos(l[1]);
		l[0] = atan2(m(2,1)/c1, m(2,2)/c1);
		l[2] = atan2(m(1,0)/c1, m(0,0)/c1);
	}
}

//-----------------------------------------------------------------------------
// convert a quaternion to Euler angles
void quat2euler(quatd& q, double l[3])
{
	mat3d Q = q.RotationMatrix();
	rot2euler(Q, l);
}
