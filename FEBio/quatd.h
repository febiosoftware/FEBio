#ifndef _QUATD_H_02082007_
#define _QUATD_H_02082007_

#include "vec3d.h"

///////////////////////////////////////////////////////////////////
// CLASS: quatd
// This class implements a quaternion
//

class quatd
{
public:
	// constructors
	quatd () { x = y = z = w = 0.f; }

	quatd( const double angle, vec3d v)
	{
		w = (double) cos(angle * 0.5);

		double sina = (double) sin(angle * 0.5);

		v.unit();
		
		x = v.x*sina;
		y = v.y*sina;
		z = v.z*sina;
	}

	quatd (vec3d& v1, vec3d& v2)
	{
		vec3d n = v1^v2;
		n.unit();

		double d = v1*v2;

		double sina = (double) sqrt((1.0-d)*0.5);
		double cosa = (double) sqrt((1.0+d)*0.5);

		w = cosa;

		x = n.x*sina;
		y = n.y*sina;
		z = n.z*sina;

	}

	quatd(const double qx, const double qy, const double qz, const double qw = 0.0)
	{
		w = qw;
		x = qx;
		y = qy;
		z = qz;
	}

	bool operator != (const quatd& q) { return ((x!=q.x) || (y!=q.y) || (z!=q.z) || (w!=q.w)); }

	// addition and substraction

	quatd operator + (const quatd& q) const
	{
		return quatd(x + q.x, y + q.y, z + q.z, w + q.w);
	}

	quatd operator - (const quatd& q) const
	{
		return quatd(x - q.x, y - q.y, z - q.z, w - q.w);
	}

	quatd& operator += (const quatd& q)
	{
		x += q.x;
		y += q.y;
		z += q.z;
		w += q.w;

		return *this;
	}

	quatd& operator -= (const quatd& q)
	{
		x -= q.x;
		y -= q.y;
		z -= q.z;
		w -= q.w;

		return *this;
	}


	// multiplication

	quatd operator * (const quatd& q) const
	{
		double qw = w*q.w - x*q.x - y*q.y - z*q.z;
		double qx = w*q.x + x*q.w + y*q.z - z*q.y;
		double qy = w*q.y + y*q.w + z*q.x - x*q.z;
		double qz = w*q.z + z*q.w + x*q.y - y*q.x;

		return quatd(qx, qy, qz, qw);
	}

	quatd& operator *= (const quatd& q)
	{
		double qw = w*q.w - x*q.x - y*q.y - z*q.z;
		double qx = w*q.x + x*q.w + y*q.z - z*q.y;
		double qy = w*q.y + y*q.w + z*q.x - x*q.z;
		double qz = w*q.z + z*q.w + x*q.y - y*q.x;

		x = qx;
		y = qy;
		z = qz;
		w = qw;

		return *this;
	}

	quatd operator*(const double a) const
	{
		return quatd(x*a, y*a, z*a, w*a);
	}

	// division

	quatd operator / (const double a) const
	{
		return quatd(x/a, y/a, z/a, w/a);
	}

	quatd& operator /= (const double a)
	{
		x /= a;
		y /= a;
		z /= a;
		w /= a;

		return *this;
	}

	// Special ops

	quatd Conjugate() const { return quatd(-x, -y, -z, w); }

	double Norm() const { return w*w + x*x + y*y + z*z; } 

	void MakeUnit() 
	{
		double N = (double) sqrt(w*w + x*x + y*y + z*z);

		if (N != 0)
		{
			x /= N;
			y /= N;
			z /= N;
			w /= N;
		}
		else w = 1.f;
	}

	quatd Inverse() const
	{
		double N = w*w + x*x + y*y + z*z;

		return quatd(-x/N, -y/N, -z/N, w/N);
	}

	double DotProduct(const quatd& q) const
	{
		return w*q.w + x*q.x + y*q.y + z*q.z;
	}

	vec3d GetVector() const
	{
		vec3d r(x,y,z);
		r.unit();
		return r;
	}

	double GetAngle() const
	{
		return (double)(acos(w)*2.0);
	}

/*	quatd& MultiplyAngle(double fa)
	{
		double angle = fa*acos(w)*2.0;

		w = cos(angle * 0.5);

		double sina = sin(angle * 0.5);

		x *= sina;
		y *= sina;
		z *= sina;
	}
*/


	// use only when *this is unit vector
	void RotateVector(vec3d& v) const
	{
		if ((w == 0) || ((x==0) && (y==0) && (z==0))) return;

		// v*q^-1
		double qw = v.x*x + v.y*y + v.z*z;
		double qx = v.x*w - v.y*z + v.z*y;
		double qy = v.y*w - v.z*x + v.x*z;
		double qz = v.z*w - v.x*y + v.y*x;

		// q* (v* q^-1)
		v.x = (double) (w*qx + x*qw + y*qz - z*qy);
		v.y = (double) (w*qy + y*qw + z*qx - x*qz);
		v.z = (double) (w*qz + z*qw + x*qy - y*qx);
	}

	// use only when *this is unit vector
	vec3d operator * (const vec3d& r)
	{
		vec3d n = r;

		// v*q^-1
		double qw = n.x*x + n.y*y + n.z*z;
		double qx = n.x*w - n.y*z + n.z*y;
		double qy = n.y*w - n.z*x + n.x*z;
		double qz = n.z*w - n.x*y + n.y*x;

		// q* (v* q^-1)
		n.x = (double) (w*qx + x*qw + y*qz - z*qy);
		n.y = (double) (w*qy + y*qw + z*qx - x*qz);
		n.z = (double) (w*qz + z*qw + x*qy - y*qx);

		return n;
	}

	void RotateVectorP(double* v, double* r) const
	{
		static double fx, fy, fz, fw;
		static double qw, qx, qy, qz;
		
		fx = (double) x;
		fy = (double) y;
		fz = (double) z;
		fw = (double) w;

		qw = v[0]*fx + v[1]*fy + v[2]*fz;
		qx = v[0]*fw - v[1]*fz + v[2]*fy;
		qy = v[1]*fw - v[2]*fx + v[0]*fz;
		qz = v[2]*fw - v[0]*fy + v[1]*fx;

		r[0] = (double) (fw*qx + fx*qw + fy*qz - fz*qy);
		r[1] = (double) (fw*qy + fy*qw + fz*qx - fx*qz);
		r[2] = (double) (fw*qz + fz*qw + fx*qy - fy*qx);
	}

public:
	double w;
	double x, y, z;
};

inline quatd operator * (const double a, const quatd& q)
{
	return q*a;
}


#endif //_QUATD_H_02082007_
