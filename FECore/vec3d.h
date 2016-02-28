#ifndef _VEC3D_H_10222006_
#define _VEC3D_H_10222006_

#include <math.h>
#include "vec2d.h"

class vec3d
{
public:
	// constructors
	vec3d() : x(0), y(0), z(0) {}
	vec3d(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
	vec3d(const vec2d& v) { x = v.r[0]; y = v.r[1]; z = 0.0; }

	// operators
	vec3d operator + (const vec3d& r) const { return vec3d(x+r.x, y+r.y, z+r.z); }
	vec3d operator - (const vec3d& r) const { return vec3d(x-r.x, y-r.y, z-r.z); }

	vec3d operator * (double a) const { return vec3d(x*a, y*a, z*a); }
	vec3d operator / (double a) const { return vec3d(x/a, y/a, z/a); }

	vec3d& operator += (const vec3d& r) { x += r.x; y += r.y; z += r.z; return (*this); }
	vec3d& operator -= (const vec3d& r) { x -= r.x; y -= r.y; z -= r.z; return (*this); }

	vec3d& operator *= (double a) { x*=a; y*=a; z*=a; return (*this); }
	vec3d& operator /= (double a) { x/=a; y/=a; z/=a; return (*this); }

	vec3d operator - () const { return vec3d(-x, -y, -z); }

	// dot product
	double operator * (const vec3d& r) const { return (x*r.x + y*r.y + z*r.z); }

	// cross product
	vec3d operator ^ (const vec3d& r) const { return vec3d(y*r.z-z*r.y,z*r.x-x*r.z,x*r.y-y*r.x); }

	// normalize the vector
	double unit()
	{
		double d = sqrt(x*x+y*y+z*z);
		if (d != 0) { x/=d; y/=d; z/=d; }
		return d;
	}

	double norm() const { return sqrt(x*x+y*y+z*z); }

public:
	double x, y, z;
};

#endif // _VEC3D_H_10222006_
