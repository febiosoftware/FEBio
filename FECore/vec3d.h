/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include <math.h>
#include "vec2d.h"

class vec3d
{
public:
	// constructors
	vec3d() : x(0), y(0), z(0) {}
	explicit vec3d(double a) : x(a), y(a), z(a) {}
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

	// return a normalized version of this vector
	vec3d normalized() const { 
		double d = sqrt(x*x + y*y + z*z); 
		d = (d == 0.0? d = 1.0 : d = 1.0/d); 
		return vec3d(x*d, y*d, z*d);
	}

	// length of vector
	double norm() const { return sqrt(x*x+y*y+z*z); }

	// length square of vector
	double norm2() const { return (x*x + y*y + z*z); }

public:
	double x, y, z;
};
