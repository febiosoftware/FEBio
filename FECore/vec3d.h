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



#pragma once
#include <math.h>
#include "vec2d.h"

//! 3D vector class with double precision components
class vec3d
{
public:
	//! Default constructor - initializes vector to (0, 0, 0)
	vec3d() : x(0), y(0), z(0) {}
	//! Constructor with single value - initializes all components to the same value
	explicit vec3d(double a) : x(a), y(a), z(a) {}
	//! Constructor with three components
	vec3d(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
	//! Constructor from 2D vector - sets z component to 0
	vec3d(const vec2d& v) { x = v.r[0]; y = v.r[1]; z = 0.0; }

	//! Vector addition operator
	vec3d operator + (const vec3d& r) const { return vec3d(x+r.x, y+r.y, z+r.z); }
	//! Vector subtraction operator
	vec3d operator - (const vec3d& r) const { return vec3d(x-r.x, y-r.y, z-r.z); }

	//! Scalar multiplication operator
	vec3d operator * (double a) const { return vec3d(x*a, y*a, z*a); }
	//! Scalar division operator
	vec3d operator / (double a) const { return vec3d(x/a, y/a, z/a); }

	//! Vector addition assignment operator
	vec3d& operator += (const vec3d& r) { x += r.x; y += r.y; z += r.z; return (*this); }
	//! Vector subtraction assignment operator
	vec3d& operator -= (const vec3d& r) { x -= r.x; y -= r.y; z -= r.z; return (*this); }

	//! Scalar multiplication assignment operator
	vec3d& operator *= (double a) { x*=a; y*=a; z*=a; return (*this); }
	//! Scalar division assignment operator
	vec3d& operator /= (double a) { x/=a; y/=a; z/=a; return (*this); }

	//! Unary negation operator
	vec3d operator - () const { return vec3d(-x, -y, -z); }

    //! Component access operator (non-const)
    double& operator() (int i)
    {
        switch(i)
        {
            case 0: {return x; break;}
            case 1: {return y; break;}
            case 2: {return z; break;}
            default: {return x; break;}
        }
    }
    
    //! Component access operator (const)
    double operator() (int i) const
    {
        switch(i)
        {
            case 0: {return x; break;}
            case 1: {return y; break;}
            case 2: {return z; break;}
            default: {return x; break;}
        }
    }
    
	//! Dot product operator
	double operator * (const vec3d& r) const { return (x*r.x + y*r.y + z*r.z); }

	//! Cross product operator
	vec3d operator ^ (const vec3d& r) const { return vec3d(y*r.z-z*r.y,z*r.x-x*r.z,x*r.y-y*r.x); }

	//! Normalize the vector in place and return original length
	double unit()
	{
		double d = sqrt(x*x+y*y+z*z);
		if (d != 0) { x/=d; y/=d; z/=d; }
		return d;
	}

	//! Return a normalized copy of this vector
	vec3d normalized() const { 
		double d = sqrt(x*x + y*y + z*z); 
		d = (d == 0.0? d = 1.0 : d = 1.0/d); 
		return vec3d(x*d, y*d, z*d);
	}

	//! Return the length (magnitude) of the vector
	double norm() const { return sqrt(x*x+y*y+z*z); }

	//! Return the squared length of the vector
	double norm2() const { return (x*x + y*y + z*z); }

public:
	//! Normalize the vector in place (FEBio Studio compatibility)
	vec3d Normalize() { unit(); return *this; }
	//! Return a normalized copy of this vector (FEBio Studio compatibility)
	vec3d Normalized() const { vec3d v(x, y, z); v.unit(); return v; }
	//! Return the length of the vector (FEBio Studio compatibility)
	double Length() const { return norm(); }
	//! Return the squared length of the vector (FEBio Studio compatibility)
	double SqrLength() const { return norm2(); }
	//! Equality comparison operator
	bool operator == (const vec3d& a) const { return ((a.x == x) && (a.y == y) && (a.z == z)); }

 public:
	//! X component of the vector
	double x, y, z;
};


//-----------------------------------------------------------------------------
//! 3D vector class with single precision (float) components
class vec3f
{
public:
	//! Default constructor - initializes vector to (0, 0, 0)
	vec3f() { x = y = z = 0; }
	//! Constructor with three float components
	vec3f(float rx, float ry, float rz) { x = rx; y = ry; z = rz; }

	//! Vector addition operator
	vec3f operator + (const vec3f& v) const { return vec3f(x + v.x, y + v.y, z + v.z); }
	//! Vector subtraction operator
	vec3f operator - (const vec3f& v) const { return vec3f(x - v.x, y - v.y, z - v.z); }
	//! Cross product operator
	vec3f operator ^ (const vec3f& v) const
	{
		return vec3f(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
	}

	//! Dot product operator
	float operator * (const vec3f& v) const { return (x * v.x + y * v.y + z * v.z); }

	//! Scalar multiplication operator
	vec3f operator * (const float g) const { return vec3f(x * g, y * g, z * g); }
	//! Scalar division operator
	vec3f operator / (const float g) const { return vec3f(x / g, y / g, z / g); }

	//! Vector addition assignment operator
	const vec3f& operator += (const vec3f& v) { x += v.x; y += v.y; z += v.z; return (*this); }
	//! Vector subtraction assignment operator
	const vec3f& operator -= (const vec3f& v) { x -= v.x; y -= v.y; z -= v.z; return (*this); }
	//! Float division assignment operator
	const vec3f& operator /= (const float& f) { x /= f; y /= f; z /= f; return (*this); }
	//! Integer division assignment operator
	const vec3f& operator /= (const int& n) { x /= n; y /= n; z /= n; return (*this); }
	//! Float multiplication assignment operator
	const vec3f& operator *= (const float& f) { x *= f; y *= f; z *= f; return (*this); }

	//! Unary negation operator
	vec3f operator - () { return vec3f(-x, -y, -z); }

	//! Return the length (magnitude) of the vector
	float Length() const { return (float)sqrt(x * x + y * y + z * z); }

	//! Return the squared length of the vector
	float SqrLength() const { return (float)(x * x + y * y + z * z); }

	//! Normalize the vector in place
	vec3f& Normalize()
	{
		float L = Length();
		if (L != 0) { x /= L; y /= L; z /= L; }

		return (*this);
	}

public:
	//! X component of the vector
	float x, y, z;
};

//! Convert vec3f to vec3d
inline vec3d to_vec3d(const vec3f& r) { return vec3d((double)r.x, (double)r.y, (double)r.z); }
//! Convert vec3d to vec3f
inline vec3f to_vec3f(const vec3d& r) { return vec3f((float)r.x, (float)r.y, (float)r.z); }


//-----------------------------------------------------------------------------
//! 3D vector class with integer components
class vec3i
{
public:
	//! Default constructor - initializes vector to (0, 0, 0)
	vec3i() { x = y = z = 0; }
	//! Constructor with three integer components
	vec3i(int X, int Y, int Z) { x = X; y = Y; z = Z;}

public:
	//! X component of the vector
	int		x, y, z;
};