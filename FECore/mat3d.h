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
#include <assert.h>
#include "vec3d.h"
#include "mat2d.h"

//-----------------------------------------------------------------------------
// useful constants for trig
#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef RAD2DEG
#define RAD2DEG (180.0/PI)
#endif

#ifndef DEG2RAD
#define DEG2RAD (PI/180.0)
#endif

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

//-----------------------------------------------------------------------------
// The following classes are defined in this file
class mat3d;	// general 3D matrix of doubles
class mat3ds;	// symmetric 3D matrix of doubles
class mat3da;	// anti-symmetric 3D matrix of doubles
class mat3dd;	// diagonal matrix of doubles

//-----------------------------------------------------------------------------
//! This class describes a diagonal matrix of doubles in 3D

class mat3dd
{
public:
	// default constructor
	mat3dd(){}

	// constructors
	explicit mat3dd(double a);
	mat3dd(double a0, double a1, double a2);

	// assignment operators
	mat3dd& operator = (const mat3dd& m);
	mat3dd& operator = (double a);

	// access operators
	double operator () (int i, int j) const;
	double& diag(int i);
	const double& diag(int i) const;

	// arithmetic operators
	mat3dd operator + (const mat3dd& m) const;
	mat3dd operator - (const mat3dd& m) const;
	mat3dd operator * (const mat3dd& m) const;
	mat3dd operator * (double a) const;
	mat3dd operator / (double a) const;

	mat3dd operator - () const;

	// arithmetic operators for mat3ds
	mat3ds operator + (const mat3ds& m) const;
	mat3ds operator - (const mat3ds& m) const;
	mat3ds operator * (const mat3ds& m) const;

	// arithmetic operators for mat3d
	mat3d operator + (const mat3d& m) const;
	mat3d operator - (const mat3d& m) const;
	mat3d operator * (const mat3d& m) const;

	// arithmetic operators for mat3da const;
	mat3d operator + (const mat3da& m) const;
	mat3d operator - (const mat3da& m) const;
	mat3d operator * (const mat3da& m) const;

	// arithmetic assignment operators
	mat3dd& operator += (const mat3dd& m);
	mat3dd& operator -= (const mat3dd& m);
	mat3dd& operator *= (const mat3dd& m);
	mat3dd& operator *= (double a);
	mat3dd& operator /= (double a);

	// matrix-vector multiplication
	vec3d operator * (const vec3d& r) const;

	// trace
	double tr() const;

	// determinant
	double det() const;

	double xx() const { return d[0]; }
	double yy() const { return d[1]; }
	double zz() const { return d[2]; }

	// TODO: Make these constexpr
	double xy() const { return 0.0; }
	double yz() const { return 0.0; }
	double xz() const { return 0.0; }

protected:
	double	d[3];	// the diagonal elements

	friend class mat3d;
	friend class mat3ds;
	friend class mat3da;
};

inline mat3dd operator * (double a, const mat3dd& d) { return d*a; }

//-----------------------------------------------------------------------------
//! This class describes a symmetric 3D matrix of doubles

class mat3ds
{
protected:
	// This enumeration can be used to remember the order
	// in which the components are stored.
	enum {
		XX = 0,
		XY = 1,
		YY = 2,
		XZ = 3,
		YZ = 4,
		ZZ = 5 };
public:
	// default constructor
	mat3ds(){}

	// constructors
	explicit mat3ds(double a);
	mat3ds(double xx, double yy, double zz, double xy, double yz, double xz);
	mat3ds(const mat3dd& d);
	mat3ds(const mat3ds& d);

	// access operators
	double& operator () (int i, int j);
	const double& operator () (int i, int j) const;

	double& xx() { return m[XX]; }
	double& yy() { return m[YY]; }
	double& zz() { return m[ZZ]; }
	double& xy() { return m[XY]; }
	double& yz() { return m[YZ]; }
	double& xz() { return m[XZ]; }

	const double& xx() const { return m[XX]; }
	const double& yy() const { return m[YY]; }
	const double& zz() const { return m[ZZ]; }
	const double& xy() const { return m[XY]; }
	const double& yz() const { return m[YZ]; }
	const double& xz() const { return m[XZ]; }

	// arithmetic operators for mat3dd objects
	mat3ds operator + (const mat3dd& d) const;
	mat3ds operator - (const mat3dd& d) const;
	mat3ds operator * (const mat3dd& d) const;

	// arithmetic operators
	mat3ds operator + (const mat3ds& t) const;
	mat3ds operator - (const mat3ds& t) const;
	mat3d  operator * (const mat3ds& t) const;
	mat3ds operator * (double g) const;
	mat3ds operator / (double g) const;

	// arithmetic operators for mat3d objects
	mat3d operator + (const mat3d& t) const;
	mat3d operator - (const mat3d& t) const;
	mat3d operator * (const mat3d& t) const;
	
    // arithmetic operators for mat3da objects
    mat3d operator + (const mat3da& t) const;
    mat3d operator - (const mat3da& t) const;
    mat3d operator * (const mat3da& t) const;
    
	// unary operators
	mat3ds operator - () const;
	
	// arithmetic assignment operators
	mat3ds& operator += (const mat3ds& t);
	mat3ds& operator -= (const mat3ds& t);
	mat3ds& operator *= (const mat3ds& t);
	mat3ds& operator *= (double g);
	mat3ds& operator /= (double g);

	// arithmetic assignment operators for mat3dd
	mat3ds& operator += (const mat3dd& d);
	mat3ds& operator -= (const mat3dd& d);

	// matrix-vector multiplication
	vec3d operator * (const vec3d& r) const;

	// trace
	double tr() const;

	// determinant
	double det() const;

	// intialize to zero
	void zero();

	// initialize to unit tensor
	void unit();

	// deviator
	mat3ds dev() const;

	// isotropic part
	mat3ds iso() const;

	// return the square 
	mat3ds sqr() const;

	// calculates the inverse
	mat3ds inverse() const;
    double invert(mat3ds& Ai);
	
	// determine eigen values and vectors
	FECORE_API void eigen(double d[3], vec3d r[3] = 0) const;
	FECORE_API void exact_eigen(double l[3]) const;
	FECORE_API void eigen2(double d[3], vec3d r[3] = 0) const;

	// L2-norm 
	double norm() const;

	// double contraction
	double dotdot(const mat3ds& S) const;

	// "effective" or von-Mises norm
	double effective_norm() const;

	// the "max shear" value
	FECORE_API double max_shear() const;

protected:
	double m[6];	// stores data in the order xx, xy, yy, xz, yz, zz

	friend class mat3dd;
	friend class mat3d;
};

inline mat3ds operator * (double a, const mat3ds& m) { return m*a; }

//-----------------------------------------------------------------------------
//! This class describes an anti-symmetric 3D matrix of doubles
//! The matrix is defined such that for a vector b the following is true:
//! A.b = a x b where A = mat3da(a).
//!
//     | 0 -z  y |   |   0  d0  d2 |
// A = | z  0 -x | = | -d0   0  d1 |
//     |-y  x  0 |   | -d2 -d1   0 |
//

class mat3da
{
public:
	// default constructor
	mat3da(){}

	// constructors
	mat3da(double xy, double yz, double xz);

	// calculates the antisymmetric matrix from a vector
	// A.b = a x b where A = mat3da(a).
	explicit mat3da(const vec3d& a);

	// access operator
	double operator () (int i, int j) const;

	double& xy() { return d[0]; }
	double& yz() { return d[1]; }
	double& xz() { return d[2]; }

	const double& xy() const { return d[0]; }
	const double& yz() const { return d[1]; }
	const double& xz() const { return d[2]; }

	mat3da operator + (const mat3da& a);
	mat3da operator - (const mat3da& a);

	mat3da operator - () const;

	mat3da operator * (double g) const;

	mat3da transpose() const;

	// matrix algebra
	mat3d operator * (const mat3d& a);

	// return the equivalent vector
	vec3d vec() const { return vec3d(-d[1], d[2], -d[0]); }

	vec3d operator * (const vec3d& a);

    // arithmetic operators for mat3ds objects
    mat3d operator + (const mat3ds& t) const;
    mat3d operator - (const mat3ds& t) const;
    
protected:
	double	d[3];	// stores xy, yz, xz

	friend class mat3dd;
	friend class mat3ds;
	friend class mat3d;
};


//-----------------------------------------------------------------------------
//! This class describes a general 3D matrix of doubles
class mat3d
{
public:
	// default constructor
	mat3d() {}

	explicit mat3d(double a);

	// constructors
	mat3d(double a00, double a01, double a02,
		  double a10, double a11, double a12,
		  double a20, double a21, double a22);

	mat3d(double m[3][3]);
	mat3d(double a[9]);

	mat3d(const mat3dd& m);
	mat3d(const mat3ds& m);
	mat3d(const mat3da& m);

	mat3d(const mat2d& m);

	mat3d(const vec3d& e1, const vec3d& e2, const vec3d& e3);

	// assignment operators
	mat3d& operator = (const mat3dd& m);
	mat3d& operator = (const mat3ds& m);
	mat3d& operator = (const mat3d& m);
	mat3d& operator = (const double m[3][3]);

	// mat3d
	mat3d operator - () 
	{
		return mat3d(-d[0][0], -d[0][1], -d[0][2], \
					 -d[1][0], -d[1][1], -d[1][2], \
					 -d[2][0], -d[2][1], -d[2][2]);
	}

	// access operators
	double& operator () (int i, int j);
	const double& operator () (int i, int j) const;
	double* operator [] (int i);
	const double* operator [] (int i) const;

	// arithmetic operators
	mat3d operator + (const mat3d& m) const;
	mat3d operator - (const mat3d& m) const;
	mat3d operator * (const mat3d& m) const;
	mat3d operator * (double a) const;
	mat3d operator / (double a) const;

	// arithmetic operators for mat3dd
	mat3d operator + (const mat3dd& m) const;
	mat3d operator - (const mat3dd& m) const;
	mat3d operator * (const mat3dd& m) const;

	// arithmetic operators for mat3ds
	mat3d operator + (const mat3ds& m) const;
	mat3d operator - (const mat3ds& m) const;
	mat3d operator * (const mat3ds& m) const;

	// arithmetic assignment operators
	mat3d& operator += (const mat3d& m);
	mat3d& operator -= (const mat3d& m);
	mat3d& operator *= (const mat3d& m);
	mat3d& operator *= (double a);
	mat3d& operator /= (double a);

	// arithmetic assignment operators for mat3dd
	mat3d& operator += (const mat3dd& m);
	mat3d& operator -= (const mat3dd& m);
	mat3d& operator *= (const mat3dd& m);

	// arithmetic assignment operators for mat3ds
	mat3d& operator += (const mat3ds& m);
	mat3d& operator -= (const mat3ds& m);
	mat3d& operator *= (const mat3ds& m);

	// matrix-vector muliplication
	vec3d operator * (const vec3d& r) const;

	// determinant
	double det() const;

	// trace
	double trace() const;

	// zero the matrix
	void zero();

	// make unit matrix
	void unit();

	// return a column vector from the matrix
	vec3d col(int j) const;

	// return a row vector from the matrix
	vec3d row(int j) const;

	// set the column of the matrix
	void setCol(int i, const vec3d& a);

	// set the row of the matrix
	void setRow(int i, const vec3d& a);

	// return the symmetric matrix 0.5*(A+A^T)
	mat3ds sym() const;

	// return the antisymmetric matrix 0.5*(A-A^T)
	mat3da skew() const;

	// calculates the inverse
	mat3d inverse() const;
    
	// inverts the matrix.
	bool invert();

	// calculates the transpose
	mat3d transpose() const;

	// calculates the transposed inverse
	mat3d transinv() const;

	// calculate the skew-symmetric matrix from a vector
	void skew(const vec3d& v);

	// calculate the exponential map
	void exp(const vec3d& v);

	// calculate the one-norm
	double norm() const;

	// double contraction
	double dotdot(const mat3d& T) const;

	// polar decomposition
	FECORE_API void right_polar(mat3d& R, mat3ds& U) const;
	FECORE_API void left_polar(mat3ds& V, mat3d& R) const;

	// return identity matrix
	static mat3d identity() { return mat3d(1,0,0, 0,1,0, 0,0,1); }

protected:
	double d[3][3];	// matrix data

	friend class mat3dd;
	friend class mat3ds;
	friend class mat3da;
};

// outer product for vectors
inline mat3d operator & (const vec3d& a, const vec3d& b)
{
	return mat3d(a.x*b.x, a.x*b.y, a.x*b.z,
				 a.y*b.x, a.y*b.y, a.y*b.z,
				 a.z*b.x, a.z*b.y, a.z*b.z);
}

inline mat3ds dyad(const vec3d& a)
{
	return mat3ds(a.x*a.x, a.y*a.y, a.z*a.z, a.x*a.y, a.y*a.z, a.x*a.z);
}

// c_ij = a_i*b_j + a_j*b_i
inline mat3ds dyads(const vec3d& a, const vec3d& b)
{
	return mat3ds(2.0*a.x*b.x, 2.0*a.y*b.y, 2.0*a.z*b.z, a.x*b.y + a.y*b.x, a.y*b.z + a.z*b.y, a.x*b.z + a.z*b.x);
}

// skew-symmetric matrix of dual vector
inline mat3d skew(const vec3d& a)
{
    return mat3d(   0, -a.z,  a.y,
                  a.z,    0, -a.x,
                 -a.y,  a.x,    0);
}

//-----------------------------------------------------------------------------
// This class stores a 2nd order diagonal tensor
class mat3fd
{
public:
	mat3fd() { x = y = z = 0.f; }
	mat3fd(float X, float Y, float Z) { x = X; y = Y; z = Z; }

public:
	float x, y, z;
};

//-----------------------------------------------------------------------------
// mat3fs stores a 2nd order symmetric tensor
//
class mat3fs
{
public:
	// constructors
	mat3fs() { x = y = z = xy = yz = xz = 0; }
	mat3fs(float fx, float fy, float fz, float fxy, float fyz, float fxz)
	{
		x = fx; y = fy; z = fz;
		xy = fxy; yz = fyz; xz = fxz;
	}

	// operators
	mat3fs& operator += (const mat3fs& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		xy += v.xy;
		yz += v.yz;
		xz += v.xz;

		return (*this);
	}

	// operators
	mat3fs& operator -= (const mat3fs& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		xy -= v.xy;
		yz -= v.yz;
		xz -= v.xz;

		return (*this);
	}

	mat3fs& operator *= (float g)
	{
		x *= g;
		y *= g;
		z *= g;
		xy *= g;
		yz *= g;
		xz *= g;

		return (*this);
	}

	mat3fs& operator /= (float g)
	{
		x /= g;
		y /= g;
		z /= g;
		xy /= g;
		yz /= g;
		xz /= g;

		return (*this);
	}

	mat3fs operator + (const mat3fs& a) { return mat3fs(x + a.x, y + a.y, z + a.z, xy + a.xy, yz + a.yz, xz + a.xz); }
	mat3fs operator - (const mat3fs& a) { return mat3fs(x - a.x, y - a.y, z - a.z, xy - a.xy, yz - a.yz, xz - a.xz); }

	mat3fs operator * (float a)
	{
		return mat3fs(x * a, y * a, z * a, xy * a, yz * a, xz * a);
	}

	mat3fs operator / (float g)
	{
		return mat3fs(x / g, y / g, z / g, xy / g, yz / g, xz / g);
	}

	vec3f operator * (vec3f& r)
	{
		return vec3f(
			x * r.x + xy * r.y + xz * r.z,
			xy * r.x + y * r.y + yz * r.z,
			xz * r.x + yz * r.y + z * r.z);
	}

	// Effective or von-mises value
	float von_mises() const
	{
		float vm;
		vm = x * x + y * y + z * z;
		vm -= x * y + y * z + x * z;
		vm += 3 * (xy * xy + yz * yz + xz * xz);
		vm = (float)sqrt(vm >= 0.0 ? vm : 0.0);
		return vm;
	}

	// principle values
	FECORE_API void Principals(float e[3]) const;

	// principle directions
	FECORE_API vec3f PrincDirection(int l);

	// deviatroric principle values
	FECORE_API void DeviatoricPrincipals(float e[3]) const;

	// max-shear value
	FECORE_API float MaxShear() const;

	// eigen-vectors and values
	FECORE_API void eigen(vec3f e[3], float l[3]) const;

	// trace
	float tr() const { return x + y + z; }

	// determinant
	float det() const { return (x * y * z + xy * yz * xz + xz * xy * yz - y * xz * xz - x * yz * yz - z * xy * xy); }

	// L2 norm
	float norm() const {
		double d = x * x + y * y + z * z + 2 * (xy * xy + yz * yz + xz * xz);
		return (float)sqrt(d);
	}

public:
	float x, y, z;
	float xy, yz, xz;
};

FECORE_API double fractional_anisotropy(const mat3fs& m);

///////////////////////////////////////////////////////////////////
// mat3f

class mat3f
{
public:
	mat3f() { zero(); }

	mat3f(float a00, float a01, float a02, float a10, float a11, float a12, float a20, float a21, float a22)
	{
		d[0][0] = a00; d[0][1] = a01; d[0][2] = a02;
		d[1][0] = a10; d[1][1] = a11; d[1][2] = a12;
		d[2][0] = a20; d[2][1] = a21; d[2][2] = a22;
	}

	mat3f(const mat3fs& a)
	{
		d[0][0] = a.x; d[0][1] = a.xy; d[0][2] = a.xz;
		d[1][0] = a.xy; d[1][1] = a.y; d[1][2] = a.yz;
		d[2][0] = a.xz; d[2][1] = a.yz; d[2][2] = a.z;
	}

	float* operator [] (int i) { return d[i]; }
	float& operator () (int i, int j) { return d[i][j]; }
	float operator () (int i, int j) const { return d[i][j]; }

	mat3f operator * (float a) const
	{
		return mat3f(\
			d[0][0] * a, d[0][1] * a, d[0][2] * a, \
			d[1][0] * a, d[1][1] * a, d[1][2] * a, \
			d[2][0] * a, d[2][1] * a, d[2][2] * a);
	}

	mat3f operator * (mat3f& m)
	{
		mat3f a;

		int k;
		for (k = 0; k < 3; k++)
		{
			a[0][0] += d[0][k] * m[k][0]; a[0][1] += d[0][k] * m[k][1]; a[0][2] += d[0][k] * m[k][2];
			a[1][0] += d[1][k] * m[k][0]; a[1][1] += d[1][k] * m[k][1]; a[1][2] += d[1][k] * m[k][2];
			a[2][0] += d[2][k] * m[k][0]; a[2][1] += d[2][k] * m[k][1]; a[2][2] += d[2][k] * m[k][2];
		}

		return a;
	}

	vec3f operator * (const vec3f& a) const
	{
		return vec3f(
			d[0][0] * a.x + d[0][1] * a.y + d[0][2] * a.z,
			d[1][0] * a.x + d[1][1] * a.y + d[1][2] * a.z,
			d[2][0] * a.x + d[2][1] * a.y + d[2][2] * a.z
			);
	}

	mat3f& operator *= (float g)
	{
		d[0][0] *= g;	d[0][1] *= g; d[0][2] *= g;
		d[1][0] *= g;	d[1][1] *= g; d[1][2] *= g;
		d[2][0] *= g;	d[2][1] *= g; d[2][2] *= g;
		return (*this);
	}

	mat3f& operator /= (float g)
	{
		d[0][0] /= g;	d[0][1] /= g; d[0][2] /= g;
		d[1][0] /= g;	d[1][1] /= g; d[1][2] /= g;
		d[2][0] /= g;	d[2][1] /= g; d[2][2] /= g;
		return (*this);
	}

	mat3f operator + (const mat3f& a) const
	{
		return mat3f( \
			d[0][0] + a.d[0][0], d[0][1] + a.d[0][1], d[0][2] + a.d[0][2], \
			d[1][0] + a.d[1][0], d[1][1] + a.d[1][1], d[1][2] + a.d[1][2], \
			d[2][0] + a.d[2][0], d[2][1] + a.d[2][1], d[2][2] + a.d[2][2]);
	}

	mat3f operator - (const mat3f& a) const
	{
		return mat3f(\
			d[0][0] - a.d[0][0], d[0][1] - a.d[0][1], d[0][2] - a.d[0][2], \
			d[1][0] - a.d[1][0], d[1][1] - a.d[1][1], d[1][2] - a.d[1][2], \
			d[2][0] - a.d[2][0], d[2][1] - a.d[2][1], d[2][2] - a.d[2][2]);
	}

	mat3f operator += (const mat3f& a)
	{
		d[0][0] += a.d[0][0]; d[0][1] += a.d[0][1]; d[0][2] += a.d[0][2];
		d[1][0] += a.d[1][0]; d[1][1] += a.d[1][1]; d[1][2] += a.d[1][2];
		d[2][0] += a.d[2][0]; d[2][1] += a.d[2][1]; d[2][2] += a.d[2][2];
		return (*this);
	}

	mat3f operator -= (const mat3f& a)
	{
		d[0][0] -= a.d[0][0]; d[0][1] -= a.d[0][1]; d[0][2] -= a.d[0][2];
		d[1][0] -= a.d[1][0]; d[1][1] -= a.d[1][1]; d[1][2] -= a.d[1][2];
		d[2][0] -= a.d[2][0]; d[2][1] -= a.d[2][1]; d[2][2] -= a.d[2][2];
		return (*this);
	}

	mat3fs sym() const
	{
		return mat3fs(d[0][0], d[1][1], d[2][2], 0.5f * (d[0][1] + d[1][0]), 0.5f * (d[1][2] + d[2][1]), 0.5f * (d[0][2] + d[2][0]));
	}

	void zero()
	{
		d[0][0] = d[0][1] = d[0][2] = 0.f;
		d[1][0] = d[1][1] = d[1][2] = 0.f;
		d[2][0] = d[2][1] = d[2][2] = 0.f;
	}

	vec3f col(int i) const
	{
		vec3f r;
		switch (i)
		{
		case 0: r.x = d[0][0]; r.y = d[1][0]; r.z = d[2][0]; break;
		case 1: r.x = d[0][1]; r.y = d[1][1]; r.z = d[2][1]; break;
		case 2: r.x = d[0][2]; r.y = d[1][2]; r.z = d[2][2]; break;
		}
		return r;
	}

	vec3f row(int i) const
	{
		vec3f r;
		switch (i)
		{
		case 0: r.x = d[0][0]; r.y = d[0][1]; r.z = d[0][2]; break;
		case 1: r.x = d[1][0]; r.y = d[1][1]; r.z = d[1][2]; break;
		case 2: r.x = d[2][0]; r.y = d[2][1]; r.z = d[2][2]; break;
		}
		return r;
	}

	mat3f transpose() const
	{
		return mat3f(
			d[0][0], d[1][0], d[2][0],
			d[0][1], d[1][1], d[2][1],
			d[0][2], d[1][2], d[2][2]
		);
	}

	// inverts the matrix.
	bool invert();

public:
	float d[3][3];
};

// The following file contains the actual definition of the class functions
#include "mat3d.hpp"
