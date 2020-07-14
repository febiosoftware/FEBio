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
// NOTE: This file is automatically included from mat3d.h
// Users should not include this file manually!

//-----------------------------------------------------------------------------
// class mat3dd : class describing diagonal 3D matrices of doubles
//-----------------------------------------------------------------------------

// constructor
inline mat3dd::mat3dd(double a) { d[0] = d[1] = d[2] = a; }
inline mat3dd::mat3dd(double a0, double a1, double a2) { d[0] = a0; d[1] = a1; d[2] = a2; }

// assignment operators
inline mat3dd& mat3dd::operator = (const mat3dd& m) { d[0] = m.d[0]; d[1] = m.d[1]; d[2] = m.d[2]; return (*this); }
inline mat3dd& mat3dd::operator = (double a) { d[0] = d[1] = d[2] = a; return (*this); }

// access operators
inline double mat3dd::operator () (int i, int j) const { return (i==j? d[i] : 0); }
inline double& mat3dd::diag(int i) { return d[i]; }
inline const double& mat3dd::diag(int i) const { return d[i]; }

// arithmetic operators
inline mat3dd mat3dd::operator + (const mat3dd& m) const { return mat3dd(d[0]+m.d[0],d[1]+m.d[1],d[2]+m.d[2]);}
inline mat3dd mat3dd::operator - (const mat3dd& m) const { return mat3dd(d[0]-m.d[0],d[1]-m.d[1],d[2]-m.d[2]);}
inline mat3dd mat3dd::operator * (const mat3dd& m) const { return mat3dd(d[0]*m.d[0],d[1]*m.d[1],d[2]*m.d[2]);}
inline mat3dd mat3dd::operator * (double a) const { return mat3dd(d[0]*a, d[1]*a, d[2]*a); }
inline mat3dd mat3dd::operator / (double a) const { a = 1./a; return mat3dd(d[0]*a, d[1]*a, d[2]*a); }
inline mat3dd mat3dd::operator - () const { return mat3dd(-d[0], -d[1], -d[2]); }

// arithmetic operators with mat3ds
inline mat3ds mat3dd::operator + (const mat3ds& m) const
{
	return mat3ds(
		d[0]+m.m[mat3ds::XX],
		d[1]+m.m[mat3ds::YY],
		d[2]+m.m[mat3ds::ZZ],
		m.m[mat3ds::XY],
		m.m[mat3ds::YZ],
		m.m[mat3ds::XZ]
		);
}

inline mat3ds mat3dd::operator - (const mat3ds& m) const
{
	return mat3ds(
		d[0]-m.m[mat3ds::XX],
		d[1]-m.m[mat3ds::YY],
		d[2]-m.m[mat3ds::ZZ],
		-m.m[mat3ds::XY],
		-m.m[mat3ds::YZ],
		-m.m[mat3ds::XZ]
		);
}


inline mat3ds mat3dd::operator * (const mat3ds& m) const
{
	return mat3ds(
		d[0]*m.m[mat3ds::XX],
		d[1]*m.m[mat3ds::YY],
		d[2]*m.m[mat3ds::ZZ],
		d[0]*m.m[mat3ds::XY],
		d[1]*m.m[mat3ds::YZ],
		d[0]*m.m[mat3ds::XZ]
	);
}

// arithmetic operators for mat3d
inline mat3d mat3dd::operator + (const mat3d& m) const
{
	return mat3d(d[0]+m.d[0][0], m.d[0][1], m.d[0][2],
				 m.d[1][0], d[1]+m.d[1][1], m.d[1][2],
				 m.d[2][0], m.d[2][1], d[2]+m.d[2][2]);
}

inline mat3d mat3dd::operator - (const mat3d& m) const
{
	return mat3d(d[0]-m.d[0][0], -m.d[0][1], -m.d[0][2],
				 -m.d[1][0], d[1]-m.d[1][1], -m.d[1][2],
				 -m.d[2][0], -m.d[2][1], d[2]-m.d[2][2]);
}

inline mat3d mat3dd::operator * (const mat3d& m) const
{
	return mat3d(d[0]*m.d[0][0], d[0]*m.d[0][1], d[0]*m.d[0][2],
				 d[1]*m.d[1][0], d[1]*m.d[1][1], d[1]*m.d[1][2],
				 d[2]*m.d[2][0], d[2]*m.d[2][1], d[2]*m.d[2][2]);
}

// arithmetic operators for mat3da
inline mat3d mat3dd::operator + (const mat3da& m) const
{
	return mat3d(   d[0],  m.d[0], m.d[2],
				 -m.d[0],    d[1], m.d[1],
				 -m.d[2], -m.d[1],   d[2]);
}

inline mat3d mat3dd::operator - (const mat3da& m) const
{
	return mat3d(  d[0],-m.d[0], -m.d[2],
				 m.d[0],   d[1], -m.d[1],
				 m.d[2], m.d[1],   d[2]);
}

inline mat3d mat3dd::operator *(const mat3da& m) const
{
	return mat3d(           0,  d[0]*m.d[0], d[0]*m.d[2],
				 -d[1]*m.d[0],            0, d[1]*m.d[1],
				 -d[2]*m.d[2], -d[2]*m.d[1],           0);
}

// arithmetic assignment operators
inline mat3dd& mat3dd::operator += (const mat3dd& m) { d[0] += m.d[0]; d[1] += m.d[1]; d[2] += m.d[2]; return (*this); }
inline mat3dd& mat3dd::operator -= (const mat3dd& m) { d[0] -= m.d[0]; d[1] -= m.d[1]; d[2] -= m.d[2]; return (*this); }
inline mat3dd& mat3dd::operator *= (const mat3dd& m) { d[0] *= m.d[0]; d[1] *= m.d[1]; d[2] *= m.d[2]; return (*this); }
inline mat3dd& mat3dd::operator *= (double a) { d[0] *= a; d[1] *= a; d[2] *= a; return (*this); }
inline mat3dd& mat3dd::operator /= (double a) { a = 1./a; d[0] *= a; d[1] *= a; d[2] *= a; return (*this); }

// matrix-vector multiplication
inline vec3d mat3dd::operator * (const vec3d& r) const { return vec3d(r.x*d[0], r.y*d[1], r.z*d[2]); }

// trace
inline double mat3dd::tr() const { return d[0]+d[1]+d[2]; }

// determinant
inline double mat3dd::det() const { return d[0]*d[1]*d[2]; }

//-----------------------------------------------------------------------------
// class mat3ds : this class describes a symmetric 3D matrix of doubles
//-----------------------------------------------------------------------------

// constructor
inline mat3ds::mat3ds(double a)
{
	m[XX] = a;
	m[XY] = a;
	m[YY] = a;
	m[XZ] = a;
	m[YZ] = a;
	m[ZZ] = a;
}

inline mat3ds::mat3ds(double xx, double yy, double zz, double xy, double yz, double xz)
{
	m[XX] = xx;
	m[XY] = xy;
	m[YY] = yy;
	m[XZ] = xz;
	m[YZ] = yz;
	m[ZZ] = zz;
}

inline mat3ds::mat3ds(const mat3dd& d)
{
	m[XX] = d.d[0];
	m[YY] = d.d[1];
	m[ZZ] = d.d[2];
	m[XY] = m[YZ] = m[XZ] = 0.;
}

inline mat3ds::mat3ds(const mat3ds& d)
{
	m[0] = d.m[0];
	m[1] = d.m[1];
	m[2] = d.m[2];
	m[3] = d.m[3];
	m[4] = d.m[4];
	m[5] = d.m[5];
}

// access operator
inline double& mat3ds::operator ()(int i, int j)
{
	const int n[] = {0, 1, 3};
	return (i<=j? m[n[j]+i] : m[n[i]+j]);
}

// access operator for const objects
inline const double& mat3ds::operator ()(int i, int j) const
{
	const int n[] = {0, 1, 3};
	return (i<=j? m[n[j]+i] : m[n[i]+j]);
}

// operator + for mat3dd objects
inline mat3ds mat3ds::operator + (const mat3dd& d) const
{
	return mat3ds(m[XX]+d.d[0], m[YY]+d.d[1], m[ZZ]+d.d[2], m[XY], m[YZ], m[XZ]);
}

// operator - for mat3dd objects
inline mat3ds mat3ds::operator - (const mat3dd& d) const
{
	return mat3ds(m[XX]-d.d[0], m[YY]-d.d[1], m[ZZ]-d.d[2], m[XY], m[YZ], m[XZ]);
}

// operator * for mat3dd objects
inline mat3ds mat3ds::operator * (const mat3dd& d) const
{
	return mat3ds(m[XX]*d.d[0], m[YY]*d.d[1], m[ZZ]*d.d[2], m[XY]*d.d[1], m[YZ]*d.d[2], m[XZ]*d.d[2]);
}


// operator +
inline mat3ds mat3ds::operator +(const mat3ds& t) const
{
	return mat3ds(m[XX]+t.m[XX], m[YY]+t.m[YY], m[ZZ]+t.m[ZZ], m[XY]+t.m[XY], m[YZ]+t.m[YZ],m[XZ]+t.m[XZ]);
}

// operator -
inline mat3ds mat3ds::operator -(const mat3ds& t) const
{
	return mat3ds(m[XX]-t.m[XX], m[YY]-t.m[YY], m[ZZ]-t.m[ZZ], m[XY]-t.m[XY], m[YZ]-t.m[YZ],m[XZ]-t.m[XZ]);
}

// operator *
inline mat3d mat3ds::operator *(const mat3ds& t) const
{
	return mat3d(
		m[XX] * t.m[XX] + m[XY] * t.m[XY] + m[XZ] * t.m[XZ],
		m[XX] * t.m[XY] + m[XY] * t.m[YY] + m[XZ] * t.m[YZ],
		m[XX] * t.m[XZ] + m[XY] * t.m[YZ] + m[XZ] * t.m[ZZ],

		m[XY] * t.m[XX] + m[YY] * t.m[XY] + m[YZ] * t.m[XZ],
		m[XY] * t.m[XY] + m[YY] * t.m[YY] + m[YZ] * t.m[YZ],
		m[XY] * t.m[XZ] + m[YY] * t.m[YZ] + m[YZ] * t.m[ZZ],

		m[XZ] * t.m[XX] + m[YZ] * t.m[XY] + m[ZZ] * t.m[XZ],
		m[XZ] * t.m[XY] + m[YZ] * t.m[YY] + m[ZZ] * t.m[YZ],
		m[XZ] * t.m[XZ] + m[YZ] * t.m[YZ] + m[ZZ] * t.m[ZZ]
	);
}

// operator *
inline mat3ds mat3ds::operator * (double g) const
{
	return mat3ds(m[XX]*g, m[YY]*g, m[ZZ]*g, m[XY]*g, m[YZ]*g, m[XZ]*g);
}

// operator /
inline mat3ds mat3ds::operator / (double g) const
{
	g = 1.0/g;
	return mat3ds(m[XX]*g, m[YY]*g, m[ZZ]*g, m[XY]*g, m[YZ]*g, m[XZ]*g);
}

// operator + for mat3d objects
inline mat3d mat3ds::operator + (const mat3d& d) const
{
	return mat3d(m[XX]+d.d[0][0], m[XY]+d.d[0][1], m[XZ]+d.d[0][2],
				 m[XY]+d.d[1][0], m[YY]+d.d[1][1], m[YZ]+d.d[1][2],
				 m[XZ]+d.d[2][0], m[YZ]+d.d[2][1], m[ZZ]+d.d[2][2]);
}

// operator - for mat3d objects
inline mat3d mat3ds::operator - (const mat3d& d) const
{
	return mat3d(m[XX]-d.d[0][0], m[XY]-d.d[0][1], m[XZ]-d.d[0][2],
				 m[XY]-d.d[1][0], m[YY]-d.d[1][1], m[YZ]-d.d[1][2],
				 m[XZ]-d.d[2][0], m[YZ]-d.d[2][1], m[ZZ]-d.d[2][2]);
}

// operator * for mat3d objects
inline mat3d mat3ds::operator * (const mat3d& d) const
{
	return mat3d(d.d[0][0]*m[XX] + d.d[1][0]*m[XY] + d.d[2][0]*m[XZ], 
				 d.d[0][1]*m[XX] + d.d[1][1]*m[XY] + d.d[2][1]*m[XZ], 
				 d.d[0][2]*m[XX] + d.d[1][2]*m[XY] + d.d[2][2]*m[XZ],
				 d.d[0][0]*m[XY] + d.d[1][0]*m[YY] + d.d[2][0]*m[YZ], 
				 d.d[0][1]*m[XY] + d.d[1][1]*m[YY] + d.d[2][1]*m[YZ], 
				 d.d[0][2]*m[XY] + d.d[1][2]*m[YY] + d.d[2][2]*m[YZ],
				 d.d[0][0]*m[XZ] + d.d[1][0]*m[YZ] + d.d[2][0]*m[ZZ], 
				 d.d[0][1]*m[XZ] + d.d[1][1]*m[YZ] + d.d[2][1]*m[ZZ], 
				 d.d[0][2]*m[XZ] + d.d[1][2]*m[YZ] + d.d[2][2]*m[ZZ]);
}


// unary operator -
inline mat3ds mat3ds::operator - () const
{
	return mat3ds(-m[XX], -m[YY], -m[ZZ], -m[XY], -m[YZ], -m[XZ]);
}

// assignment operator +=
inline mat3ds& mat3ds::operator += (const mat3ds& t)
{
	m[XX] += t.m[XX]; m[YY] += t.m[YY]; m[ZZ] += t.m[ZZ];
	m[XY] += t.m[XY]; m[YZ] += t.m[YZ]; m[XZ] += t.m[XZ];
	return (*this);
}

// assignment operator -=
inline mat3ds& mat3ds::operator -= (const mat3ds& t)
{
	m[XX] -= t.m[XX]; m[YY] -= t.m[YY]; m[ZZ] -= t.m[ZZ];
	m[XY] -= t.m[XY]; m[YZ] -= t.m[YZ]; m[XZ] -= t.m[XZ];
	return (*this);
}

// assignment operator *=
inline mat3ds& mat3ds::operator *= (const mat3ds& t)
{
	double xx = m[XX]*t.m[XX]+m[XY]*t.m[XY]+m[XZ]*t.m[XZ];
	double yy = m[XY]*t.m[XY]+m[YY]*t.m[YY]+m[YZ]*t.m[YZ];
	double zz = m[XZ]*t.m[XZ]+m[YZ]*t.m[YZ]+m[ZZ]*t.m[ZZ];
	double xy = m[XX]*t.m[XY]+m[XY]*t.m[YY]+m[XZ]*t.m[YZ];
	double yz = m[XY]*t.m[XZ]+m[YY]*t.m[YZ]+m[YZ]*t.m[ZZ];
	double xz = m[XX]*t.m[XZ]+m[XY]*t.m[YZ]+m[XZ]*t.m[ZZ];

	m[XX] = xx; m[YY] = yy; m[ZZ] = zz;
	m[XY] = xy; m[YZ] = yz; m[XZ] = xz;

	return (*this);
}

// assignment operator *=
inline mat3ds& mat3ds::operator *= (double g)
{
	m[XX] *= g; m[YY] *= g; m[ZZ] *= g;
	m[XY] *= g; m[YZ] *= g; m[XZ] *= g;
	return (*this);
}

// assignment operator /=
inline mat3ds& mat3ds::operator /= (double g)
{
	g = 1./g;
	m[XX] *= g; m[YY] *= g; m[ZZ] *= g;
	m[XY] *= g; m[YZ] *= g; m[XZ] *= g;
	return (*this);
}

// arithmetic assignment operators for mat3dd
inline mat3ds& mat3ds::operator += (const mat3dd& d)
{
	m[XX] += d.d[0];
	m[YY] += d.d[1];
	m[ZZ] += d.d[2];
	return (*this);
}

inline mat3ds& mat3ds::operator -= (const mat3dd& d)
{
	m[XX] -= d.d[0];
	m[YY] -= d.d[1];
	m[ZZ] -= d.d[2];
	return (*this);
}

// matrix-vector multiplication
inline vec3d mat3ds::operator* (const vec3d& r) const
{
	return vec3d(
		r.x*m[XX]+r.y*m[XY]+r.z*m[XZ],
		r.x*m[XY]+r.y*m[YY]+r.z*m[YZ],
		r.x*m[XZ]+r.y*m[YZ]+r.z*m[ZZ]
	);
}

// trace
inline double mat3ds::tr() const
{
	return m[XX]+m[YY]+m[ZZ];
}

// determinant
inline double mat3ds::det() const
{
	return (m[XX]*(m[YY]*m[ZZ] - m[YZ]*m[YZ])
		  + m[XY]*(m[YZ]*m[XZ] - m[ZZ]*m[XY])
		  + m[XZ]*(m[XY]*m[YZ] - m[YY]*m[XZ]));
}

// zero
inline void mat3ds::zero()
{
	m[0] = m[1] = m[2] = m[3] = m[4] = m[5] = 0;
}

// unit tensor
inline void mat3ds::unit()
{
	m[XX] = m[YY] = m[ZZ] = 1.0;
	m[XY] = m[YZ] = m[XZ] = 0.0;
}

// deviator
inline mat3ds mat3ds::dev() const
{
	double t = (m[XX]+m[YY]+m[ZZ])/3.0;
	return mat3ds(m[XX]-t, m[YY]-t, m[ZZ]-t, m[XY], m[YZ], m[XZ]);
}

// isotropic part
inline mat3ds mat3ds::iso() const
{
	double t = (m[XX]+m[YY]+m[ZZ])/3.0;
	return mat3ds(t, t, t, 0, 0, 0);
}

// return the square 
inline mat3ds mat3ds::sqr() const
{
	return mat3ds(
		m[XX] * m[XX] + m[XY] * m[XY] + m[XZ] * m[XZ],
		m[XY] * m[XY] + m[YY] * m[YY] + m[YZ] * m[YZ],
		m[XZ] * m[XZ] + m[YZ] * m[YZ] + m[ZZ] * m[ZZ],
		m[XX] * m[XY] + m[XY] * m[YY] + m[XZ] * m[YZ], 
		m[XY] * m[XZ] + m[YY] * m[YZ] + m[YZ] * m[ZZ],
		m[XX] * m[XZ] + m[XY] * m[YZ] + m[XZ] * m[ZZ]
	);
}

// inverse
inline mat3ds mat3ds::inverse() const
{
	double D = det();
	assert(D != 0);
	D = 1/D;
	
	return mat3ds(D*(m[YY]*m[ZZ]-m[YZ]*m[YZ]), 
				  D*(m[XX]*m[ZZ]-m[XZ]*m[XZ]), 
				  D*(m[XX]*m[YY]-m[XY]*m[XY]),
				  D*(m[XZ]*m[YZ]-m[XY]*m[ZZ]), 
				  D*(m[XY]*m[XZ]-m[XX]*m[YZ]), 
				  D*(m[XY]*m[YZ]-m[YY]*m[XZ]));
}

// invert
inline double mat3ds::invert(mat3ds& Ai)
{
    double D = det();
    if (D != 0) {
        double Di = 1/D;
        
        Ai = mat3ds(Di*(m[YY]*m[ZZ]-m[YZ]*m[YZ]),
                    Di*(m[XX]*m[ZZ]-m[XZ]*m[XZ]),
                    Di*(m[XX]*m[YY]-m[XY]*m[XY]),
                    Di*(m[XZ]*m[YZ]-m[XY]*m[ZZ]),
                    Di*(m[XY]*m[XZ]-m[XX]*m[YZ]),
                    Di*(m[XY]*m[YZ]-m[YY]*m[XZ]));
    }
    return D;
}

// L2-norm
inline double mat3ds::norm() const
{ 
	double D = m[XX]*m[XX] + m[YY]*m[YY] + m[ZZ]*m[ZZ] + 2*(m[XY]*m[XY] + m[YZ]*m[YZ] + m[XZ]*m[XZ]);
	return sqrt(D); 
}

// double contraction
inline double mat3ds::dotdot(const mat3ds& B) const
{
	const double* n = B.m;
	return m[XX]*n[XX] + m[YY]*n[YY] + m[ZZ]*n[ZZ] + 2.0*(m[XY]*n[XY] + m[YZ]*n[YZ] + m[XZ]*n[XZ]);
}

// Effective or von-mises value
inline double mat3ds::effective_norm() const
{
	double vm;
	vm = m[XX] * m[XX] + m[YY] * m[YY] + m[ZZ] * m[ZZ];
	vm -= m[XX] * m[YY] + m[YY] * m[ZZ] + m[XX] * m[ZZ];
	vm += 3 * (m[XY] * m[XY] + m[YZ] * m[YZ] + m[XZ] * m[XZ]);
	vm = sqrt(vm >= 0.0 ? vm : 0.0);
	return vm;
}

//-----------------------------------------------------------------------------
// class mat3da : anti-symmetric 3D matrix of doubles
//-----------------------------------------------------------------------------

// constructors
inline mat3da::mat3da(double xy, double yz, double xz)
{
	d[0] = xy; d[1] = yz; d[2] = xz;
}

// calculates the antisymmetric matrix from a vector such that for any b,
// A.b = a x b where A = mat3da(a).
inline mat3da::mat3da(const vec3d& a)
{
	d[0] = -a.z; d[1] = -a.x; d[2] = a.y;
}

// access operator
inline double mat3da::operator ()(int i, int j) const
{
	return (i==j? 0 : (i<j? d[((j-1)<<1)-i] : -d[((i-1)<<1)-j]));
}

// scalar multiplication
inline mat3da mat3da::operator * (double g) const
{
	return mat3da(d[0]*g, d[1]*g, d[2]*g);
}

inline mat3da mat3da::operator + (const mat3da& a)
{
	return mat3da(d[0]+a.d[0], d[1]+a.d[1], d[2]+a.d[2]);
}

inline mat3da mat3da::operator - (const mat3da& a)
{
	return mat3da(d[0]-a.d[0], d[1]-a.d[1], d[2]-a.d[2]);
}

inline mat3da mat3da::operator - () const
{
	return mat3da(-d[0], -d[1], -d[2]);
}

inline mat3da mat3da::transpose() const
{
	return mat3da(-d[0], -d[1], -d[2]);
}

// matrix multiplication
inline mat3d mat3da::operator * (const mat3d& m)
{
	return mat3d(
		 d[0]*m.d[1][0] + d[2]*m.d[2][0],  d[0]*m.d[1][1] + d[2]*m.d[2][1],  d[0]*m.d[1][2] + d[2]*m.d[2][2],
		-d[0]*m.d[0][0] + d[1]*m.d[2][0], -d[0]*m.d[0][1] + d[1]*m.d[2][1], -d[0]*m.d[0][2] + d[1]*m.d[2][2],
		-d[2]*m.d[0][0] - d[1]*m.d[1][0], -d[2]*m.d[0][1] - d[1]*m.d[1][1], -d[2]*m.d[0][2] - d[1]*m.d[1][2]
	);
}

inline vec3d mat3da::operator * (const vec3d& a)
{
	return vec3d(
		 d[0] * a.y + d[2] * a.z, \
		-d[0] * a.x + d[1] * a.z, \
		-d[2] * a.x - d[1] * a.y);
}

//-----------------------------------------------------------------------------
// class mat3d : general 3D matrix of doubles
//-----------------------------------------------------------------------------

// constructors

inline mat3d::mat3d(double a)
{
	d[0][0] = a; d[0][1] = a; d[0][2] = a;
	d[1][0] = a; d[1][1] = a; d[1][2] = a;
	d[2][0] = a; d[2][1] = a; d[2][2] = a;
}

inline mat3d::mat3d(double a00, double a01, double a02,
					double a10, double a11, double a12,
					double a20, double a21, double a22)
{
	d[0][0] = a00; d[0][1] = a01; d[0][2] = a02;
	d[1][0] = a10; d[1][1] = a11; d[1][2] = a12;
	d[2][0] = a20; d[2][1] = a21; d[2][2] = a22;
}

inline mat3d::mat3d(double m[3][3])
{
	d[0][0] = m[0][0]; d[0][1] = m[0][1]; d[0][2] = m[0][2];
	d[1][0] = m[1][0]; d[1][1] = m[1][1]; d[1][2] = m[1][2];
	d[2][0] = m[2][0]; d[2][1] = m[2][1]; d[2][2] = m[2][2];
}

inline mat3d::mat3d(const mat3dd& m)
{
	d[0][0] = m.d[0]; d[1][1] = m.d[1]; d[2][2] = m.d[2];
	d[0][1] = d[1][0] = 0;
	d[1][2] = d[2][1] = 0;
	d[0][2] = d[2][0] = 0;
}

inline mat3d::mat3d(const mat3ds& m)
{
	d[0][0] = m.m[mat3ds::XX];
	d[1][1] = m.m[mat3ds::YY];
	d[2][2] = m.m[mat3ds::ZZ];
	d[0][1] = d[1][0] = m.m[mat3ds::XY];
	d[1][2] = d[2][1] = m.m[mat3ds::YZ];
	d[0][2] = d[2][0] = m.m[mat3ds::XZ];
}

inline mat3d::mat3d(const mat3da& m)
{
	d[0][0] = d[1][1] = d[2][2] = 0;
	d[0][1] = m.d[0]; d[1][0] = -m.d[0];
	d[1][2] = m.d[1]; d[2][1] = -m.d[1];
	d[0][2] = m.d[2]; d[2][0] = -m.d[2];
}

inline mat3d::mat3d(const mat2d& m)
{
	d[0][0] = m(0,0); d[0][1] = m(0,1); d[0][2] = 0.0;
	d[1][0] = m(1,0); d[1][1] = m(1,1); d[1][2] = 0.0;
	d[2][0] = d[2][1] = 0.0;
	d[2][2] = 0.0;	// Should I set this to 1.0 instead? that way det(), inverse() etc. remain valid for the mat3d
}

inline mat3d::mat3d(const vec3d& e1, const vec3d& e2, const vec3d& e3)
{
	d[0][0] = e1.x; d[0][1] = e2.x; d[0][2] = e3.x;
	d[1][0] = e1.y; d[1][1] = e2.y; d[1][2] = e3.y;
	d[2][0] = e1.z; d[2][1] = e2.z; d[2][2] = e3.z;
}

// assignment operators
inline mat3d& mat3d::operator = (const mat3dd& m)
{
	d[0][0] = m.d[0];
	d[1][1] = m.d[1];
	d[2][2] = m.d[2];
	d[0][1] = d[1][0] = 0;
	d[1][2] = d[2][1] = 0;
	d[0][2] = d[2][0] = 0;
	return (*this);
}

inline mat3d& mat3d::operator = (const mat3ds& m)
{
	d[0][0] = m.m[mat3ds::XX];
	d[1][1] = m.m[mat3ds::YY];
	d[2][2] = m.m[mat3ds::ZZ];
	d[0][1] = d[1][0] = m.m[mat3ds::XY];
	d[1][2] = d[2][1] = m.m[mat3ds::YZ];
	d[0][2] = d[2][0] = m.m[mat3ds::XZ];
	return (*this);
}

inline mat3d& mat3d::operator = (const double m[3][3])
{
	d[0][0] = m[0][0]; d[0][1] = m[0][1]; d[0][2] = m[0][2];
	d[1][0] = m[1][0]; d[1][1] = m[1][1]; d[1][2] = m[1][2];
	d[2][0] = m[2][0]; d[2][1] = m[2][1]; d[2][2] = m[2][2];
	return (*this);
}

inline mat3d& mat3d::operator = (const mat3d& m)
{
	d[0][0] = m.d[0][0]; d[0][1] = m.d[0][1]; d[0][2] = m.d[0][2];
	d[1][0] = m.d[1][0]; d[1][1] = m.d[1][1]; d[1][2] = m.d[1][2];
	d[2][0] = m.d[2][0]; d[2][1] = m.d[2][1]; d[2][2] = m.d[2][2];
	return (*this);
}

// access operator
inline double& mat3d::operator () (int i, int j) { return d[i][j]; }
inline const double& mat3d::operator () (int i, int j) const { return d[i][j]; }
inline double* mat3d::operator [] (int i) { return d[i]; }
inline const double* mat3d::operator [] (int i) const { return d[i]; }

// arithmetic operators
inline mat3d mat3d::operator + (const mat3d& m) const
{
	return mat3d( d[0][0]+m.d[0][0], d[0][1]+m.d[0][1], d[0][2]+m.d[0][2],
				  d[1][0]+m.d[1][0], d[1][1]+m.d[1][1], d[1][2]+m.d[1][2],
				  d[2][0]+m.d[2][0], d[2][1]+m.d[2][1], d[2][2]+m.d[2][2]);
}

inline mat3d mat3d::operator - (const mat3d& m) const
{
	return mat3d( d[0][0]-m.d[0][0], d[0][1]-m.d[0][1], d[0][2]-m.d[0][2],
				  d[1][0]-m.d[1][0], d[1][1]-m.d[1][1], d[1][2]-m.d[1][2],
				  d[2][0]-m.d[2][0], d[2][1]-m.d[2][1], d[2][2]-m.d[2][2]);
}

inline mat3d mat3d::operator * (const mat3d& m) const
{
	return mat3d(d[0][0]*m.d[0][0]+d[0][1]*m.d[1][0]+d[0][2]*m.d[2][0],
				 d[0][0]*m.d[0][1]+d[0][1]*m.d[1][1]+d[0][2]*m.d[2][1],
				 d[0][0]*m.d[0][2]+d[0][1]*m.d[1][2]+d[0][2]*m.d[2][2],
				 d[1][0]*m.d[0][0]+d[1][1]*m.d[1][0]+d[1][2]*m.d[2][0],
				 d[1][0]*m.d[0][1]+d[1][1]*m.d[1][1]+d[1][2]*m.d[2][1],
				 d[1][0]*m.d[0][2]+d[1][1]*m.d[1][2]+d[1][2]*m.d[2][2],
				 d[2][0]*m.d[0][0]+d[2][1]*m.d[1][0]+d[2][2]*m.d[2][0],
				 d[2][0]*m.d[0][1]+d[2][1]*m.d[1][1]+d[2][2]*m.d[2][1],
				 d[2][0]*m.d[0][2]+d[2][1]*m.d[1][2]+d[2][2]*m.d[2][2]);
}

inline mat3d mat3d::operator * (double a) const
{
	return mat3d(d[0][0]*a, d[0][1]*a, d[0][2]*a,
				 d[1][0]*a, d[1][1]*a, d[1][2]*a,
				 d[2][0]*a, d[2][1]*a, d[2][2]*a);
}

inline mat3d mat3d::operator / (double a) const
{
	a = 1./a;
	return mat3d(d[0][0]*a, d[0][1]*a, d[0][2]*a,
				 d[1][0]*a, d[1][1]*a, d[1][2]*a,
				 d[2][0]*a, d[2][1]*a, d[2][2]*a);
}

// arithmetic operators for mat3dd
inline mat3d mat3d::operator + (const mat3dd& m) const
{
	return mat3d( d[0][0]+m.d[0], d[0][1], d[0][2],
				  d[1][0], d[1][1]+m.d[1], d[1][2],
				  d[2][0], d[2][1], d[2][2]+m.d[2]);
}

inline mat3d mat3d::operator - (const mat3dd& m) const
{
	return mat3d( d[0][0]-m.d[0], d[0][1], d[0][2],
				  d[1][0], d[1][1]-m.d[1], d[1][2],
				  d[2][0], d[2][1], d[2][2]-m.d[2]);
}

inline mat3d mat3d::operator * (const mat3dd& m) const
{
	return mat3d( d[0][0]*m.d[0], d[0][1]*m.d[1], d[0][2]*m.d[2],
				  d[1][0]*m.d[0], d[1][1]*m.d[1], d[1][2]*m.d[2],
				  d[2][0]*m.d[0], d[2][1]*m.d[1], d[2][2]*m.d[2]);
}

// arithmetic operators for mat3ds
inline mat3d mat3d::operator + (const mat3ds& m) const
{
	return mat3d(d[0][0]+m.m[m.XX], d[0][1]+m.m[m.XY], d[0][2]+m.m[m.XZ],
				 d[1][0]+m.m[m.XY], d[1][1]+m.m[m.YY], d[1][2]+m.m[m.YZ],
				 d[2][0]+m.m[m.XZ], d[2][1]+m.m[m.YZ], d[2][2]+m.m[m.ZZ]);
}

inline mat3d mat3d::operator - (const mat3ds& m) const
{
	return mat3d(d[0][0]-m.m[m.XX], d[0][1]-m.m[m.XY], d[0][2]-m.m[m.XZ],
				 d[1][0]-m.m[m.XY], d[1][1]-m.m[m.YY], d[1][2]-m.m[m.YZ],
				 d[2][0]-m.m[m.XZ], d[2][1]-m.m[m.YZ], d[2][2]-m.m[m.ZZ]);
}

inline mat3d mat3d::operator * (const mat3ds& m) const
{
	return mat3d(
		d[0][0]*m.m[m.XX] + d[0][1]*m.m[m.XY] + d[0][2]*m.m[m.XZ],
		d[0][0]*m.m[m.XY] + d[0][1]*m.m[m.YY] + d[0][2]*m.m[m.YZ],
		d[0][0]*m.m[m.XZ] + d[0][1]*m.m[m.YZ] + d[0][2]*m.m[m.ZZ],
		d[1][0]*m.m[m.XX] + d[1][1]*m.m[m.XY] + d[1][2]*m.m[m.XZ],
		d[1][0]*m.m[m.XY] + d[1][1]*m.m[m.YY] + d[1][2]*m.m[m.YZ],
		d[1][0]*m.m[m.XZ] + d[1][1]*m.m[m.YZ] + d[1][2]*m.m[m.ZZ],
		d[2][0]*m.m[m.XX] + d[2][1]*m.m[m.XY] + d[2][2]*m.m[m.XZ],
		d[2][0]*m.m[m.XY] + d[2][1]*m.m[m.YY] + d[2][2]*m.m[m.YZ],
		d[2][0]*m.m[m.XZ] + d[2][1]*m.m[m.YZ] + d[2][2]*m.m[m.ZZ]);
}

// arithmetic assignment operators
inline mat3d& mat3d::operator += (const mat3d& m)
{
	d[0][0] += m.d[0][0]; d[0][1] += m.d[0][1]; d[0][2] += m.d[0][2];
	d[1][0] += m.d[1][0]; d[1][1] += m.d[1][1]; d[1][2] += m.d[1][2];
	d[2][0] += m.d[2][0]; d[2][1] += m.d[2][1]; d[2][2] += m.d[2][2];
	return (*this);
}

inline mat3d& mat3d::operator -= (const mat3d& m)
{
	d[0][0] -= m.d[0][0]; d[0][1] -= m.d[0][1]; d[0][2] -= m.d[0][2];
	d[1][0] -= m.d[1][0]; d[1][1] -= m.d[1][1]; d[1][2] -= m.d[1][2];
	d[2][0] -= m.d[2][0]; d[2][1] -= m.d[2][1]; d[2][2] -= m.d[2][2];
	return (*this);
}

inline mat3d& mat3d::operator *= (const mat3d& m)
{
	double d00 = d[0][0]*m.d[0][0]+d[0][1]*m.d[1][0]+d[0][2]*m.d[2][0];
	double d01 = d[0][0]*m.d[0][1]+d[0][1]*m.d[1][1]+d[0][2]*m.d[2][1];
	double d02 = d[0][0]*m.d[0][2]+d[0][1]*m.d[1][2]+d[0][2]*m.d[2][2];
	double d10 = d[1][0]*m.d[0][0]+d[1][1]*m.d[1][0]+d[1][2]*m.d[2][0];
	double d11 = d[1][0]*m.d[0][1]+d[1][1]*m.d[1][1]+d[1][2]*m.d[2][1];
	double d12 = d[1][0]*m.d[0][2]+d[1][1]*m.d[1][2]+d[1][2]*m.d[2][2];
	double d20 = d[2][0]*m.d[0][0]+d[2][1]*m.d[1][0]+d[2][2]*m.d[2][0];
	double d21 = d[2][0]*m.d[0][1]+d[2][1]*m.d[1][1]+d[2][2]*m.d[2][1];
	double d22 = d[2][0]*m.d[0][2]+d[2][1]*m.d[1][2]+d[2][2]*m.d[2][2];

	d[0][0] = d00; d[0][1] = d01; d[0][2] = d02;
	d[1][0] = d10; d[1][1] = d11; d[1][2] = d12;
	d[2][0] = d20; d[2][1] = d21; d[2][2] = d22;

	return (*this);
}

inline mat3d& mat3d::operator *= (double a)
{
	d[0][0]*=a; d[0][1]*=a; d[0][2]*=a;
	d[1][0]*=a; d[1][1]*=a; d[1][2]*=a;
	d[2][0]*=a; d[2][1]*=a; d[2][2]*=a;
	return (*this);
}

inline mat3d& mat3d::operator /= (double a)
{
	a = 1./a;
	d[0][0]*=a; d[0][1]*=a; d[0][2]*=a;
	d[1][0]*=a; d[1][1]*=a; d[1][2]*=a;
	d[2][0]*=a; d[2][1]*=a; d[2][2]*=a;
	return (*this);
}

// arithmetic assignment operators for mat3dd
inline mat3d& mat3d::operator += (const mat3dd& m)
{
	d[0][0] += m.d[0];
	d[1][1] += m.d[1];
	d[2][2] += m.d[2];
	return (*this);
}

inline mat3d& mat3d::operator -= (const mat3dd& m)
{
	d[0][0] -= m.d[0];
	d[1][1] -= m.d[1];
	d[2][2] -= m.d[2];
	return (*this);
}

inline mat3d& mat3d::operator *= (const mat3dd& m)
{
	d[0][0] *= m.d[0]; d[0][1] *= m.d[1]; d[0][2] *= m.d[2];
	d[1][0] *= m.d[0]; d[1][1] *= m.d[1]; d[1][2] *= m.d[2];
	d[2][0] *= m.d[0]; d[2][1] *= m.d[1]; d[2][2] *= m.d[2];
	return (*this);
}

// arithmetic operators for mat3ds
inline mat3d& mat3d::operator += (const mat3ds& m)
{
	d[0][0] += m.m[m.XX]; d[0][1] += m.m[m.XY]; d[0][2] += m.m[m.XZ];
	d[1][0] += m.m[m.XY]; d[1][1] += m.m[m.YY]; d[1][2] += m.m[m.YZ];
	d[2][0] += m.m[m.XZ]; d[2][1] += m.m[m.YZ]; d[2][2] += m.m[m.ZZ];
	return (*this);
}

inline mat3d& mat3d::operator -= (const mat3ds& m)
{
	d[0][0] -= m.m[m.XX]; d[0][1] -= m.m[m.XY]; d[0][2] -= m.m[m.XZ];
	d[1][0] -= m.m[m.XY]; d[1][1] -= m.m[m.YY]; d[1][2] -= m.m[m.YZ];
	d[2][0] -= m.m[m.XZ]; d[2][1] -= m.m[m.YZ]; d[2][2] -= m.m[m.ZZ];
	return (*this);
}

inline mat3d& mat3d::operator *= (const mat3ds& m)
{
	double d00 = d[0][0]*m.m[m.XX]+d[0][1]*m.m[m.XY]+d[0][2]*m.m[m.XZ];
	double d01 = d[0][0]*m.m[m.XY]+d[0][1]*m.m[m.YY]+d[0][2]*m.m[m.YZ];
	double d02 = d[0][0]*m.m[m.XZ]+d[0][1]*m.m[m.YZ]+d[0][2]*m.m[m.ZZ];
	double d10 = d[1][0]*m.m[m.XX]+d[1][1]*m.m[m.XY]+d[1][2]*m.m[m.XZ];
	double d11 = d[1][0]*m.m[m.XY]+d[1][1]*m.m[m.YY]+d[1][2]*m.m[m.YZ];
	double d12 = d[1][0]*m.m[m.XZ]+d[1][1]*m.m[m.YZ]+d[1][2]*m.m[m.ZZ];
	double d20 = d[2][0]*m.m[m.XX]+d[2][1]*m.m[m.XY]+d[2][2]*m.m[m.XZ];
	double d21 = d[2][0]*m.m[m.XY]+d[2][1]*m.m[m.YY]+d[2][2]*m.m[m.YZ];
	double d22 = d[2][0]*m.m[m.XZ]+d[2][1]*m.m[m.YZ]+d[2][2]*m.m[m.ZZ];

	d[0][0] = d00; d[0][1] = d01; d[0][2] = d02;
	d[1][0] = d10; d[1][1] = d11; d[1][2] = d12;
	d[2][0] = d20; d[2][1] = d21; d[2][2] = d22;

	return (*this);
}


// matrix-vector multiplication
inline vec3d mat3d::operator * (const vec3d& r) const
{
	return vec3d(d[0][0]*r.x+d[0][1]*r.y+d[0][2]*r.z,
				 d[1][0]*r.x+d[1][1]*r.y+d[1][2]*r.z,
				 d[2][0]*r.x+d[2][1]*r.y+d[2][2]*r.z);
}

// determinant
inline double mat3d::det() const
{
	return (d[0][0]*(d[1][1]*d[2][2] - d[1][2]*d[2][1])
		  + d[0][1]*(d[1][2]*d[2][0] - d[2][2]*d[1][0])
		  + d[0][2]*(d[1][0]*d[2][1] - d[1][1]*d[2][0]));
}

// trace
inline double mat3d::trace() const { return d[0][0]+d[1][1]+d[2][2]; }

inline void mat3d::unit()
{
	d[0][0] = d[1][1] = d[2][2] = 1;
	d[0][1] = d[1][0] = 0;
	d[0][2] = d[2][0] = 0;
	d[1][2] = d[2][1] = 0;
}

// zero the matrix
inline void mat3d::zero()
{
	d[0][0] = d[0][1] = d[0][2] = 0;
	d[1][0] = d[1][1] = d[1][2] = 0;
	d[2][0] = d[2][1] = d[2][2] = 0;
}

// return a column vector from the matrix
inline vec3d mat3d::col(int j) const
{
	return vec3d(d[0][j], d[1][j], d[2][j]);
}

// return a row vector from the matrix
inline vec3d mat3d::row(int j) const
{
	return vec3d(d[j][0], d[j][1], d[j][2]);
}

// set the column of the matrix
inline void mat3d::setCol(int i, const vec3d& a)
{
	d[0][i] = a.x;
	d[1][i] = a.y;
	d[2][i] = a.z;
}

// set the row of the matrix
inline void mat3d::setRow(int i, const vec3d& a)
{
	d[i][0] = a.x;
	d[i][1] = a.y;
	d[i][2] = a.z;
}

// return the symmetric matrix 0.5*(A+A^T)
inline mat3ds mat3d::sym() const
{
	return mat3ds(
		d[0][0],
		d[1][1],
		d[2][2],
		0.5*(d[0][1]+d[1][0]),
		0.5*(d[1][2]+d[2][1]),
		0.5*(d[0][2]+d[2][0]));
}

// return the anti-symmetric matrix 0.5*(A - A^T)
inline mat3da mat3d::skew() const
{
	return mat3da(
		0.5*(d[0][1] - d[1][0]),
		0.5*(d[1][2] - d[2][1]),
		0.5*(d[0][2] - d[2][0]));
}

// return the inverse matrix
inline mat3d mat3d::inverse() const
{
	double D = det();
	assert(D != 0);
	D = 1/D;

	return mat3d(D*(d[1][1]*d[2][2] - d[1][2]*d[2][1]),
				 D*(d[0][2]*d[2][1] - d[0][1]*d[2][2]),
				 D*(d[0][1]*d[1][2] - d[1][1]*d[0][2]),
				 D*(d[1][2]*d[2][0] - d[1][0]*d[2][2]),
				 D*(d[0][0]*d[2][2] - d[0][2]*d[2][0]),
				 D*(d[0][2]*d[1][0] - d[0][0]*d[1][2]),
				 D*(d[1][0]*d[2][1] - d[1][1]*d[2][0]),
				 D*(d[0][1]*d[2][0] - d[0][0]*d[2][1]),
				 D*(d[0][0]*d[1][1] - d[0][1]*d[1][0]));
}

// return the inverse matrix
inline double mat3d::invert(mat3d& Ai)
{
    double D = det();
    if (D != 0) {
        double Di = 1/D;
        Ai = mat3d(Di*(d[1][1]*d[2][2] - d[1][2]*d[2][1]),
                   Di*(d[0][2]*d[2][1] - d[0][1]*d[2][2]),
                   Di*(d[0][1]*d[1][2] - d[1][1]*d[0][2]),
                   Di*(d[1][2]*d[2][0] - d[1][0]*d[2][2]),
                   Di*(d[0][0]*d[2][2] - d[0][2]*d[2][0]),
                   Di*(d[0][2]*d[1][0] - d[0][0]*d[1][2]),
                   Di*(d[1][0]*d[2][1] - d[1][1]*d[2][0]),
                   Di*(d[0][1]*d[2][0] - d[0][0]*d[2][1]),
                   Di*(d[0][0]*d[1][1] - d[0][1]*d[1][0]));
    }
    return D;
}

// return the transpose matrix
inline mat3d mat3d::transpose() const
{
	return mat3d(d[0][0], d[1][0], d[2][0],
				 d[0][1], d[1][1], d[2][1],
				 d[0][2], d[1][2], d[2][2]);
}

// return the transposed inverse matrix
inline mat3d mat3d::transinv() const
{
	double D = det();
	assert(D != 0);
	D = 1/D;

	return mat3d(D*(d[1][1]*d[2][2] - d[1][2]*d[2][1]), // xx
				 D*(d[1][2]*d[2][0] - d[1][0]*d[2][2]), // yx
				 D*(d[1][0]*d[2][1] - d[1][1]*d[2][0]), // zx
				 D*(d[0][2]*d[2][1] - d[0][1]*d[2][2]), // xy
				 D*(d[0][0]*d[2][2] - d[0][2]*d[2][0]), // yy
				 D*(d[0][1]*d[2][0] - d[0][0]*d[2][1]), // zy
				 D*(d[0][1]*d[1][2] - d[1][1]*d[0][2]), // xz
				 D*(d[0][2]*d[1][0] - d[0][0]*d[1][2]), // yz
				 D*(d[0][0]*d[1][1] - d[0][1]*d[1][0])); // zz
}

// calculate the skew symmetric matrix from a vector
inline void mat3d::skew(const vec3d& v)
{
	d[0][0] =    0; d[0][1] = -v.z; d[0][2] =  v.y;
	d[1][0] =  v.z; d[1][1] =    0; d[1][2] = -v.x;
	d[2][0] = -v.y; d[2][1] =  v.x; d[2][2] =    0;
}

// calculate the one-norm (max of absolute column-sum)
inline double mat3d::norm() const
{
	double s, sc;
	sc = fabs(d[0][0]) + fabs(d[1][0]) + fabs(d[2][0]); s = sc;
	sc = fabs(d[0][1]) + fabs(d[1][1]) + fabs(d[2][1]); if (sc > s) s = sc;
	sc = fabs(d[0][2]) + fabs(d[1][2]) + fabs(d[2][2]); if (sc > s) s = sc;
	return s;
}

// double contraction
inline double mat3d::dotdot(const mat3d& T) const
{
	return (T.d[0][0]*d[0][0] + T.d[0][1]*d[0][1] + T.d[0][2]*d[0][2] + T.d[1][0]*d[1][0] + T.d[1][1]*d[1][1] + T.d[1][2]*d[1][2] + T.d[2][0]*d[2][0] + T.d[2][1]*d[2][1] + T.d[2][2]*d[2][2]);
}
