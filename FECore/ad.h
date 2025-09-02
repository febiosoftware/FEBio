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
#include "mat3d.h"
#include "tens4d.h"
#include "quatd.h"
#include "fecore_api.h"
#include <functional>

namespace ad {

	struct number {
		double r, dr;
		number() : r(0.0), dr(0.0) {}
		number(double v, double dv = 0) : r(v), dr(dv) {}
		void operator = (const double& a) { r = a; dr = 0.0; }
	};

	// negation
	inline number operator - (const number& a)
	{
		return number(-a.r, -a.dr);
	}

	// addition
	inline number operator + (const number& a, const number& b)
	{
		return number(a.r + b.r, a.dr + b.dr);
	}

	inline number operator + (double a, const number& b)
	{
		return number(a + b.r, b.dr);
	}

	inline number operator + (const number& a, double b)
	{
		return number(a.r + b, a.dr);
	}

	// subtraction
	inline number operator - (const number& a, const number& b)
	{
		return number(a.r - b.r, a.dr - b.dr);
	}

	inline number operator - (double a, const number& b)
	{
		return number(a - b.r, -b.dr);
	}

	inline number operator - (const number& a, double b)
	{
		return number(a.r - b, a.dr);
	}

	// multiplication
	inline number operator * (const number& a, const number& b)
	{
		return number(a.r * b.r, a.dr * b.r + a.r * b.dr);
	}

	inline number operator * (double a, const number& b)
	{
		return number(a * b.r, a * b.dr);
	}

	inline number operator * (const number& a, double b)
	{
		return number(a.r * b, a.dr * b);
	}

	// division
	inline number operator / (const number& a, const number& b)
	{
		return number(a.r / b.r, (a.dr - a.r * b.dr / b.r) / b.r);
	}

	inline number operator / (double a, const number& b)
	{
		return number(a / b.r, -a * b.dr / (b.r * b.r));
	}

	inline number operator / (const number& a, double b)
	{
		return number(a.r / b, a.dr / b);
	}

	// math functions
	inline number log(const number& a)
	{
		return number(::log(a.r), a.dr / a.r);
	}

	inline number sqrt(const number& a)
	{
		double s = ::sqrt(a.r);
		return number(s, 0.5 * a.dr / s);
	}

	inline number exp(const number& a)
	{
		double e = ::exp(a.r);
		return number(e, e * a.dr);
	}

	inline number pow(const number& a, double e)
	{
		if (e == 0.0) return number(1.0);
		double b = ::pow(a.r, e - 1.0);
		return number(a.r * b, e * b * a.dr);
	}

	inline number sin(const number& a)
	{
		return number(::sin(a.r), ::cos(a.r) * a.dr);
	}

	inline number cos(const number& a)
	{
		return number(::cos(a.r), -::sin(a.r) * a.dr);
	}

	inline number cosh(const number& a)
	{
		return number(::cosh(a.r), ::sinh(a.r) * a.dr);
	}

	inline number sinh(const number& a)
	{
		return number(::sinh(a.r), ::cosh(a.r) * a.dr);
	}

	struct vec3d
	{
		number x, y, z;
		vec3d() {}
		vec3d(const ::vec3d& a) : x(a.x), y(a.y), z(a.z) {}
		vec3d(double a, double b, double c) : x(a), y(b), z(c) {}
		vec3d(const number& a, const number& b, const number& c) : x(a), y(b), z(c) {}

		::vec3d values() const { return ::vec3d(x.r, y.r, z.r); }
		::vec3d partials() const { return ::vec3d(x.dr, y.dr, z.dr); }

		number& operator [] (size_t n) { return (&x)[n]; }

		number length() const { return sqrt(x * x + y * y + z * z); }
	};

	inline number operator * (const vec3d& a, const vec3d& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	inline vec3d operator * (double a, const vec3d& b)
	{
		return vec3d(a * b.x, a * b.y, a * b.z);
	}

	inline vec3d operator * (const number& a, const vec3d& b)
	{
		return vec3d(a * b.x, a * b.y, a * b.z);
	}

	inline vec3d operator * (const vec3d& a, double b)
	{
		return vec3d(a.x * b, a.y * b, a.z * b);
	}

	inline vec3d operator * (const vec3d& a, const number& b)
	{
		return vec3d(a.x * b, a.y * b, a.z * b);
	}

	inline vec3d operator + (const vec3d& a, const vec3d& b)
	{
		return vec3d(a.x + b.x, a.y + b.y, a.z + b.z);
	}

	inline vec3d operator - (const vec3d& a, const vec3d& b)
	{
		return vec3d(a.x - b.x, a.y - b.y, a.z - b.z);
	}

	inline vec3d operator - (const vec3d& a)
	{
		return vec3d(-a.x, -a.y, -a.z);
	}

	inline vec3d cross(const vec3d& a, const vec3d& b)
	{
		return vec3d(a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x);
	}

	// alternative cross product operator. 
	// Make sure to include parentheses when using this operator since ^ has a lower precedence than some other operators.
	// e.g. use (a ^ b) * c instead of a ^ b * c
	inline vec3d operator ^ (const vec3d& a, const vec3d& b)
	{
		return cross(a, b);
	}

	FECORE_API double Evaluate(std::function<number(ad::vec3d&)> W, const ::vec3d& a);
	FECORE_API::vec3d Grad(std::function<number(ad::vec3d&)> W, const ::vec3d& a);
	FECORE_API::mat3d Grad(std::function<ad::vec3d(ad::vec3d&)> F, const ::vec3d& a);

	struct mat3d
	{
		number m[3][3];
		number* operator [] (size_t n) { return m[n]; }
		const number* operator [] (size_t n) const { return m[n]; }

		mat3d() {}
		mat3d(const ::mat3d& A)
		{
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j) m[i][j] = A[i][j];
		}
		mat3d(const number& a11, const number& a12, const number& a13,
			const number& a21, const number& a22, const number& a23,
			const number& a31, const number& a32, const number& a33)
		{
			m[0][0] = a11; m[0][1] = a12; m[0][2] = a13;
			m[1][0] = a21; m[1][1] = a22; m[1][2] = a23;
			m[2][0] = a31; m[2][1] = a32; m[2][2] = a33;
		}

		mat3d(double d)
		{
			m[0][0].r = d; m[0][1].r = 0.0; m[0][2].r = 0.0;
			m[1][0].r = 0.0; m[1][1].r = d; m[1][2].r = 0.0;
			m[2][0].r = 0.0; m[2][1].r = 0.0; m[2][2].r = d;
		}
		::mat3d values() const
		{
			return ::mat3d(m[0][0].r, m[0][1].r, m[0][2].r,
				m[1][0].r, m[1][1].r, m[1][2].r,
				m[2][0].r, m[2][1].r, m[2][2].r);
		}
		::mat3d partials() const
		{
			return ::mat3d(m[0][0].dr, m[0][1].dr, m[0][2].dr,
				m[1][0].dr, m[1][1].dr, m[1][2].dr,
				m[2][0].dr, m[2][1].dr, m[2][2].dr);
		}
	};

	inline mat3d operator * (const mat3d& A, const mat3d& B)
	{
		mat3d C(0.0);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					C[i][j] = C[i][j] + A[i][k] * B[k][j];
		return C;
	}

	inline vec3d operator * (const mat3d& A, const vec3d& v)
	{
		return vec3d(
			A[0][0] * v.x + A[0][1] * v.y + A[0][2] * v.z,
			A[1][0] * v.x + A[1][1] * v.y + A[1][2] * v.z,
			A[2][0] * v.x + A[2][1] * v.y + A[2][2] * v.z
		);
	}

	inline mat3d dyad(const vec3d& a)
	{
		return mat3d(
			a.x * a.x, a.x * a.y, a.x * a.z,
			a.y * a.x, a.y * a.y, a.y * a.z,
			a.z * a.x, a.z * a.y, a.z * a.z
		);
	}

	inline mat3d dyad(const vec3d& a, const vec3d& b)
	{
		return mat3d(
			a.x * b.x, a.x * b.y, a.x * b.z,
			a.y * b.x, a.y * b.y, a.y * b.z,
			a.z * b.x, a.z * b.y, a.z * b.z
		);
	}

	struct mat3ds
	{
		// This enumeration can be used to remember the order
		// in which the components are stored.
		enum {
			XX = 0,
			XY = 1,
			YY = 2,
			XZ = 3,
			YZ = 4,
			ZZ = 5
		};

		number m[6]; // {xx,xy,yy,xz,yz,zz}
		number& operator [] (size_t n) { return m[n]; }
		mat3ds() {}
		mat3ds(
			const number& xx,
			const number& yy,
			const number& zz,
			const number& xy,
			const number& yz,
			const number& xz)
		{
			m[XX] = xx;
			m[YY] = yy;
			m[ZZ] = zz;
			m[XY] = xy;
			m[YZ] = yz;
			m[XZ] = xz;
		}

		mat3ds(const ::mat3ds& C)
		{
			m[0] = C.xx();
			m[1] = C.xy();
			m[2] = C.yy();
			m[3] = C.xz();
			m[4] = C.yz();
			m[5] = C.zz();
		}

		mat3ds(double d)
		{
			m[XX].r = d;
			m[YY].r = d;
			m[ZZ].r = d;
		}

		number& xx() { return m[XX]; }
		number& yy() { return m[YY]; }
		number& zz() { return m[ZZ]; }
		number& xy() { return m[XY]; }
		number& yz() { return m[YZ]; }
		number& xz() { return m[XZ]; }

		const number& xx() const { return m[XX]; }
		const number& yy() const { return m[YY]; }
		const number& zz() const { return m[ZZ]; }
		const number& xy() const { return m[XY]; }
		const number& yz() const { return m[YZ]; }
		const number& xz() const { return m[XZ]; }

		::mat3ds values() const
		{
			return ::mat3ds(m[XX].r, m[YY].r, m[ZZ].r, m[XY].r, m[YZ].r, m[XZ].r);
		}

		::mat3ds partials() const
		{
			return ::mat3ds(m[XX].dr, m[YY].dr, m[ZZ].dr, m[XY].dr, m[YZ].dr, m[XZ].dr);
		}

		// functions
		number tr() const { return m[XX] + m[YY] + m[ZZ]; }

		number det() const {
			return (m[XX] * (m[YY] * m[ZZ] - m[YZ] * m[YZ])
				+ m[XY] * (m[YZ] * m[XZ] - m[ZZ] * m[XY])
				+ m[XZ] * (m[XY] * m[YZ] - m[YY] * m[XZ]));
		}

		// double contraction
		number dotdot(const mat3ds& B) const
		{
			const number* n = B.m;
			return m[XX] * n[XX] + m[YY] * n[YY] + m[ZZ] * n[ZZ] + 2.0 * (m[XY] * n[XY] + m[YZ] * n[YZ] + m[XZ] * n[XZ]);
		}

		mat3ds inverse() const
		{
			number Di = 1.0 / det();

			return mat3ds(
				Di * (m[YY] * m[ZZ] - m[YZ] * m[YZ]),
				Di * (m[XX] * m[ZZ] - m[XZ] * m[XZ]),
				Di * (m[XX] * m[YY] - m[XY] * m[XY]),
				Di * (m[XZ] * m[YZ] - m[XY] * m[ZZ]),
				Di * (m[XY] * m[XZ] - m[XX] * m[YZ]),
				Di * (m[XY] * m[YZ] - m[YY] * m[XZ]));
		}

		// return the square 
		mat3ds sqr() const
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
	};

	// arithmetic operations
	inline mat3ds operator + (const mat3ds& A, const mat3ds& B)
	{
		return mat3ds(
			A.xx() + B.xx(),
			A.yy() + B.yy(),
			A.zz() + B.zz(),
			A.xy() + B.xy(),
			A.yz() + B.yz(),
			A.xz() + B.xz()
		);
	}

	inline mat3ds operator - (const mat3ds& A, const mat3ds& B)
	{
		return mat3ds(
			A.xx() - B.xx(),
			A.yy() - B.yy(),
			A.zz() - B.zz(),
			A.xy() - B.xy(),
			A.yz() - B.yz(),
			A.xz() - B.xz()
		);
	}

	inline mat3ds operator * (const mat3ds& A, double b)
	{
		return mat3ds(
			A.xx() * b,
			A.yy() * b,
			A.zz() * b,
			A.xy() * b,
			A.yz() * b,
			A.xz() * b
		);
	}

	inline mat3ds operator * (double a, const mat3ds& B)
	{
		return mat3ds(
			B.xx() * a,
			B.yy() * a,
			B.zz() * a,
			B.xy() * a,
			B.yz() * a,
			B.xz() * a
		);
	}

	inline mat3ds operator * (const mat3ds& A, const number& b)
	{
		return mat3ds(
			A.xx() * b,
			A.yy() * b,
			A.zz() * b,
			A.xy() * b,
			A.yz() * b,
			A.xz() * b
		);
	}

	FECORE_API double Evaluate(std::function<number(mat3ds& C)> W, const ::mat3ds& C);
	FECORE_API::mat3ds Derive(std::function<number(mat3ds& C)> W, const ::mat3ds& C);
	FECORE_API::tens4ds Derive(std::function<mat3ds(mat3ds& C)> S, const ::mat3ds& C);
}

// This namespace contains classes and functions that can be used to evaluate directional 
// derivatives for vec3d and quatd variables automatically. 
// To use this, first create some variables, either dd::vec3d or dd::quat3d.
// For example,
// 
//  dd::vec3d r(1,0,0);
//  dd::quatd q(1,0,0);
//
// then, create functions using these variables. For example,
//
//  auto F = [&]() { 
//    ::vec3d z0(0, 1, 0);
//    return r + q*z0;
// };
//
// Make sure to capture by reference!
// Finally, to evaluate the directional derivatives: 
//
//  mat3d dFr = dd::D(F, r);
//  mat3d dFq = dd::D(F, q);
//
namespace dd {

	// this class defines a number that can be used in calculations.
	// Never initialize directly. Instead, use vec3d operators to construct numbers.
	// For example.
	// dd::vec3d r(1,0,0);
	// dd::number l = dd::sqrt(r*r);
	struct number {
		double v = 0;
		::vec3d dv = ::vec3d(0.0, 0.0, 0.0);
	};

	inline dd::number operator - (double a, const dd::number& n)
	{
		return { a - n.v, -n.dv };
	}

	inline dd::number operator * (double a, const dd::number& n)
	{
		return { a * n.v, n.dv * a };
	}

	inline dd::number operator * (const dd::number& n, double a)
	{
		return { a * n.v, n.dv * a };
	}

	inline dd::number operator / (double a, const dd::number& n)
	{
		return { a / n.v, n.dv * (-a / (n.v * n.v)) };
	}

	inline dd::number sqrt(const dd::number& a)
	{
		return { ::sqrt(a.v), a.dv * (0.5 / ::sqrt(a.v)) };
	};

	// represents vector variables.
	struct vec3d {
		::vec3d v = ::vec3d(0, 0, 0);
		::mat3d dv = ::mat3d(0.0);

		void activate(bool b) { dv = (b ? ::mat3d::identity() : ::mat3d(0.0)); }

		vec3d(::vec3d a) : v(a), dv(::mat3d(0.0)) {}
		vec3d(::vec3d a, ::mat3d da) : v(a), dv(da) {}
		vec3d(double x, double y, double z) : v(::vec3d(x, y, z)), dv(::mat3d(0.0)) {}
	};

	inline dd::number operator * (const dd::vec3d& a, const dd::vec3d& b)
	{
		return dd::number{ a.v * b.v, b.dv.transpose() * a.v + a.dv.transpose() * b.v };
	}

	inline dd::vec3d operator * (const dd::vec3d& r, const dd::number& a)
	{
		return { r.v * a.v, (r.v & a.dv) + r.dv*a.v };
	}

	inline dd::vec3d operator - (const dd::vec3d& a)
	{
		return { -a.v, -a.dv };
	}

	inline dd::vec3d operator + (const dd::vec3d& a, const dd::vec3d& b)
	{
		return { a.v + b.v, a.dv + b.dv };
	}

	inline dd::vec3d operator - (const dd::vec3d& a, const dd::vec3d& b)
	{
		return { a.v - b.v, a.dv - b.dv };
	}

	inline dd::vec3d operator ^ (const dd::vec3d& a, const dd::vec3d& b)
	{
		::mat3da ahat(a.v);
		::mat3da bhat(b.v);
		return { (a.v ^ b.v), ahat * b.dv - bhat * a.dv };
	}

	inline dd::vec3d operator * (const dd::vec3d& a, double s)
	{
		return { a.v * s, a.dv * s };
	}

	// represents rotational variables
	struct quatd {
		::quatd q;
		::mat3d dq;

		quatd(::quatd r) { q = r; dq = ::mat3d(0.0); }
		quatd(::vec3d r) { q = ::quatd(r); dq = ::mat3d(0.0); }
		quatd(double x, double y, double z) { q = ::quatd(::vec3d(x,y,z)); dq = ::mat3d(0.0); }
		quatd(double x, double y, double z, double w) { q = ::quatd(x,y,z,w); dq = ::mat3d(0.0); }

		void activate(bool b) { dq = (b ? ::mat3d::identity() : ::mat3d(0.0)); }

		dd::vec3d operator * (const ::vec3d& a)
		{
			::vec3d qa = q * a;
			return { qa, ::mat3da(-qa) * dq };
		}
	};

	// Use this function to return the directional derivative of "F" w.r.t. "r".
	template <typename T>
	::mat3d D(std::function<dd::vec3d()> F, T& r)
	{
		r.activate(true);
		::mat3d dF = F().dv;
		r.activate(false);
		return dF;
	}
}
