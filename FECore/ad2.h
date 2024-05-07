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
#include "tens4d.h"
#include "fecore_api.h"
#include <functional>

namespace ad2 {

	struct number {
		double r, d1, d2, dd;
		number() : r(0.0), d1(0.0), d2(0.0), dd(0.0) {}
		number(double v, double v1 = 0, double v2=0, double ddv = 0) : r(v), d1(v1), d2(v2), dd(ddv) {}
		void operator = (const double& a) { r = a; d1 = d2 = dd = 0.0; }
	};

	// addition
	inline number operator + (const number& a, const number& b)
	{
		return number(a.r + b.r, a.d1 + b.d1, a.d2 + b.d2, a.dd + b.dd);
	}

	inline number operator + (double a, const number& b)
	{
		return number(a + b.r, b.d1, b.d2, b.dd);
	}

	inline number operator + (const number& a, double b)
	{
		return number(a.r + b, a.d1, a.d2, a.dd);
	}

	// subtraction
	inline number operator - (const number& a, const number& b)
	{
		return number(a.r - b.r, a.d1 - b.d1, a.d2 - b.d2, a.dd - b.dd);
	}

	inline number operator - (double a, const number& b)
	{
		return number(a - b.r, -b.d1, -b.d2, -b.dd);
	}

	inline number operator - (const number& a, double b)
	{
		return number(a.r - b, a.d1, a.d2, a.dd);
	}

	// multiplication
	inline number operator * (const number& a, const number& b)
	{
		return number(
			a.r * b.r,
			a.d1 * b.r + a.r * b.d1,
			a.d2 * b.r + a.r * b.d2,
			a.dd * b.r + a.d1 * b.d2 + a.d2 * b.d1 + a.r * b.dd
			);
	}

	inline number operator * (double a, const number& b)
	{
		return number(
			a * b.r,
			a * b.d1,
			a * b.d2,
			a * b.dd
		);
	}

	inline number operator * (const number& a, double b)
	{
		return number(
			a.r * b,
			a.d1 * b,
			a.d2 * b,
			a.dd * b
		);
	}

	// invert
	inline number inv(const number& a)
	{
		return number(
			1.0/a.r,
			-a.d1/(a.r*a.r),
			-a.d2/(a.r*a.r),
			((2.0 * a.d1 * a.d2) / a.r - a.dd)/(a.r*a.r)
		);
	}

	// division
	inline number operator / (const number& a, const number& b)
	{
		return a*inv(b);
	}

	inline number operator / (double a, const number& b)
	{
		return a * inv(b);
	}

	inline number operator / (const number& a, double b)
	{
		return number(
			a.r / b,
			a.d1 / b,
			a.d2 / b,
			a.dd / b
		);
	}

	// math functions
	inline number log(const number& a)
	{
		return number(
			::log(a.r), 
			a.d1 / a.r,
			a.d2 / a.r,
			(a.dd - a.d1*a.d2/a.r)/a.r
			);
	}

	inline number sqrt(const number& a)
	{
		double s = ::sqrt(a.r);
		return number(
			s, 
			0.5 * a.d1 / s,
			0.5 * a.d2 / s,
			0.5 * a.dd / s - 0.25*a.d1*a.d2/(s*s*s)
			);
	}

	inline number exp(const number& a)
	{
		double e = ::exp(a.r);
		return number(
			e, 
			e * a.d1,
			e * a.d2,
			e * (a.d1*a.d2 + a.dd)
			);
	}

	inline number pow(const number& a, double e)
	{
		if (e == 0.0) return number(1.0);

		double b = ::pow(a.r, e - 2.0);
		return number(
			b * a.r * a.r, 
			e * a.r * b * a.d1,
			e * a.r * b * a.d2,
			e * b * ((e - 1.0)* a.d1 * a.d2 + a.r * a.dd)
		);
	}

	inline number sin(const number& a)
	{
		double sa = ::sin(a.r);
		double ca = ::cos(a.r);
		return number(
			sa,
			ca * a.d1,
			ca * a.d2,
			ca * a.dd - sa * a.d1 * a.d2
		);
	}

	inline number cos(const number& a)
	{
		double sa = ::sin(a.r);
		double ca = ::cos(a.r);
		return number(
			ca,
			-sa*a.d1,
			-sa*a.d2,
			-sa*a.dd - ca*a.d1*a.d2
			);
	}

	inline number cosh(const number& a)
	{
		double sa = ::sinh(a.r);
		double ca = ::cosh(a.r);
		return number(
			ca,
			sa * a.d1,
			sa * a.d2,
			ca * a.d1 * a.d2 + sa * a.dd
			);
	}

	inline number sinh(const number& a)
	{
		double sa = ::sinh(a.r);
		double ca = ::cosh(a.r);
		return number(
			sa,
			ca * a.d1,
			ca * a.d2,
			sa * a.d1 * a.d2 + ca * a.dd
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
	FECORE_API ::mat3ds Derive(std::function<number(mat3ds& C)> W, const ::mat3ds& C);
	FECORE_API ::tens4ds Derive2(std::function<number(mat3ds& C)> W, const ::mat3ds& C);
}
