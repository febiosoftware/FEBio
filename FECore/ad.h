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
#include "fecore_api.h"
#include <functional>

namespace ad {

	struct number {
		double r, dr;
		number() : r(0.0), dr(0.0) {}
		number(double v, double dv = 0) : r(v), dr(dv) {}
		void operator = (const double& a) { r = a; dr = 0.0; }
	};

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
		return number(a.r * b.r , a.dr * b.r + a.r * b.dr);
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
		return number(a.r / b.r, (a.dr - a.r* b.dr/ b.r) / b.r);
	}

	inline number operator / (double a, const number& b)
	{
		return number(a / b.r,  - a * b.dr / (b.r* b.r));
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
			A.xx()*b,
			A.yy()*b,
			A.zz()*b,
			A.xy()*b,
			A.yz()*b,
			A.xz()*b
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
	FECORE_API ::tens4ds Derive(std::function<mat3ds(mat3ds& C)> S, const ::mat3ds& C);
}
