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

namespace fad {

	struct number {
		double real, partial;
		number() : real(0.0), partial(0.0) {}
		number(double real, double partial = 0) : real(real), partial(partial) {}
		void operator = (const double& a) { real = a; }
	};

	number operator + (const number& a, const number& b)
	{
		return number(a.real + b.real, a.partial + b.partial);
	}

	number operator - (const number& a, const number& b)
	{
		return number(a.real - b.real, a.partial - b.partial);
	}

	number operator * (const number& a, const number& b)
	{
		return number(a.real * b.real, a.partial * b.real + a.real * b.partial);
	}

	number operator / (const number& a, const number& b)
	{
		return number(a.real / b.real, a.partial / b.real - a.real * b.partial / (b.real * b.real));
	}

	number Ln(const number& a)
	{
		return number(log(a.real), a.partial / a.real);
	}

	number Sqrt(const number& a)
	{
		double s = sqrt(a.real);
		return number(s, 0.5 * a.partial / s);
	}

	struct Mat3ds
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
		Mat3ds() {}
		Mat3ds(const mat3ds& C)
		{
			m[0] = C.xx();
			m[1] = C.xy();
			m[2] = C.yy();
			m[3] = C.xz();
			m[4] = C.yz();
			m[5] = C.zz();
		}

		number Trace() const { return m[XX] + m[YY] + m[ZZ]; }

		number Det() const {
			return (m[XX] * (m[YY] * m[ZZ] - m[YZ] * m[YZ])
				+ m[XY] * (m[YZ] * m[XZ] - m[ZZ] * m[XY])
				+ m[XZ] * (m[XY] * m[YZ] - m[YY] * m[XZ]));
		}
	};

	typedef number (*ScalarFunctionOfMat3ds)(const Mat3ds& C);
	typedef Mat3ds (*Mat3dsFunctionOfMat3ds)(const Mat3ds& C);

	FECORE_API double Evaluate(ScalarFunctionOfMat3ds, const mat3ds& C);

	FECORE_API mat3ds Derive(ScalarFunctionOfMat3ds, const mat3ds& C);
}
