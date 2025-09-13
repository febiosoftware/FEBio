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
#include "ad.h"

double ad::Evaluate(std::function<ad::number(vec3d&)> W, const ::vec3d& a)
{
	ad::vec3d da(a);
	return W(da).r;
}

::vec3d ad::Grad(std::function<ad::number(vec3d&)> W, const ::vec3d& a)
{
	ad::vec3d da(a);
	da.x.dr = 1; da.y.dr = 0; da.z.dr = 0; double dx = W(da).dr;
	da.x.dr = 0; da.y.dr = 1; da.z.dr = 0; double dy = W(da).dr;
	da.x.dr = 0; da.y.dr = 0; da.z.dr = 1; double dz = W(da).dr;
	return ::vec3d(dx, dy, dz);
}

::mat3d ad::Grad(std::function<ad::vec3d(ad::vec3d&)> F, const ::vec3d& a)
{
	ad::vec3d da(a);
	double D[3][3] = { 0 };
	for (int i = 0; i < 3; ++i)
	{
		da[i].dr = 1;
		::vec3d Fi = F(da).partials();
		da[i].dr = 0;
		D[0][i] = Fi.x;
		D[1][i] = Fi.y;
		D[2][i] = Fi.z;
	}
	return ::mat3d(D);
}

double ad::Evaluate(std::function<ad::number(ad::mat3ds& C)> W, const ::mat3ds& C)
{
	ad::mat3ds dC(C);
	return W(dC).r;
}

::mat3ds ad::Derive(std::function<ad::number(ad::mat3ds& C)> W, const ::mat3ds& C)
{
	ad::mat3ds dC(C);
	double S[6] = { 0.0 };
	for (int i = 0; i < 6; ++i)
	{
		dC[i].dr = 1;
		S[i] = W(dC).dr;
		dC[i].dr = 0;
	}
	return ::mat3ds(S[0], S[2], S[5], 0.5 * S[1], 0.5 * S[4], 0.5 * S[3]);
}

::tens4ds ad::Derive(std::function<mat3ds(mat3ds& C)> S, const ::mat3ds& C)
{
	ad::mat3ds dC(C);

	constexpr int l[6] = { 0, 2, 5, 1, 4, 3 };
	double D[6][6] = { 0 };
	for (int i = 0; i < 6; ++i)
	{
		int n = l[i];
		dC[n].dr = 1;
		::mat3ds Si = S(dC).partials();
		dC[n].dr = 0;

		D[0][i] = Si.xx();
		D[1][i] = Si.yy();
		D[2][i] = Si.zz();
		D[3][i] = 0.5 * Si.xy();
		D[4][i] = 0.5 * Si.yz();
		D[5][i] = 0.5 * Si.xz();
	}

	return tens4ds(D);
}
