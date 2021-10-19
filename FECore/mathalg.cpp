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
#include "stdafx.h"
#include "mathalg.h"

mat3ds Log(const mat3ds& p, const mat3ds& X)
{
	double l[3], s[3];
	vec3d u[3], v[3];

	// evaluate eigen-decomposition of p
	p.eigen(l, u);
	mat3d U(u[0], u[1], u[2]);
	mat3dd L(l[0], l[1], l[2]);
	mat3dd rootL(sqrt(l[0]), sqrt(l[1]), sqrt(l[2]));

	mat3d G = U * rootL;
	mat3d Gi = G.inverse();

	mat3ds Y = (Gi * X*Gi.transpose()).sym();

	Y.eigen(s, v);
	mat3d V = mat3d(v[0], v[1], v[2]);

	mat3d GV = G * V;

	mat3dd logS(log(s[0]), log(s[1]), log(s[2]));

	mat3d LogX = (GV)*logS*(GV.transpose());

	return LogX.sym();
}

mat3ds Exp(const mat3ds& p, const mat3ds& X)
{
	double l[3], s[3];
	vec3d u[3], v[3];

	// evaluate eigen-decomposition of p
	p.eigen(l, u);
	mat3d U(u[0], u[1], u[2]);
	mat3dd L(l[0], l[1], l[2]);
	mat3dd rootL(sqrt(l[0]), sqrt(l[1]), sqrt(l[2]));

	mat3d G = U * rootL;
	mat3d Gi = G.inverse();

	mat3ds Y = (Gi * X*Gi.transpose()).sym();

	Y.eigen(s, v);
	mat3d V = mat3d(v[0], v[1], v[2]);

	mat3d GV = G * V;

	mat3dd expS(exp(s[0]), exp(s[1]), exp(s[2]));

	mat3d ExpX = (GV)*expS*(GV.transpose());

	return ExpX.sym();
}

mat3ds weightedAverageStructureTensor(mat3ds* d, double* w, int n)
{
	const double eps = 1.0e-9;

	mat3ds mu(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
	double tau = 1.0;

	double normXi = 0.0, normXip = 0.0;
	mat3ds Xi, Xip;
	int nc = 0;
	do
	{
		Xip = Xi;
		normXip = normXi;

		Xi = weightedAverage<mat3ds>(d, w, n, [&](const mat3ds& a) {
			return Log(mu, a); 
		});

		mu = Exp(mu*tau, Xi);

		normXi = Xi.norm();

		if ((nc != 0) && (normXi > normXip))
		{
			Xi = Xip;
			tau *= 0.5;
		}
		nc++;
	} 
	while (normXi > eps);

	return mu;
}
