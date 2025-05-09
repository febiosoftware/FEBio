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
#include "FERadialBodyForce.h"
#include "FEElasticMaterial.h"

BEGIN_FECORE_CLASS(FERadialBodyForce, FEBodyForce);
	ADD_PARAMETER(p, "p");
	ADD_PARAMETER(c, "c");
	ADD_PARAMETER(s[0], "sx");
	ADD_PARAMETER(s[1], "sy");
	ADD_PARAMETER(s[2], "sz");
END_FECORE_CLASS();

FERadialBodyForce::FERadialBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
	p = 0;
	c = vec3d(0, 0, 0);
	s[0] = 1;
	s[1] = 1;
	s[2] = 1;
}

vec3d FERadialBodyForce::force(FEMaterialPoint& mp)
{
	vec3d& x = mp.m_rt - c;
	if (p == 0)
	{
		vec3d f(s[0] * x.x, s[1] * x.y, s[2] * x.z);
		return f;
	}
	else
	{
		double R = x.norm();
		double Rp = pow(R, p);
		vec3d f(s[0] * Rp * x.x, s[1] * Rp * x.y, s[2] * Rp * x.z);
		return f;
	}
}

mat3d FERadialBodyForce::stiffness(FEMaterialPoint& mp)
{ 
	mat3d K; K.zero();
	if (p == 0)
	{
		K[0][0] = -s[0];
		K[1][1] = -s[1];
		K[2][2] = -s[2];
	}
	else
	{
		vec3d& x = mp.m_rt - c;
		mat3d XX = (x & x);

		double R = x.norm();
		double Rp = -pow(R, p);
		double Rpm2 = -pow(R, p-2 );

		K[0][0] = s[0] * p * Rpm2 * XX[0][0] + s[0] * Rp;
		K[0][1] = s[0] * p * Rpm2 * XX[0][1];
		K[0][2] = s[0] * p * Rpm2 * XX[0][2];

		K[1][0] = s[1] * p * Rpm2 * XX[1][0];
		K[1][1] = s[1] * p * Rpm2 * XX[1][1] + s[1] * Rp;
		K[1][2] = s[1] * p * Rpm2 * XX[1][2];

		K[2][0] = s[2] * p * Rpm2 * XX[2][0];
		K[2][1] = s[2] * p * Rpm2 * XX[2][1];
		K[2][2] = s[2] * p * Rpm2 * XX[2][2] + s[2] * Rp;
	}
	return K; 
}

double FERadialBodyForce::divforce(FEMaterialPoint& mp)
{
	// TODO: implement this
	assert(false);
	return 0;
}
