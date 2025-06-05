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
#include "FEAxialBodyForce.h"
#include "FEElasticMaterial.h"

BEGIN_FECORE_CLASS(FEAxialBodyForce, FEBodyForce);
	ADD_PARAMETER(s, "scale");
	ADD_PARAMETER(p, "p");
	ADD_PARAMETER(c, "c");
	ADD_PARAMETER(n, "axis");
END_FECORE_CLASS();

FEAxialBodyForce::FEAxialBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
	p = 0;
	c = vec3d(0, 0, 0);
	n = vec3d(0, 0, 1);
	s = 1;
}

vec3d FEAxialBodyForce::force(FEMaterialPoint& mp)
{
	vec3d x = mp.m_rt;
	vec3d q = x - c;
	vec3d r = q - n*((q*n)/(n*n));

	if (p == 0)
	{
		vec3d f = r * s;
		return f;
	}
	else
	{
		double R = r.norm();
		double Rp = pow(R, p);
		vec3d f = r *(s*Rp);
		return f;
	}
}

mat3d FEAxialBodyForce::stiffness(FEMaterialPoint& mp)
{ 
	mat3d K; K.zero();
	if (p == 0)
	{
		K[0][0] = -s;
		K[1][1] = -s;
		K[2][2] = -s;
	}
	else
	{
		vec3d x = mp.m_rt;
		vec3d q = x - c;
		vec3d r = q - n * ((q * n) / (n * n));

		double R = x.norm();
		double Rp = pow(R, p);
		double Rpm2 = pow(R, p-2 );

		K[0][0] = -s*(p * Rpm2 * r.x*r.x + Rp);
		K[0][1] = -s*(p * Rpm2 * r.x*r.y);
		K[0][2] = -s*(p * Rpm2 * r.x*r.z);

		K[1][0] = -s*(p * Rpm2 * r.y*r.x);
		K[1][1] = -s*(p * Rpm2 * r.y*r.y + Rp);
		K[1][2] = -s*(p * Rpm2 * r.y*r.z);

		K[2][0] = -s*(p * Rpm2 * r.z*r.x);
		K[2][1] = -s*(p * Rpm2 * r.z*r.y);
		K[2][2] = -s*(p * Rpm2 * r.z*r.z + Rp);
	}
	return K; 
}

double FEAxialBodyForce::divforce(FEMaterialPoint& mp)
{
	// TODO: implement this
	assert(false);
	return 0;
}
