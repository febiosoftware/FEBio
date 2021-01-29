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
#include "stdafx.h"
#include "FETorsionalSpring.h"

BEGIN_FECORE_CLASS(FETorsionalSpring, FEDiscreteElasticMaterial)
	ADD_PARAMETER(m_r, "r");
END_FECORE_CLASS();

FETorsionalSpring::FETorsionalSpring(FEModel* fem) : FEDiscreteElasticMaterial(fem)
{
	m_r = 0.0;
}

// evaluate the force at a discrete element
vec3d FETorsionalSpring::Force(FEDiscreteMaterialPoint& mp)
{
	vec3d e0 = mp.m_dr0; e0.unit();
	vec3d et = mp.m_drt; 
	double l = et.unit();

	double w = e0 * et;

	mat3ds exe = dyad(et);
	mat3dd I(1.0);

	mat3ds P = (I - exe) / l;
	vec3d t = P * e0;

	double F = -m_r*w;

	return t*F;
}

// evaluate the stiffness at a discrete element (= dF / dr)
mat3d FETorsionalSpring::Stiffness(FEDiscreteMaterialPoint& mp)
{
	vec3d e0 = mp.m_dr0; e0.unit();
	vec3d et = mp.m_drt;
	double l = et.unit();

	double w = e0 * et;

	mat3ds exe = dyad(et);
	mat3dd I(1.0);

	mat3ds P = (I - exe) / l;
	vec3d t = P * e0;

	mat3ds Q = dyads(e0, et) + I*w - exe*(3.0*w);

	mat3ds T = dyad(t);

	return T*(m_r) + Q*(m_r*w);
}

double FETorsionalSpring::StrainEnergy(FEDiscreteMaterialPoint& mp)
{
	vec3d e0 = mp.m_dr0; e0.unit();
	vec3d et = mp.m_drt; 
	double l = et.unit();

	double w = e0 * et;

	double E = m_r*(1-w*w)/2;

	return E;

}