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
#include "FEMembraneMaterial.h"

//-----------------------------------------------------------------------------
// convet the displacement gradient to a strain
void FEMembraneMaterialPoint::strain(double e[3])
{
	const double h1[6] = {1, 0, 0, 0, 0, 0};
	const double h2[6] = {0, 0, 0, 0, 1, 0};
	const double h3[6] = {0, 1, 0, 1, 0, 0};

	const double H1[6][6] = {
		{1, 0, 0, 0, 0, 0},
		{0, 1, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0}};

	const double H2[6][6] = {
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 0, 0},
		{0, 0, 0, 0, 1, 0},
		{0, 0, 0, 0, 0, 1}};

	const double H3[6][6] = {
		{0, 0, 0, 1, 0, 0},
		{0, 0, 0, 0, 1, 0},
		{0, 0, 0, 0, 0, 1},
		{1, 0, 0, 0, 0, 0},
		{0, 1, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0}};

	e[0] = e[1] = e[2] = 0.0;
	for (int i=0; i<6; ++i)
	{
		e[0] += h1[i]*g[i];
		e[1] += h2[i]*g[i];
		e[2] += h3[i]*g[i];
		for (int j=0; j<6; ++j)
		{
			e[0] += 0.5*g[i]*H1[i][j]*g[j];
			e[1] += 0.5*g[i]*H2[i][j]*g[j];
			e[2] += 0.5*g[i]*H3[i][j]*g[j];
		}
	}
}

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEElasticMembrane, FEMembraneMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
void FEElasticMembrane::Stress(FEMaterialPoint &mp, double s[])
{
	FEMembraneMaterialPoint& pt = *mp.ExtractData<FEMembraneMaterialPoint>();

	// get the material point strain
	double e[3];
	pt.strain(e);

	// elasticity tangent
	double E = 1.0;
	double v = 0.35;
	double d = 1.0/(1.0 - v*v);
	double D[3][3] = {
		{d*E, d*v*E, 0},
		{d*v*E, d*E, 0},
		{0, 0, d*E*(1.0-v)*0.5}};

	// calculate 2D PK2-stress
	s[0] = D[0][0]*e[0] + D[0][1]*e[1] + D[0][2]*e[2];
	s[1] = D[1][0]*e[0] + D[1][1]*e[1] + D[1][2]*e[2];
	s[2] = D[2][0]*e[0] + D[2][1]*e[1] + D[2][2]*e[2];
}

//-----------------------------------------------------------------------------
void FEElasticMembrane::Tangent(FEMaterialPoint &mp, double D[][3])
{
	FEMembraneMaterialPoint& pt = *mp.ExtractData<FEMembraneMaterialPoint>();

	// elasticity tangent
	double E = 1.0;
	double v = 0.35;
	double d = 1.0/(1.0 - v*v);
	D[0][0] =   d*E; D[0][1] = d*v*E; D[0][2] = 0;
	D[1][0] = d*v*E; D[1][1] =   d*E; D[1][2] = 0;
	D[2][0] =     0; D[2][1] =     0; D[2][2] = d*E*(1.0-v)*0.5;
}
