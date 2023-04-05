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
#include "FEElasticBeamMaterial.h"

BEGIN_FECORE_CLASS(FEElasticBeamMaterial, FEMaterial)
	ADD_PARAMETER(m_density, "density");
	ADD_PARAMETER(m_E , "E" );
	ADD_PARAMETER(m_G , "G" );
	ADD_PARAMETER(m_A , "A" );
	ADD_PARAMETER(m_A1, "A1");
	ADD_PARAMETER(m_A2, "A2");
	ADD_PARAMETER(m_I1, "I1");
	ADD_PARAMETER(m_I2, "I2");
	ADD_PARAMETER(m_J , "J" );
END_FECORE_CLASS();

FEElasticBeamMaterial::FEElasticBeamMaterial(FEModel* fem) : FEMaterial(fem)
{
	m_density = 1.0;
	m_A = m_A1 = m_A2 = 0.0;
	m_G = m_E = 0.0;
	m_I1 = m_I2 = m_J = 0.0;
}

void FEElasticBeamMaterial::Stress(FEElasticBeamMaterialPoint& mp)
{
	// TODO: calculate these
	vec3d gamma;
	vec3d kappa;
	mat3d R;

	// material vectors
	vec3d N = vec3d(m_G * m_A1 * gamma.x, m_G*m_A2*gamma.y, m_E*m_A*gamma.z);
	vec3d M = vec3d(m_E * m_I1 * kappa.x, m_E*m_I2*kappa.y, m_G*m_J*gamma.z);

	// spatial vectors
	mp.m_t = R * N;
	mp.m_m = R * M;
}

void FEElasticBeamMaterial::Tangent(FEElasticBeamMaterialPoint& mp, matrix& C)
{
	C.resize(6, 6);
	C.zero();
	C[0][0] = m_G * m_A1;
	C[1][1] = m_G * m_A2;
	C[2][2] = m_E * m_A;
	C[3][3] = m_E * m_I1;
	C[4][4] = m_E * m_I2;
	C[5][5] = m_G * m_J;
}

FEMaterialPointData* FEElasticBeamMaterial::CreateMaterialPointData()
{
	return new FEElasticBeamMaterialPoint();
}
