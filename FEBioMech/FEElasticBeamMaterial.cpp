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

FEElasticBeamMaterialPoint::FEElasticBeamMaterialPoint()
{

}

void FEElasticBeamMaterialPoint::Init()
{

}

void FEElasticBeamMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	m_Rp = m_Rt;
}

//=======================================================================================
BEGIN_FECORE_CLASS(FEElasticBeamMaterial, FEMaterial)
	ADD_PARAMETER(m_density, "density");
	ADD_PARAMETER(m_E , "E" );
	ADD_PARAMETER(m_G , "G" );
	ADD_PARAMETER(m_A , "A" );
	ADD_PARAMETER(m_A1, "A1");
	ADD_PARAMETER(m_A2, "A2");
	ADD_PARAMETER(m_I1, "I1");
	ADD_PARAMETER(m_I2, "I2");
END_FECORE_CLASS();

FEElasticBeamMaterial::FEElasticBeamMaterial(FEModel* fem) : FEMaterial(fem)
{
	m_density = 1.0;
	m_A = m_A1 = m_A2 = 0.0;
	m_G = m_E = 0.0;
	m_I1 = m_I2 = 0;
}

void FEElasticBeamMaterial::Stress(FEElasticBeamMaterialPoint& mp)
{
	vec3d gamma = mp.m_Gamma;
	vec3d kappa = mp.m_Kappa;
	quatd R = mp.m_Rt;
	double J = m_I1 + m_I2;

	// material vectors
	vec3d N = vec3d(m_G * m_A1 * gamma.x, m_G*m_A2*gamma.y, m_E*m_A*gamma.z);
	vec3d M = vec3d(m_E * m_I1 * kappa.x, m_E*m_I2*kappa.y, m_G*  J*kappa.z);

	// spatial vectors
	mp.m_t = R * N;
	mp.m_m = R * M;
}

void FEElasticBeamMaterial::Tangent(FEElasticBeamMaterialPoint& mp, matrix& C)
{
	double J = m_I1 + m_I2;

	mat3dd C1(m_G * m_A1, m_G * m_A2, m_E * m_A);
	mat3dd C2(m_E * m_I1, m_E * m_I2, m_G *   J);
	quatd R = mp.m_Rt;
	mat3d Q = R.RotationMatrix();

	mat3d c1 = Q * C1;
	mat3d c2 = Q * C2;

	C.resize(6, 6);
	C.zero();
	C.add(0, 0, c1);
	C.add(3, 3, c2);
}

FEMaterialPointData* FEElasticBeamMaterial::CreateMaterialPointData()
{
	return new FEElasticBeamMaterialPoint();
}
