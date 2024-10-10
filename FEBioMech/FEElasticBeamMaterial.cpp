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
	m_Ri = quatd(0, 0, 0);

	m_vp = m_vt;
	m_ap = m_at;
	m_wp = m_wt;
	m_alp = m_alt;
}

//=======================================================================================
BEGIN_FECORE_CLASS(FEElasticBeamMaterial, FEMaterial)
	ADD_PARAMETER(m_density, "density")->setUnits(UNIT_DENSITY);
	ADD_PARAMETER(m_E , "E" )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_G , "G" )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_A , "A" )->setUnits(UNIT_AREA);
	ADD_PARAMETER(m_A1, "A1")->setUnits(UNIT_AREA);
	ADD_PARAMETER(m_A2, "A2")->setUnits(UNIT_AREA);
	ADD_PARAMETER(m_I1, "I1");
	ADD_PARAMETER(m_I2, "I2");
END_FECORE_CLASS();

FEElasticBeamMaterial::FEElasticBeamMaterial(FEModel* fem) : FEMaterial(fem)
{
	m_density = 1.0;
	m_A = m_A1 = m_A2 = 0.0;
	m_G = m_E = 0.0;
	m_I1 = m_I2 = 0;

	AddDomainParameter(new FEBeamStress());
}

void FEElasticBeamMaterial::Stress(FEElasticBeamMaterialPoint& mp)
{
	mat3d Q = mp.m_Q;
	mat3d Qt = Q.transpose();

	vec3d gamma = Qt * mp.m_Gamma;
	vec3d kappa = Qt * mp.m_Kappa;
	double J = m_I1 + m_I2;

	// material vectors
	vec3d N = Q*vec3d(m_G * m_A1 * gamma.x, m_G*m_A2*gamma.y, m_E*m_A*gamma.z);
	vec3d M = Q*vec3d(m_E * m_I1 * kappa.x, m_E*m_I2*kappa.y, m_G*  J*kappa.z);

	// spatial vectors
	quatd R = mp.m_Rt;
	mp.m_t = R * N;
	mp.m_m = R * M;

	// Cauchy stress
	vec3d E1 = Q.col(0);
	vec3d E2 = Q.col(1);
	vec3d E3 = Q.col(2);
	vec3d t1 = R * E1;
	vec3d t2 = R * E2;
	vec3d t3 = R * E3;
	vec3d t = mp.m_t;
	mp.m_s = (dyad(t1) * (t * t1) + dyad(t2) * (t * t2) + dyad(t3) * (t * t3));
}

void FEElasticBeamMaterial::Tangent(FEElasticBeamMaterialPoint& mp, matrix& C)
{
	mat3d E = mp.m_Q;
	mat3d Et = E.transpose();

	double J = m_I1 + m_I2;

	mat3dd D1(m_G * m_A1, m_G * m_A2, m_E * m_A);
	mat3dd D2(m_E * m_I1, m_E * m_I2, m_G *   J);

	mat3d C1 = E * D1 * Et;
	mat3d C2 = E * D2 * Et;

	quatd R = mp.m_Rt;
	mat3d Q = R.RotationMatrix();
	mat3d Qt = Q.transpose();

	mat3d c1 = Q * C1 * Qt;
	mat3d c2 = Q * C2 * Qt;

	C.resize(6, 6);
	C.zero();
	C.add(0, 0, c1);
	C.add(3, 3, c2);
}

FEMaterialPointData* FEElasticBeamMaterial::CreateMaterialPointData()
{
	return new FEElasticBeamMaterialPoint();
}

FEBeamStress::FEBeamStress() : FEDomainParameter("stress") {}

FEParamValue FEBeamStress::value(FEMaterialPoint& mp)
{
	FEElasticBeamMaterialPoint& pt = *mp.ExtractData<FEElasticBeamMaterialPoint>();
	return pt.m_s;
}
