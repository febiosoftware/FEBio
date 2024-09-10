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
#include <FECore/FEMaterial.h>
#include "febiomech_api.h"

class FEElasticBeamMaterialPoint : public FEMaterialPointData
{
public:
	FEElasticBeamMaterialPoint();

	void Init() override;

	//! The Update function is used to update material point data
	//! Note that this gets called at the start of the time step during PreSolveUpdate
	void Update(const FETimeInfo& timeInfo) override;

public:
	vec3d	m_t;	// stress vector
	vec3d	m_m;	// moment vector

	// local coordinate system
	mat3d	m_Q;

	// strain measures
	vec3d	m_G0;
	vec3d	m_Gamma;
	vec3d	m_Kappa;

	// rotation information
	quatd	m_Rp;	// rotation at previous time step
	quatd	m_Ri;	// increment at current time step
	quatd	m_Rt;	// current rotation

	vec3d	m_k;	// spatial curvature
	vec3d	m_kn;	// spatial curvature at current increment (temp storage)

	// dynamics
	vec3d	m_vt, m_vp;	// linear velocity at current and previous time
	vec3d	m_at, m_ap;	// linear acceleration at current and previous time
	vec3d	m_wt, m_wp;	// (spatial) angular velocity at current and previous time
	vec3d	m_alt, m_alp;	// (spatial) angular acceleration at current and previous time

	vec3d	m_dpt;	// rate of linear momentum (current)
	vec3d	m_dht;	// rate of angular momentum (current)

	// for output
	mat3ds m_s; // Cauchy stress tensor in global coordinates
};

class FEElasticBeamMaterial : public FEMaterial
{
public:
	FEElasticBeamMaterial(FEModel* fem);

	void Stress(FEElasticBeamMaterialPoint& mp);

	void Tangent(FEElasticBeamMaterialPoint& mp, matrix& C);

	FEMaterialPointData* CreateMaterialPointData() override;

public:
	double	m_density;
	double	m_A, m_A1, m_A2;
	double	m_E, m_G;
	double	m_I1, m_I2;

	DECLARE_FECORE_CLASS();
};

class FEBIOMECH_API FEBeamStress : public FEDomainParameter
{
public:
	FEBeamStress();
	FEParamValue value(FEMaterialPoint& mp) override;
};
