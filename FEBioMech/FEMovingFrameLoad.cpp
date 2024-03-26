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
#include "FEMovingFrameLoad.h"
#include "FEElasticMaterialPoint.h"
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEMovingFrameLoad, FEBodyForce)
	ADD_PARAMETER(m_w.x, "wx")->setLongName("x-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_w.y, "wy")->setLongName("y-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_w.z, "wz")->setLongName("z-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_a.x, "ax")->setLongName("x-linear acceleration")->setUnits(UNIT_ACCELERATION);
	ADD_PARAMETER(m_a.y, "ay")->setLongName("y-linear acceleration")->setUnits(UNIT_ACCELERATION);
	ADD_PARAMETER(m_a.z, "az")->setLongName("z-linear acceleration")->setUnits(UNIT_ACCELERATION);
END_FECORE_CLASS();

FEMovingFrameLoad::FEMovingFrameLoad(FEModel* fem) : FEBodyForce(fem)
{
	m_w = vec3d(0, 0, 0);
	m_al = vec3d(0, 0, 0);
	m_q = quatd(0, vec3d(0, 0, 1));
}

bool FEMovingFrameLoad::Init()
{
	m_wp = m_w;
	return FEBodyForce::Init();
}

void FEMovingFrameLoad::PrepStep()
{
	double dt = GetTimeInfo().timeIncrement;
	m_al = (m_w - m_wp) / dt;
	m_wp = m_w;

	quatd qn = m_q;
	quatd wn(m_w.x, m_w.y, m_w.z, 0.0);
	m_q = qn + (wn * qn) * (0.5 * dt);
	m_q.MakeUnit();

	vec3d R = m_q.GetRotationVector();

	feLogDebug("al : %lg, %lg, %lg\n", m_al.x, m_al.y, m_al.z);
	feLogDebug("R : %lg, %lg, %lg\n", R.x, R.y, R.z);
}

vec3d FEMovingFrameLoad::force(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();

	vec3d r = pt.m_rt;
	vec3d v = ep.m_v;
	vec3d f = m_a + (m_al ^ r) + (m_w ^ (m_w ^ r)) + (m_w ^ v) * 2.0;
	return f;
}

mat3d FEMovingFrameLoad::stiffness(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();
	vec3d v = ep.m_v;

	const FETimeInfo& tp = GetTimeInfo();
	double gamma = tp.gamma;
	double dt = tp.timeIncrement;

	mat3da Sw(m_w);
	mat3da A(m_al);
	mat3d S2 = Sw * Sw;
	mat3d K = S2 + A + Sw * (4.0 * gamma / dt);

	// NOTE: we need a negative sign, since the external load stiffness needs to be subtracted from the global stiffness
	return -K;
}
