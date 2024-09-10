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
	ADD_PARAMETER(m_wt.x, "wx")->setLongName("x-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_wt.y, "wy")->setLongName("y-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_wt.z, "wz")->setLongName("z-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_at.x, "ax")->setLongName("x-linear acceleration")->setUnits(UNIT_ACCELERATION);
	ADD_PARAMETER(m_at.y, "ay")->setLongName("y-linear acceleration")->setUnits(UNIT_ACCELERATION);
	ADD_PARAMETER(m_at.z, "az")->setLongName("z-linear acceleration")->setUnits(UNIT_ACCELERATION);

	ADD_PARAMETER(m_wt, "wt")->setLongName("Angular velocity in fixed frame")->setUnits(UNIT_ANGULAR_VELOCITY)->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_Wt, "Wt")->setLongName("Angular velocity in moving frame")->setUnits(UNIT_ANGULAR_VELOCITY)->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_rt, "rt")->setLongName("Rotation vector in fixed frame")->setUnits(UNIT_RADIAN)->SetFlags(FE_PARAM_HIDDEN);
END_FECORE_CLASS();

FEMovingFrameLoad::FEMovingFrameLoad(FEModel* fem) : FEBodyForce(fem)
{
}

void FEMovingFrameLoad::Activate()
{
	m_wp = m_wt;
	m_ap = m_at;
	m_alp = m_alt;
	m_qp = m_qt;
	FEBodyForce::Activate();
}

void FEMovingFrameLoad::PrepStep()
{
	const FETimeInfo& tp = GetTimeInfo();
	double dt = tp.timeIncrement;
	double alpha = tp.alphaf;

	// evaluate quantities at current time (t) and intermediate time
	m_a = m_at * alpha + m_ap * (1.0 - alpha);
	m_ap = m_at;

	m_alp = m_alt;
	m_alt = (m_wt - m_wp) / dt;
	m_al = m_alt * alpha + m_alp * (1.0 - alpha);

	m_w = m_wt * alpha + m_wp * (1.0 - alpha);

	quatd wa(m_w.x, m_w.y, m_w.z, 0.0);
	quatd wt(m_wt.x, m_wt.y, m_wt.z, 0.0);

	m_qt = m_qp + (wt * m_qp) * (0.5 * dt); m_qt.MakeUnit();
	m_q = m_qp + (wa * m_qp) * (0.5 * dt * alpha); m_q.MakeUnit();

	m_qp = m_qt;
	m_wp = m_wt;

	// output quantities
	quatd qi = m_qt.Inverse();
	m_Wt = qi * m_wt;
	m_rt = m_qt.GetRotationVector();

	vec3d R = m_q.GetRotationVector();

	feLogDebug("al : %lg, %lg, %lg\n", m_al.x, m_al.y, m_al.z);
	feLogDebug("R : %lg, %lg, %lg\n", R.x, R.y, R.z);
}

vec3d FEMovingFrameLoad::force(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();

	const FETimeInfo& tp = GetTimeInfo();
	double alpha = tp.alphaf;

	// NOTE: position and velocity at material point are already evaluated at intermediate time
	vec3d r = pt.m_rt;
	vec3d v = ep.m_v;

	quatd qi = m_q.Inverse();
	vec3d w = qi * m_w;
	vec3d al = qi * m_al;

	vec3d f = qi*m_a + (al ^ r) + (w ^ (w ^ r)) + (w ^ v) * 2.0;
	return f;
}

mat3d FEMovingFrameLoad::stiffness(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();
	vec3d v = ep.m_v;

	const FETimeInfo& tp = GetTimeInfo();
	double gamma = tp.gamma;
	double dt = tp.timeIncrement;
	double alpha = tp.alphaf;

	quatd qi = m_q.Inverse();
	vec3d w = qi * m_w;
	vec3d al = qi * m_al;

	mat3da Sw(w);
	mat3da A(al);
	mat3d S2 = Sw * Sw;
	mat3d K = (S2 + A + Sw * (4.0 * gamma / dt))*alpha;

	return K;
}
