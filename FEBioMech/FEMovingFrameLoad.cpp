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
	ADD_PARAMETER(m_wi.x, "wx")->setLongName("x-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_wi.y, "wy")->setLongName("y-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_wi.z, "wz")->setLongName("z-angular velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_at.x, "ax")->setLongName("x-linear acceleration")->setUnits(UNIT_ACCELERATION);
	ADD_PARAMETER(m_at.y, "ay")->setLongName("y-linear acceleration")->setUnits(UNIT_ACCELERATION);
	ADD_PARAMETER(m_at.z, "az")->setLongName("z-linear acceleration")->setUnits(UNIT_ACCELERATION);
	ADD_PARAMETER(m_frame, "frame")->setEnums("fixed\0moving\0");

	ADD_PARAMETER(m_wt, "wt")->setLongName("Angular velocity in fixed frame")->setUnits(UNIT_ANGULAR_VELOCITY)->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_Wt, "Wt")->setLongName("Angular velocity in moving frame")->setUnits(UNIT_ANGULAR_VELOCITY)->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_rt, "rt")->setLongName("Rotation vector in fixed frame")->setUnits(UNIT_RADIAN)->SetFlags(FE_PARAM_HIDDEN);
END_FECORE_CLASS();

FEMovingFrameLoad::FEMovingFrameLoad(FEModel* fem) : FEBodyForce(fem)
{
	m_frame = 0; // assume fixed frame angular velocity
}

void FEMovingFrameLoad::Activate()
{
	m_ap = m_at;
	m_wp = m_wt;
	m_Wp = m_Wt;
	m_alp = m_alt;
	m_Alp = m_Alt;
	m_qp = m_qt;
	FEBodyForce::Activate();
}

void FEMovingFrameLoad::PrepStep()
{
	const FETimeInfo& tp = GetTimeInfo();
	double dt = tp.timeIncrement;
	double alpha = tp.alphaf;

	// evaluate quantities at current time (t) and intermediate time
	if (m_frame == 0)
	{
		m_wt = m_wi; // input is ang vel in fixed frame

		// current quantities
		m_alp = m_alt;
		m_alt = (m_wt - m_wp) / dt;
		quatd wt(m_wt.x, m_wt.y, m_wt.z, 0.0);
		m_qt = m_qp + (wt * m_qp) * (0.5 * dt); m_qt.MakeUnit();

		// intermediate quantities
		m_al = m_alt * alpha + m_alp * (1.0 - alpha);
		m_w = m_wt * alpha + m_wp * (1.0 - alpha);
		quatd wa(m_w.x, m_w.y, m_w.z, 0.0);
		m_q = m_qp + (wa * m_qp) * (0.5 * dt * alpha); m_q.MakeUnit();

		m_qp = m_qt;
		m_wp = m_wt;

		// quantities in body frame
		quatd qti = m_qt.Inverse();
		m_Wt = qti * m_wt;
		m_rt = m_qt.GetRotationVector();

		quatd qi = m_q.Inverse();
		m_Al = qi * m_al;
		m_W = qi * m_w;
	}
	else
	{
		m_Wt = m_wi; // input is ang vel in body frame

		// current quantities
		m_Alp = m_Alt;
		m_Alt = (m_Wt - m_Wp) / dt;
		quatd Wt(m_Wt.x, m_Wt.y, m_Wt.z, 0.0);
		m_qt = m_qp + (m_qp * Wt) * (0.5 * dt); m_qt.MakeUnit();

		// intermediate quantities
		m_Al = m_Alt * alpha + m_Alp * (1.0 - alpha);
		m_W = m_Wt * alpha + m_Wp * (1.0 - alpha);
		quatd Wa(m_W.x, m_W.y, m_W.z, 0.0);
		m_q = m_qp + (m_qp * Wa) * (0.5 * dt * alpha); m_q.MakeUnit();

		m_qp = m_qt;
		m_Wp = m_Wt;

		// quantities in fixed frame
		m_wt = m_qt * m_Wt;
		m_rt = m_qt.GetRotationVector();
		m_alp = m_alt;
		m_alt = m_qt * m_Alt;

		m_al = m_q * m_Al;
		m_w  = m_q * m_W;
	}

	// linear acceleration (always assumed to be given in the fixed frame)
	m_a = m_at * alpha + m_ap * (1.0 - alpha);
	m_ap = m_at;
	quatd qi = m_q.Inverse();
	m_A = qi * m_a;

	FEBodyForce::PrepStep();
}

vec3d FEMovingFrameLoad::force(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();

	const FETimeInfo& tp = GetTimeInfo();
	double alpha = tp.alphaf;

	// NOTE: position and velocity at material point are already evaluated at intermediate time
	vec3d r = pt.m_rt;
	vec3d v = ep.m_v;

	vec3d f = m_A + (m_Al ^ r) + (m_W ^ (m_W ^ r)) + (m_W ^ v) * 2.0;
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

	mat3da Sw(m_W);
	mat3da A(m_Al);
	mat3d S2 = Sw * Sw;
	mat3d K = (S2 + A + Sw * (4.0 * gamma / dt))*alpha;

	return K;
}
