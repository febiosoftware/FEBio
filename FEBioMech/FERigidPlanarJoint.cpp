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
#include "FERigidPlanarJoint.h"
#include "FERigidBody.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>
#include <FECore/ad.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidPlanarJoint, FERigidConnector);
	ADD_PARAMETER(m_laugon, "laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0LAGMULT\0");
	ADD_PARAMETER(m_atol, "tolerance"     );
	ADD_PARAMETER(m_gtol, "gaptol"        );
	ADD_PARAMETER(m_qtol, "angtol"        );
	ADD_PARAMETER(m_eps , "force_penalty" );
	ADD_PARAMETER(m_ups , "moment_penalty");
	ADD_PARAMETER(m_q0  , "joint_origin"  );
	ADD_PARAMETER(m_e0[0], "rotation_axis" );
	ADD_PARAMETER(m_e0[1], "translation_axis_1");
	ADD_PARAMETER(m_naugmin, "minaug"        );
	ADD_PARAMETER(m_naugmax, "maxaug"        );
	ADD_PARAMETER(m_bqx , "prescribed_rotation");
	ADD_PARAMETER(m_qpx , "rotation"      );
	ADD_PARAMETER(m_bdy , "prescribed_translation_1");
	ADD_PARAMETER(m_dpy , "translation_1" );
	ADD_PARAMETER(m_bdz , "prescribed_translation_2");
	ADD_PARAMETER(m_dpz , "translation_2" );
	ADD_PARAMETER(m_bautopen, "auto_penalty");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidPlanarJoint::FERigidPlanarJoint(FEModel* pfem) : FERigidConnector(pfem)
{
    m_nID = m_ncount++;

	m_laugon = FECore::AUGLAG_METHOD; // for backward compatibility
    m_atol = 0;
    m_gtol = 0;
    m_qtol = 0;
    m_naugmin = 0;
    m_naugmax = 10;
    m_qpx = 0;
    m_dpy = 0;
    m_dpz = 0;
    m_bqx = false;
    m_bdy = false;
    m_bdz = false;
	m_bautopen = false;
    m_eps = m_ups = 1.0;
}

//-----------------------------------------------------------------------------
//! initial position
vec3d FERigidPlanarJoint::InitialPosition() const
{
    return m_q0;
}

//-----------------------------------------------------------------------------
//! current position
vec3d FERigidPlanarJoint::Position() const
{
    FERigidBody& RBa = *m_rbA;
    vec3d qa = m_qa0;
    RBa.GetRotation().RotateVector(qa);
    return RBa.m_rt + qa;
}

//-----------------------------------------------------------------------------
//! current orientation
quatd FERigidPlanarJoint::Orientation() const
{
	quatd Q0(vec3d(0, 0, 1), m_e0[0]);
	FERigidBody& RBa = *m_rbA;
	return RBa.GetRotation()*Q0;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidPlanarJoint::Init()
{
    // initialize joint basis
    m_e0[0].unit();
    m_e0[2] = m_e0[0] ^ m_e0[1]; m_e0[2].unit();
    m_e0[1] = m_e0[2] ^ m_e0[0]; m_e0[1].unit();
    
    // reset force
    m_F = vec3d(0,0,0); m_L = vec3d(0,0,0);
    m_M = vec3d(0,0,0); m_U = vec3d(0,0,0);
    
	// base class first
	if (FERigidConnector::Init() == false) return false;

    m_qa0 = m_q0 - m_rbA->m_r0;
    m_qb0 = m_q0 - m_rbB->m_r0;
    
    m_ea0[0] = m_e0[0]; m_ea0[1] = m_e0[1]; m_ea0[2] = m_e0[2];
    m_eb0[0] = m_e0[0]; m_eb0[1] = m_e0[1]; m_eb0[2] = m_e0[2];
    
    return true;
}

//-----------------------------------------------------------------------------
void FERigidPlanarJoint::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
	ar & m_qa0 & m_qb0;
	ar & m_L & m_U;
	ar & m_e0;
	ar & m_ea0;
	ar & m_eb0;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidPlanarJoint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<double> fa(6);
    vector<double> fb(6);
    
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    double alpha = tp.alphaf;
    
    // body A
    vec3d ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
	eat[0] = m_ea0[0]; RBa.GetRotation().RotateVector(eat[0]);
	eat[1] = m_ea0[1]; RBa.GetRotation().RotateVector(eat[1]);
	eat[2] = m_ea0[2]; RBa.GetRotation().RotateVector(eat[2]);
    eap[0] = m_ea0[0]; RBa.m_qp.RotateVector(eap[0]);
    eap[1] = m_ea0[1]; RBa.m_qp.RotateVector(eap[1]);
    eap[2] = m_ea0[2]; RBa.m_qp.RotateVector(eap[2]);
    ea[0] = eat[0]*alpha + eap[0]*(1-alpha);
    ea[1] = eat[1]*alpha + eap[1]*(1-alpha);
    ea[2] = eat[2]*alpha + eap[2]*(1-alpha);
    
    // body b
    vec3d rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
	ebt[0] = m_eb0[0]; RBb.GetRotation().RotateVector(ebt[0]);
	ebt[1] = m_eb0[1]; RBb.GetRotation().RotateVector(ebt[1]);
	ebt[2] = m_eb0[2]; RBb.GetRotation().RotateVector(ebt[2]);
    ebp[0] = m_eb0[0]; RBb.m_qp.RotateVector(ebp[0]);
    ebp[1] = m_eb0[1]; RBb.m_qp.RotateVector(ebp[1]);
    ebp[2] = m_eb0[2]; RBb.m_qp.RotateVector(ebp[2]);
    eb[0] = ebt[0]*alpha + ebp[0]*(1-alpha);
    eb[1] = ebt[1]*alpha + ebp[1]*(1-alpha);
    eb[2] = ebt[2]*alpha + ebp[2]*(1-alpha);

	if (m_laugon != FECore::LAGMULT_METHOD)
	{
		mat3ds P;
		vec3d p;
		if (m_bdy && m_bdz) {
			P = mat3dd(1);
			p = ea[1] * m_dpy + ea[2] * m_dpz;
		}
		else if (m_bdy) {
			P = mat3dd(1) - dyad(ea[2]);
			p = ea[1] * m_dpy;
		}
		else if (m_bdz) {
			P = mat3dd(1) - dyad(ea[1]);
			p = ea[2] * m_dpz;
		}
		else {
			P = dyad(ea[0]);
			p = vec3d(0, 0, 0);
		}
		vec3d c = P * (rb + zb - ra - za) - p;
		m_F = m_L + c * m_eps;

		vec3d ksi;
		if (m_bqx) {
			quatd q = (alpha * RBb.GetRotation() + (1 - alpha) * RBb.m_qp) * (alpha * RBa.GetRotation() + (1 - alpha) * RBa.m_qp).Inverse();
			quatd a(m_qpx, ea[0]);
			quatd r = a * q.Inverse();
			r.MakeUnit();
			ksi = r.GetVector() * r.GetAngle();
		}
		else
		{
			ksi = (ea[0] ^ eb[0]) / 2;
		}
		m_M = m_U + ksi * m_ups;

		fa[0] = m_F.x;
		fa[1] = m_F.y;
		fa[2] = m_F.z;

		fa[3] = za.y * m_F.z - za.z * m_F.y + m_M.x;
		fa[4] = za.z * m_F.x - za.x * m_F.z + m_M.y;
		fa[5] = za.x * m_F.y - za.y * m_F.x + m_M.z;

		fb[0] = -m_F.x;
		fb[1] = -m_F.y;
		fb[2] = -m_F.z;

		fb[3] = -zb.y * m_F.z + zb.z * m_F.y - m_M.x;
		fb[4] = -zb.z * m_F.x + zb.x * m_F.z - m_M.y;
		fb[5] = -zb.x * m_F.y + zb.y * m_F.x - m_M.z;

		for (int i = 0; i < 6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
		for (int i = 0; i < 6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];
	}
	else
	{
		vec3d d = rb + zb - ra - za;
		mat3dd I(1.0);
		mat3d P = dyad(ea[0]);
		mat3d PdT = I * (ea[0] * d) + (d & ea[0]);

		vec3d F = P*m_F;
		vec3d Ma = (za ^ F) - (ea[0] ^ (PdT*m_F)) + (ea[0]^m_M);
		vec3d Mb = (zb ^ F) + (eb[0] ^ m_M);

		fa[0] = -F.x;
		fa[1] = -F.y;
		fa[2] = -F.z;

		fa[3] = -Ma.x;
		fa[4] = -Ma.y;
		fa[5] = -Ma.z;

		fb[0] = F.x;
		fb[1] = F.y;
		fb[2] = F.z;

		fb[3] = Mb.x;
		fb[4] = Mb.y;
		fb[5] = Mb.z;

		for (int i = 0; i < 6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
		for (int i = 0; i < 6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];

		// translational constraint
		vec3d c = ea[0]*(d*ea[0]);
		R[m_LM[0]] += c.x;
		R[m_LM[1]] += c.y;
		R[m_LM[2]] += c.z;

		// rotational constraint
		vec3d ksi = eb[0] - ea[0];
		R[m_LM[3]] += ksi.x;
		R[m_LM[4]] += ksi.y;
		R[m_LM[5]] += ksi.z;
	}
    
    RBa.m_Fr -= vec3d(fa[0],fa[1],fa[2]);
    RBa.m_Mr -= vec3d(fa[3],fa[4],fa[5]);
    RBb.m_Fr -= vec3d(fb[0],fb[1],fb[2]);
    RBb.m_Mr -= vec3d(fb[3],fb[4],fb[5]);
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidPlanarJoint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    double alpha = tp.alphaf;
    
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
    vector<int> LM(12);
	FEElementMatrix ke; ke.resize(12, 12);
    ke.zero();
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    // body A
	quatd Qat = RBa.GetRotation();
	quatd Qap = RBa.m_qp;
	vec3d rat = RBa.m_rt;
	vec3d rap = RBa.m_rp;
	vec3d ra = rat * alpha + rap * (1 - alpha);
	vec3d zat = Qat * m_qa0;
	vec3d zap = Qap * m_qa0;
	vec3d za = zat * alpha + zap * (1 - alpha);
	eat[0] = m_ea0[0]; RBa.GetRotation().RotateVector(eat[0]);
	eat[1] = m_ea0[1]; RBa.GetRotation().RotateVector(eat[1]);
	eat[2] = m_ea0[2]; RBa.GetRotation().RotateVector(eat[2]);
    eap[0] = m_ea0[0]; RBa.m_qp.RotateVector(eap[0]);
    eap[1] = m_ea0[1]; RBa.m_qp.RotateVector(eap[1]);
    eap[2] = m_ea0[2]; RBa.m_qp.RotateVector(eap[2]);
    ea[0] = eat[0]*alpha + eap[0]*(1-alpha);
    ea[1] = eat[1]*alpha + eap[1]*(1-alpha);
    ea[2] = eat[2]*alpha + eap[2]*(1-alpha);
    mat3d zahat; zahat.skew(za);
    mat3d zathat; zathat.skew(zat);
    
    // body b
	quatd Qbt = RBb.GetRotation();
	quatd Qbp = RBb.m_qp;
	vec3d rbt = RBb.m_rt;
	vec3d rbp = RBb.m_rp;
	vec3d rb = rbt * alpha + rbp * (1 - alpha);
	vec3d zbt = Qbt * m_qb0;
	vec3d zbp = Qbp * m_qb0;
	vec3d zb = zbt * alpha + zbp * (1 - alpha);
	ebt[0] = m_eb0[0]; RBb.GetRotation().RotateVector(ebt[0]);
	ebt[1] = m_eb0[1]; RBb.GetRotation().RotateVector(ebt[1]);
	ebt[2] = m_eb0[2]; RBb.GetRotation().RotateVector(ebt[2]);
    ebp[0] = m_eb0[0]; RBb.m_qp.RotateVector(ebp[0]);
    ebp[1] = m_eb0[1]; RBb.m_qp.RotateVector(ebp[1]);
    ebp[2] = m_eb0[2]; RBb.m_qp.RotateVector(ebp[2]);
    eb[0] = ebt[0]*alpha + ebp[0]*(1-alpha);
    eb[1] = ebt[1]*alpha + ebp[1]*(1-alpha);
    eb[2] = ebt[2]*alpha + ebp[2]*(1-alpha);
    mat3d zbhat; zbhat.skew(zb);
    mat3d zbthat; zbthat.skew(zbt);
    
    mat3d eahat[3], ebhat[3], eathat[3], ebthat[3];
    for (int j=0; j<3; ++j) {
        eahat[j] = skew(ea[j]);
        ebhat[j] = skew(eb[j]);
        eathat[j] = skew(eat[j]);
        ebthat[j] = skew(ebt[j]);
    }

	if (m_laugon != FECore::LAGMULT_METHOD)
	{
		mat3ds P;
		vec3d p;
		mat3d Q, Wba, Wab;
		vec3d d = rb + zb - ra - za;
		if (m_bdy && m_bdz) {
			P = mat3dd(1);
			p = ea[1] * m_dpy + ea[2] * m_dpz;
			Q = mat3dd(0);
		}
		else if (m_bdy) {
			P = mat3dd(1) - dyad(ea[2]);
			p = ea[1] * m_dpy;
			Q = ((ea[2] & d) + mat3dd(1) * (ea[2] * d)) * eathat[2];
		}
		else if (m_bdz) {
			P = mat3dd(1) - dyad(ea[1]);
			p = ea[2] * m_dpz;
			Q = ((ea[1] & d) + mat3dd(1) * (ea[1] * d)) * eathat[1];
		}
		else {
			P = dyad(ea[0]);
			p = vec3d(0, 0, 0);
			Q = ((ea[0] & d) + mat3dd(1) * (ea[0] * d)) * eathat[0] * (-1);
		}
		vec3d c = P * d - p;
		m_F = m_L + c * m_eps;

		vec3d ksi;
		if (m_bqx) {
			quatd q = (alpha * RBb.GetRotation() + (1 - alpha) * RBb.m_qp) * (alpha * RBa.GetRotation() + (1 - alpha) * RBa.m_qp).Inverse();
			quatd a(m_qpx, ea[0]);
			quatd r = a * q.Inverse();
			r.MakeUnit();
			ksi = r.GetVector() * r.GetAngle();
			quatd qa = RBa.GetRotation() * (alpha * RBa.GetRotation() + (1 - alpha) * RBa.m_qp).Inverse();
			quatd qb = RBb.GetRotation() * (alpha * RBb.GetRotation() + (1 - alpha) * RBb.m_qp).Inverse();
			qa.MakeUnit();
			qb.MakeUnit();
			mat3d Qa = qa.RotationMatrix();
			mat3d Qb = qb.RotationMatrix();
			mat3d A = a.RotationMatrix();
			mat3d R = r.RotationMatrix();
			mat3dd I(1);
			Wba = A * (I * Qa.trace() - Qa) / 2;
			Wab = R * (I * Qb.trace() - Qb) / 2;
		}
		else
		{
			ksi = (ea[0] ^ eb[0]) / 2;
			Wba = (ebhat[0] * eathat[0]) / 2;
			Wab = (eahat[0] * ebthat[0]) / 2;
		}
		m_M = m_U + ksi * m_ups;

		mat3da Fhat(m_F);

		mat3d K;

		// (1,1)
		K = P*(alpha*m_eps);
		ke[0][0] = K[0][0]; ke[0][1] = K[0][1]; ke[0][2] = K[0][2];
		ke[1][0] = K[1][0]; ke[1][1] = K[1][1]; ke[1][2] = K[1][2];
		ke[2][0] = K[2][0]; ke[2][1] = K[2][1]; ke[2][2] = K[2][2];

		// (1,2)
		K = (P*zathat+Q)*(-m_eps*alpha);
		ke[0][3] = K[0][0]; ke[0][4] = K[0][1]; ke[0][5] = K[0][2];
		ke[1][3] = K[1][0]; ke[1][4] = K[1][1]; ke[1][5] = K[1][2];
		ke[2][3] = K[2][0]; ke[2][4] = K[2][1]; ke[2][5] = K[2][2];

		// (1,3)
		K = P*(-alpha*m_eps);
		ke[0][6] = K[0][0]; ke[0][7] = K[0][1]; ke[0][8] = K[0][2];
		ke[1][6] = K[1][0]; ke[1][7] = K[1][1]; ke[1][8] = K[1][2];
		ke[2][6] = K[2][0]; ke[2][7] = K[2][1]; ke[2][8] = K[2][2];

		// (1,4)
		K = P*zbthat*(alpha*m_eps);
		ke[0][9] = K[0][0]; ke[0][10] = K[0][1]; ke[0][11] = K[0][2];
		ke[1][9] = K[1][0]; ke[1][10] = K[1][1]; ke[1][11] = K[1][2];
		ke[2][9] = K[2][0]; ke[2][10] = K[2][1]; ke[2][11] = K[2][2];

		// (2,1)
		K = zahat*P*(alpha*m_eps);
		ke[3][0] = K[0][0]; ke[3][1] = K[0][1]; ke[3][2] = K[0][2];
		ke[4][0] = K[1][0]; ke[4][1] = K[1][1]; ke[4][2] = K[1][2];
		ke[5][0] = K[2][0]; ke[5][1] = K[2][1]; ke[5][2] = K[2][2];

		// (2,2)
		K = (zahat * (P * zathat + Q) * m_eps + Fhat * zathat + Wba * m_ups) * (-alpha);
		ke[3][3] = K[0][0]; ke[3][4] = K[0][1]; ke[3][5] = K[0][2];
		ke[4][3] = K[1][0]; ke[4][4] = K[1][1]; ke[4][5] = K[1][2];
		ke[5][3] = K[2][0]; ke[5][4] = K[2][1]; ke[5][5] = K[2][2];

		// (2,3)
		K = zahat*P*(-alpha*m_eps);
		ke[3][6] = K[0][0]; ke[3][7] = K[0][1]; ke[3][8] = K[0][2];
		ke[4][6] = K[1][0]; ke[4][7] = K[1][1]; ke[4][8] = K[1][2];
		ke[5][6] = K[2][0]; ke[5][7] = K[2][1]; ke[5][8] = K[2][2];

		// (2,4)
		K = (zahat*P*zbthat*m_eps + Wab*m_ups)*alpha;
		ke[3][9] = K[0][0]; ke[3][10] = K[0][1]; ke[3][11] = K[0][2];
		ke[4][9] = K[1][0]; ke[4][10] = K[1][1]; ke[4][11] = K[1][2];
		ke[5][9] = K[2][0]; ke[5][10] = K[2][1]; ke[5][11] = K[2][2];


		// (3,1)
		K = P*(-alpha*m_eps);
		ke[6][0] = K[0][0]; ke[6][1] = K[0][1]; ke[6][2] = K[0][2];
		ke[7][0] = K[1][0]; ke[7][1] = K[1][1]; ke[7][2] = K[1][2];
		ke[8][0] = K[2][0]; ke[8][1] = K[2][1]; ke[8][2] = K[2][2];

		// (3,2)
		K = (P*zathat+Q)*(m_eps*alpha);
		ke[6][3] = K[0][0]; ke[6][4] = K[0][1]; ke[6][5] = K[0][2];
		ke[7][3] = K[1][0]; ke[7][4] = K[1][1]; ke[7][5] = K[1][2];
		ke[8][3] = K[2][0]; ke[8][4] = K[2][1]; ke[8][5] = K[2][2];

		// (3,3)
		K = P*(alpha*m_eps);
		ke[6][6] = K[0][0]; ke[6][7] = K[0][1]; ke[6][8] = K[0][2];
		ke[7][6] = K[1][0]; ke[7][7] = K[1][1]; ke[7][8] = K[1][2];
		ke[8][6] = K[2][0]; ke[8][7] = K[2][1]; ke[8][8] = K[2][2];

		// (3,4)
		K = P*zbthat*(-alpha*m_eps);
		ke[6][9] = K[0][0]; ke[6][10] = K[0][1]; ke[6][11] = K[0][2];
		ke[7][9] = K[1][0]; ke[7][10] = K[1][1]; ke[7][11] = K[1][2];
		ke[8][9] = K[2][0]; ke[8][10] = K[2][1]; ke[8][11] = K[2][2];


		// (4,1)
		K = zbhat*P*(-alpha*m_eps);
		ke[9 ][0] = K[0][0]; ke[ 9][1] = K[0][1]; ke[ 9][2] = K[0][2];
		ke[10][0] = K[1][0]; ke[10][1] = K[1][1]; ke[10][2] = K[1][2];
		ke[11][0] = K[2][0]; ke[11][1] = K[2][1]; ke[11][2] = K[2][2];

		// (4,2)
		K = (zbhat*(P*zathat+Q)*m_eps + Wba*m_ups)*alpha;
		ke[9 ][3] = K[0][0]; ke[ 9][4] = K[0][1]; ke[ 9][5] = K[0][2];
		ke[10][3] = K[1][0]; ke[10][4] = K[1][1]; ke[10][5] = K[1][2];
		ke[11][3] = K[2][0]; ke[11][4] = K[2][1]; ke[11][5] = K[2][2];

		// (4,3)
		K = zbhat * P * (alpha * m_eps);
		ke[9 ][6] = K[0][0]; ke[ 9][7] = K[0][1]; ke[ 9][8] = K[0][2];
		ke[10][6] = K[1][0]; ke[10][7] = K[1][1]; ke[10][8] = K[1][2];
		ke[11][6] = K[2][0]; ke[11][7] = K[2][1]; ke[11][8] = K[2][2];

		// (4,4)
		K = (zbhat*P*zbthat*m_eps - Fhat * zbthat + Wab * m_ups) * (-alpha);
		ke[9 ][9] = K[0][0]; ke[ 9][10] = K[0][1]; ke[ 9][11] = K[0][2];
		ke[10][9] = K[1][0]; ke[10][10] = K[1][1]; ke[10][11] = K[1][2];
		ke[11][9] = K[2][0]; ke[11][10] = K[2][1]; ke[11][11] = K[2][2];

		for (int j = 0; j < 6; ++j)
		{
			LM[j] = RBa.m_LM[j];
			LM[j + 6] = RBb.m_LM[j];
		}
		ke.SetIndices(LM);
		LS.Assemble(ke);
	}
	else
	{
		dd::vec3d pa(ra);
		dd::quatd qa(Qat);
		dd::vec3d pb(rb);
		dd::quatd qb(Qbt);
		dd::vec3d lam(m_F);
		dd::vec3d mu(m_M);

		auto F = [&]() {
			dd::vec3d a = qa * m_ea0[0];
			return a*(a*lam);
			};

		auto c = [&]() {
			dd::vec3d za = qa * m_qa0;
			dd::vec3d zb = qb * m_qb0;
			dd::vec3d d = pb + zb - pa - za;
			dd::vec3d a = qa * m_ea0[0];
			return a * (d * a);
			};

		auto m1 = [&]() {
			dd::vec3d a = qa * m_ea0[0];
			dd::vec3d b = qb * m_eb0[0];
			return b - a;
			};

		auto Ma = [&] {
			dd::vec3d za = qa * m_qa0;
			dd::vec3d ea = qa * m_ea0[0];
			dd::vec3d d = pb + zb - pa - za;
			dd::vec3d F = ea * (lam * ea);
			return (za ^ F) - (ea ^ (d*(ea*lam) + lam*(ea*d))) + (ea ^ mu);
			};

		auto Mb = [&] {
			dd::vec3d zb = qb * m_qb0;
			dd::vec3d ea = qa * m_ea0[0];
			dd::vec3d eb = qb * m_eb0[0];
			dd::vec3d F = ea * (lam * ea);
			return (zb ^ F) + (eb ^ mu);
			};

		mat3d K00 = -dd::D( F, pa);
		mat3d K10 = -dd::D(Ma, pa);
		mat3d K20 =  dd::D( F, pa);
		mat3d K30 =  dd::D(Mb, pa);
		mat3d K40 =  dd::D( c, pa);
		mat3d K50 =  dd::D(m1, pa);

		mat3d K01 = -dd::D( F, qa);
		mat3d K11 = -dd::D(Ma, qa);
		mat3d K21 =  dd::D( F, qa);
		mat3d K31 =  dd::D(Mb, qa);
		mat3d K41 =  dd::D( c, qa);
		mat3d K51 =  dd::D(m1, qa);

		mat3d K02 = -dd::D( F, pb);
		mat3d K12 = -dd::D(Ma, pb);
		mat3d K22 =  dd::D( F, pb);
		mat3d K32 =  dd::D(Mb, pb);
		mat3d K42 =  dd::D( c, pb);
		mat3d K52 =  dd::D(m1, pb);

		mat3d K03 = -dd::D( F, qb);
		mat3d K13 = -dd::D(Ma, qb);
		mat3d K23 =  dd::D( F, qb);
		mat3d K33 =  dd::D(Mb, qb);
		mat3d K43 =  dd::D( c, qb);
		mat3d K53 =  dd::D(m1, qb);

		mat3d K04 = -dd::D( F, lam);
		mat3d K14 = -dd::D(Ma, lam);
		mat3d K24 =  dd::D( F, lam);
		mat3d K34 =  dd::D(Mb, lam);
		mat3d K44 =  dd::D( c, lam);
		mat3d K54 =  dd::D(m1, lam);

		mat3d K05 = -dd::D( F, mu);
		mat3d K15 = -dd::D(Ma, mu);
		mat3d K25 =  dd::D( F, mu);
		mat3d K35 =  dd::D(Mb, mu);
		mat3d K45 =  dd::D( c, mu);
		mat3d K55 =  dd::D(m1, mu);

		matrix kd(18, 18); kd.zero();
		kd.sub( 0, 0, K00); kd.sub( 0, 3, K01); kd.sub( 0, 6, K02); kd.sub( 0, 9, K03); kd.sub( 0, 12, K04); kd.sub( 0, 15, K05);
		kd.sub( 3, 0, K10); kd.sub( 3, 3, K11); kd.sub( 3, 6, K12); kd.sub( 3, 9, K13); kd.sub( 3, 12, K14); kd.sub( 3, 15, K15);
		kd.sub( 6, 0, K20); kd.sub( 6, 3, K21); kd.sub( 6, 6, K22); kd.sub( 6, 9, K23); kd.sub( 6, 12, K24); kd.sub( 6, 15, K25);
		kd.sub( 9, 0, K30); kd.sub( 9, 3, K31); kd.sub( 9, 6, K32); kd.sub( 9, 9, K33); kd.sub( 9, 12, K34); kd.sub( 9, 15, K35);
		kd.sub(12, 0, K40); kd.sub(12, 3, K41); kd.sub(12, 6, K42); kd.sub(12, 9, K43); kd.sub(12, 12, K44); kd.sub(12, 15, K45);
		kd.sub(15, 0, K50); kd.sub(15, 3, K51); kd.sub(15, 6, K52); kd.sub(15, 9, K53); kd.sub(15, 12, K54); kd.sub(15, 15, K55);

		FEElementMatrix ke(18, 18);
		ke = kd;
		vector<int> LM(18);
		UnpackLM(LM);
		ke.SetIndices(LM);
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
bool FERigidPlanarJoint::Augment(int naug, const FETimeInfo& tp)
{
	if (m_laugon != FECore::AUGLAG_METHOD) return true;

    vec3d ra, rb, qa, qb, c, ksi, Lm;
    vec3d za, zb;
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    double normF0, normF1;
    vec3d Um;
    double normM0, normM1;
    bool bconv = true;
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    double alpha = tp.alphaf;
    
    ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*alpha + zap*(1-alpha);
	eat[0] = m_ea0[0]; RBa.GetRotation().RotateVector(eat[0]);
	eat[1] = m_ea0[1]; RBa.GetRotation().RotateVector(eat[1]);
	eat[2] = m_ea0[2]; RBa.GetRotation().RotateVector(eat[2]);
    eap[0] = m_ea0[0]; RBa.m_qp.RotateVector(eap[0]);
    eap[1] = m_ea0[1]; RBa.m_qp.RotateVector(eap[1]);
    eap[2] = m_ea0[2]; RBa.m_qp.RotateVector(eap[2]);
    ea[0] = eat[0]*alpha + eap[0]*(1-alpha);
    ea[1] = eat[1]*alpha + eap[1]*(1-alpha);
    ea[2] = eat[2]*alpha + eap[2]*(1-alpha);
    
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*alpha + zbp*(1-alpha);
	ebt[0] = m_eb0[0]; RBb.GetRotation().RotateVector(ebt[0]);
	ebt[1] = m_eb0[1]; RBb.GetRotation().RotateVector(ebt[1]);
	ebt[2] = m_eb0[2]; RBb.GetRotation().RotateVector(ebt[2]);
    ebp[0] = m_eb0[0]; RBb.m_qp.RotateVector(ebp[0]);
    ebp[1] = m_eb0[1]; RBb.m_qp.RotateVector(ebp[1]);
    ebp[2] = m_eb0[2]; RBb.m_qp.RotateVector(ebp[2]);
    eb[0] = ebt[0]*alpha + ebp[0]*(1-alpha);
    eb[1] = ebt[1]*alpha + ebp[1]*(1-alpha);
    eb[2] = ebt[2]*alpha + ebp[2]*(1-alpha);
    
    mat3ds P;
    vec3d p;
    if (m_bdy && m_bdz) {
        P = mat3dd(1);
        p = ea[1]*m_dpy + ea[2]*m_dpz;
    }
    else if (m_bdy) {
        P = mat3dd(1) - dyad(ea[2]);
        p = ea[1]*m_dpy;
    }
    else if (m_bdz) {
        P = mat3dd(1) - dyad(ea[1]);
        p = ea[2]*m_dpz;
    }
    else {
        P = dyad(ea[0]);
        p = vec3d(0, 0, 0);
    }
    c = P*(rb + zb - ra - za) - p;
    
    if (m_bqx) {
		quatd q = (alpha*RBb.GetRotation() + (1 - alpha)*RBb.m_qp)*(alpha*RBa.GetRotation() + (1 - alpha)*RBa.m_qp).Inverse();
        quatd a(m_qpx,ea[0]);
        quatd r = a*q.Inverse();
        r.MakeUnit();
        ksi = r.GetVector()*r.GetAngle();
    }
    else
    {
        ksi = (ea[0] ^ eb[0])/2;
    }
    
    normF0 = sqrt(m_L*m_L);
    
    // calculate trial multiplier
    Lm = m_L + c*m_eps;
    
    normF1 = sqrt(Lm*Lm);
    
    normM0 = sqrt(m_U*m_U);
    
    // calculate trial multiplier
    Um = m_U + ksi*m_ups;
    
    normM1 = sqrt(Um*Um);
    
    // check convergence of constraints
    feLog(" rigid connector # %d (planar joint)\n", m_nID+1);
    feLog("                  CURRENT        REQUIRED\n");
    double pctn = 0;
    double gap = c.norm();
    double qap = ksi.norm();
    if (fabs(normF1) > 1e-10) pctn = fabs((normF1 - normF0)/normF1);
    if (m_atol) feLog("    force : %15le %15le\n", pctn, m_atol);
    else        feLog("    force : %15le        ***\n", pctn);
    if (m_gtol) feLog("    gap   : %15le %15le\n", gap, m_gtol);
    else        feLog("    gap   : %15le        ***\n", gap);
    double qctn = 0;
    if (fabs(normM1) > 1e-10) qctn = fabs((normM1 - normM0)/normM1);
    if (m_atol) feLog("    moment: %15le %15le\n", qctn, m_atol);
    else        feLog("    moment: %15le        ***\n", qctn);
    if (m_qtol) feLog("    angle : %15le %15le\n", qap, m_qtol);
    else        feLog("    angle : %15le        ***\n", qap);
    
    if (m_atol && ((pctn >= m_atol) || (qctn >= m_atol))) bconv = false;
    if (m_gtol && (gap >= m_gtol)) bconv = false;
    if (m_qtol && (qap >= m_qtol)) bconv = false;
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    if (!bconv)
    {
        // update multipliers
        m_L = Lm;
        m_U = Um;
    }
    
    // auto-penalty update (works only with gaptol and angtol)
	if (m_bautopen) 
	{
		if (m_gtol && (gap > m_gtol)) {
			m_eps = fmax(gap / m_gtol, 100)*m_eps;
			feLog("    force_penalty :         %15le\n", m_eps);
		}
		if (m_qtol && (qap > m_qtol)) {
			m_ups = fmax(qap / m_qtol, 100)*m_ups;
			feLog("    moment_penalty :        %15le\n", m_ups);
		}
	}

    return bconv;
}

//-----------------------------------------------------------------------------
void FERigidPlanarJoint::Update()
{
	if (m_laugon == FECore::LAGMULT_METHOD) return;

    vec3d ra, rb;
    vec3d za, zb;
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    const FETimeInfo& tp = GetTimeInfo();
    double alpha = tp.alphaf;
    
    ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*alpha + zap*(1-alpha);
	eat[0] = m_ea0[0]; RBa.GetRotation().RotateVector(eat[0]);
	eat[1] = m_ea0[1]; RBa.GetRotation().RotateVector(eat[1]);
	eat[2] = m_ea0[2]; RBa.GetRotation().RotateVector(eat[2]);
    eap[0] = m_ea0[0]; RBa.m_qp.RotateVector(eap[0]);
    eap[1] = m_ea0[1]; RBa.m_qp.RotateVector(eap[1]);
    eap[2] = m_ea0[2]; RBa.m_qp.RotateVector(eap[2]);
    ea[0] = eat[0]*alpha + eap[0]*(1-alpha);
    ea[1] = eat[1]*alpha + eap[1]*(1-alpha);
    ea[2] = eat[2]*alpha + eap[2]*(1-alpha);
    
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*alpha + zbp*(1-alpha);
	ebt[0] = m_eb0[0]; RBb.GetRotation().RotateVector(ebt[0]);
	ebt[1] = m_eb0[1]; RBb.GetRotation().RotateVector(ebt[1]);
	ebt[2] = m_eb0[2]; RBb.GetRotation().RotateVector(ebt[2]);
    ebp[0] = m_eb0[0]; RBb.m_qp.RotateVector(ebp[0]);
    ebp[1] = m_eb0[1]; RBb.m_qp.RotateVector(ebp[1]);
    ebp[2] = m_eb0[2]; RBb.m_qp.RotateVector(ebp[2]);
    eb[0] = ebt[0]*alpha + ebp[0]*(1-alpha);
    eb[1] = ebt[1]*alpha + ebp[1]*(1-alpha);
    eb[2] = ebt[2]*alpha + ebp[2]*(1-alpha);
    
    mat3ds P;
    vec3d p;
    if (m_bdy && m_bdz) {
        P = mat3dd(1);
        p = ea[1]*m_dpy + ea[2]*m_dpz;
    }
    else if (m_bdy) {
        P = mat3dd(1) - dyad(ea[2]);
        p = ea[1]*m_dpy;
    }
    else if (m_bdz) {
        P = mat3dd(1) - dyad(ea[1]);
        p = ea[2]*m_dpz;
    }
    else {
        P = dyad(ea[0]);
        p = vec3d(0, 0, 0);
    }
    vec3d c = P*(rb + zb - ra - za) - p;
    m_F = m_L + c*m_eps;
    
    vec3d ksi;
    if (m_bqx) {
		quatd q = (alpha*RBb.GetRotation() + (1 - alpha)*RBb.m_qp)*(alpha*RBa.GetRotation() + (1 - alpha)*RBa.m_qp).Inverse();
        quatd a(m_qpx,ea[0]);
        quatd r = a*q.Inverse();
        r.MakeUnit();
        ksi = r.GetVector()*r.GetAngle();
    }
    else
    {
        ksi = (ea[0] ^ eb[0])/2;
    }
    m_M = m_U + ksi*m_ups;
    
}

//-----------------------------------------------------------------------------
void FERigidPlanarJoint::Reset()
{
    m_F = vec3d(0,0,0);
    m_L = vec3d(0,0,0);
    m_M = vec3d(0,0,0);
    m_U = vec3d(0,0,0);
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    m_qa0 = m_q0 - RBa.m_r0;
    m_qb0 = m_q0 - RBb.m_r0;
    
    m_ea0[0] = m_e0[0]; m_ea0[1] = m_e0[1]; m_ea0[2] = m_e0[2];
    m_eb0[0] = m_e0[0]; m_eb0[1] = m_e0[1]; m_eb0[2] = m_e0[2];
}

//-----------------------------------------------------------------------------
vec3d FERigidPlanarJoint::RelativeTranslation(const bool global)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;
    
    // body A
    vec3d ra = RBa.m_rt;
    vec3d za = m_qa0; RBa.GetRotation().RotateVector(za);
    
    // body B
    vec3d rb = RBb.m_rt;
    vec3d zb = m_qb0; RBb.GetRotation().RotateVector(zb);

    // relative translation in global coordinate system
    vec3d x = rb + zb - ra - za;
    
    if (global) return x;

    // evaluate local basis for body A
    vec3d ea[3];
    ea[0] = m_ea0[0]; RBa.GetRotation().RotateVector(ea[0]);
    ea[1] = m_ea0[1]; RBa.GetRotation().RotateVector(ea[1]);
    ea[2] = m_ea0[2]; RBa.GetRotation().RotateVector(ea[2]);

    // project relative translation onto local basis
    vec3d y(x*ea[0], x*ea[1], x*ea[2]);
    
    return y;
}

//-----------------------------------------------------------------------------
vec3d FERigidPlanarJoint::RelativeRotation(const bool global)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;
    
    // get relative rotation
    quatd Q = RBb.GetRotation()*RBa.GetRotation().Inverse(); Q.MakeUnit();
    
    // relative rotation vector
    vec3d q = Q.GetRotationVector();
    
    if (global) return q;
    
    // evaluate local basis for body A
    vec3d ea[3];
    ea[0] = m_ea0[0]; RBa.GetRotation().RotateVector(ea[0]);
    ea[1] = m_ea0[1]; RBa.GetRotation().RotateVector(ea[1]);
    ea[2] = m_ea0[2]; RBa.GetRotation().RotateVector(ea[2]);

    // project relative rotation onto local basis
    vec3d y(q*ea[0], q*ea[1], q*ea[2]);
    
    return y;
}

// allocate equations
int FERigidPlanarJoint::InitEquations(int neq)
{
	const int LMeq = 6;
	m_LM.resize(LMeq, -1);
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		for (int i = 0; i < LMeq; ++i) m_LM[i] = neq + i;
		return LMeq;
	}
	else return 0;
}

// Build the matrix profile
void FERigidPlanarJoint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> lm;
	UnpackLM(lm);

	// add it to the pile
	M.build_add(lm);
}

void FERigidPlanarJoint::UnpackLM(vector<int>& lm)
{
	// add the dofs of rigid body A
	lm.resize(18);
	lm[0] = m_rbA->m_LM[0];
	lm[1] = m_rbA->m_LM[1];
	lm[2] = m_rbA->m_LM[2];
	lm[3] = m_rbA->m_LM[3];
	lm[4] = m_rbA->m_LM[4];
	lm[5] = m_rbA->m_LM[5];

	// add the dofs of rigid body B
	lm[6] = m_rbB->m_LM[0];
	lm[7] = m_rbB->m_LM[1];
	lm[8] = m_rbB->m_LM[2];
	lm[9] = m_rbB->m_LM[3];
	lm[10] = m_rbB->m_LM[4];
	lm[11] = m_rbB->m_LM[5];

	// add the LM equations
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		for (int i = 0; i < m_LM.size(); ++i) lm[12 + i] = m_LM[i];
	}
}

void FERigidPlanarJoint::Update(const std::vector<double>& ui)
{
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		m_F.x += ui[m_LM[0]];
		m_F.y += ui[m_LM[1]];
		m_F.z += ui[m_LM[2]];

		m_M.x += ui[m_LM[3]];
		m_M.y += ui[m_LM[4]];
		m_M.z += ui[m_LM[5]];
	}
}
