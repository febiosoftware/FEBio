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
#include "FERigidRevoluteJoint.h"
#include "FERigidBody.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidRevoluteJoint, FERigidConnector);
	ADD_PARAMETER(m_laugon, "laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0LAGMULT\0");
	ADD_PARAMETER(m_atol, "tolerance"     );
    ADD_PARAMETER(m_gtol, "gaptol"        );
    ADD_PARAMETER(m_qtol, "angtol"        );
    ADD_PARAMETER(m_eps , "force_penalty" );
    ADD_PARAMETER(m_ups , "moment_penalty");
    ADD_PARAMETER(m_cps , "force_damping" );
    ADD_PARAMETER(m_rps , "moment_damping");
    ADD_PARAMETER(m_q0  , "joint_origin"  );
    ADD_PARAMETER(m_e0[0], "rotation_axis" );
    ADD_PARAMETER(m_e0[1], "transverse_axis");
    ADD_PARAMETER(m_naugmin, "minaug"        );
    ADD_PARAMETER(m_naugmax, "maxaug"        );
    ADD_PARAMETER(m_bq  , "prescribed_rotation");
    ADD_PARAMETER(m_qp  , "rotation"      );
    ADD_PARAMETER(m_Mp  , "moment"        );
	ADD_PARAMETER(m_bautopen, "auto_penalty");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidRevoluteJoint::FERigidRevoluteJoint(FEModel* pfem) : FERigidConnector(pfem)
{
    m_nID = m_ncount++;

	m_laugon = FECore::AUGLAG_METHOD; // for backward compatibility
	m_atol = 0;
    m_gtol = 0;
    m_qtol = 0;
    m_naugmin = 0;
    m_naugmax = 10;
    m_qp = 0;
    m_Mp = 0;
    m_bq = false;
    m_eps = m_ups = 1.0;
    m_cps = m_rps = 0.0;
    m_e0[0] = vec3d(0,0,1);
    m_e0[1] = vec3d(1,0,0);
	m_bautopen = false;
}

//-----------------------------------------------------------------------------
FERigidRevoluteJoint::~FERigidRevoluteJoint()
{

}

//-----------------------------------------------------------------------------
//! initial position 
vec3d FERigidRevoluteJoint::InitialPosition() const
{
	return m_q0;
}

//-----------------------------------------------------------------------------
//! current position
vec3d FERigidRevoluteJoint::Position() const
{
	FERigidBody& RBa = *m_rbA;
	vec3d qa = m_qa0;
	RBa.GetRotation().RotateVector(qa);
	return RBa.m_rt + qa;
}

//-----------------------------------------------------------------------------
//! current axis
quatd FERigidRevoluteJoint::Orientation() const
{
	quatd Q0(vec3d(0, 0, 1), m_e0[0]);
	FERigidBody& RBa = *m_rbA;
	return RBa.GetRotation()*Q0;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidRevoluteJoint::Init()
{
    if (m_bq && (m_Mp != 0)) {
        feLogError("Rotation and moment cannot be prescribed simultaneously in rigid connector %d (revolute joint)\n", m_nID+1);
        return false;
    }
    
    // initialize joint basis
    m_e0[0].unit();
    m_e0[2] = m_e0[0] ^ m_e0[1]; m_e0[2].unit();
    m_e0[1] = m_e0[2] ^ m_e0[0];
    
    // reset force
    m_F = vec3d(0,0,0); m_L = vec3d(0,0,0); m_Fp = vec3d(0,0,0);
    m_M = vec3d(0,0,0); m_U = vec3d(0,0,0); m_Up = vec3d(0,0,0);
    
	// base class first
	if (FERigidConnector::Init() == false) return false;

    m_qa0 = m_q0 - m_rbA->m_r0;
    m_qb0 = m_q0 - m_rbB->m_r0;
    
    m_ea0[0] = m_e0[0]; m_ea0[1] = m_e0[1]; m_ea0[2] = m_e0[2];
    m_eb0[0] = m_e0[0]; m_eb0[1] = m_e0[1]; m_eb0[2] = m_e0[2];
    
    return true;
}

//-----------------------------------------------------------------------------
void FERigidRevoluteJoint::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
    ar & m_qa0 & m_qb0;
	ar & m_L & m_U;
	ar & m_e0;
    ar & m_ea0;
    ar & m_eb0;
	ar & m_LM;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidRevoluteJoint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
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
		vec3d c = rb + zb - ra - za;
		m_F = m_L + c * m_eps;

		vec3d ksi;
		if (m_bq) {
			quatd q = (alpha * RBb.GetRotation() + (1 - alpha) * RBb.m_qp) * (alpha * RBa.GetRotation() + (1 - alpha) * RBa.m_qp).Inverse();
			quatd a(m_qp, ea[0]);
			quatd r = a * q.Inverse();
			r.MakeUnit();
			ksi = r.GetVector() * r.GetAngle();
		}
		else
			ksi = (ea[0] ^ eb[0]) / 2;
		m_M = m_U + ksi * m_ups + ea[0] * m_Mp;

		// add damping
		if (m_cps > 0) {
			// body A
			vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
			vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
			vec3d va = vat * alpha + vap * (1 - alpha);

			// body b
			vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
			vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
			vec3d vb = vbt * alpha + vbp * (1 - alpha);

			m_F += (vb - va)*m_cps;
		}
		if (m_rps > 0) {
			// body A
			vec3d wa = RBa.m_wt * alpha + RBa.m_wp * (1 - alpha);

			// body b
			vec3d wb = RBb.m_wt * alpha + RBb.m_wp * (1 - alpha);

			mat3ds P = m_bq ? mat3dd(1) : mat3dd(1) - dyad(ea[0]);

			m_M += P * (wb - wa) * m_rps;
		}

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
		vec3d Ma = (za ^ m_F) + (ea[0] ^ m_M);
		vec3d Mb = (zb ^ m_F) + (eb[0] ^ m_M);

		fa[0] = -m_F.x;
		fa[1] = -m_F.y;
		fa[2] = -m_F.z;
		fa[3] = -Ma.x;
		fa[4] = -Ma.y;
		fa[5] = -Ma.z;

		fb[0] = m_F.x;
		fb[1] = m_F.y;
		fb[2] = m_F.z;
		fb[3] = Mb.x;
		fb[4] = Mb.y;
		fb[5] = Mb.z;

		for (int i = 0; i < 6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
		for (int i = 0; i < 6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];

		// translational constraint
		vec3d c = rb + zb - ra - za;
		R[m_LM[0]] += c.x;
		R[m_LM[1]] += c.y;
		R[m_LM[2]] += c.z;

		// rotational constraint
		vec3d ksi1 = eb[0] - ea[0];
		R[m_LM[3]] += ksi1.x;
		R[m_LM[4]] += ksi1.y;
		R[m_LM[5]] += ksi1.z;
	}
    
    RBa.m_Fr -= vec3d(fa[0],fa[1],fa[2]);
    RBa.m_Mr -= vec3d(fa[3],fa[4],fa[5]);
    RBb.m_Fr -= vec3d(fb[0],fb[1],fb[2]);
    RBb.m_Mr -= vec3d(fb[3],fb[4],fb[5]);
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidRevoluteJoint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	double alpha = tp.alphaf;
    double beta  = tp.beta;
    double gamma = tp.gamma;
    
    // get time increment
    double dt = tp.timeIncrement;
    
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
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
    mat3d zbhat; zbhat.skew(zb);
    mat3d zbthat; zbthat.skew(zbt);

	mat3d eahat[3], ebhat[3], eathat[3], ebthat[3];
	for (int j = 0; j < 3; ++j) {
		eahat[j] = skew(ea[j]);
		ebhat[j] = skew(eb[j]);
		eathat[j] = skew(eat[j]);
		ebthat[j] = skew(ebt[j]);
	}

	if (m_laugon != FECore::LAGMULT_METHOD)
	{
		vec3d c = rb + zb - ra - za;
		m_F = m_L + c * m_eps;
		mat3dd I(1);

		vec3d ksi;
		quatd q, a, r;
		if (m_bq) {
			q = (alpha * RBb.GetRotation() + (1 - alpha) * RBb.m_qp) * (alpha * RBa.GetRotation() + (1 - alpha) * RBa.m_qp).Inverse();
			a = quatd(m_qp, ea[0]);
			r = a * q.Inverse();
			r.MakeUnit();
			ksi = r.GetVector() * r.GetAngle();
		}
		else
			ksi = (ea[0] ^ eb[0]) / 2;
		m_M = m_U + ksi * m_ups + ea[0] * m_Mp;

		mat3d K, Wba, Wab;
		Wba = (ebhat[0] * eathat[0]) / 2;
		Wab = (eahat[0] * ebthat[0]) / 2;
		if (m_bq) {
			quatd qa = RBa.GetRotation() * (alpha * RBa.GetRotation() + (1 - alpha) * RBa.m_qp).Inverse();
			quatd qb = RBb.GetRotation() * (alpha * RBb.GetRotation() + (1 - alpha) * RBb.m_qp).Inverse();
			qa.MakeUnit();
			qb.MakeUnit();
			mat3d Qa = qa.RotationMatrix();
			mat3d Qb = qb.RotationMatrix();
			mat3d A = a.RotationMatrix();
			mat3d R = r.RotationMatrix();
			mat3dd I(1);
			Wba = A *(I * Qa.trace() - Qa) / 2;
			Wab = R *(I * Qb.trace() - Qb) / 2;
		}

		// add damping
		mat3ds A;
		mat3d Ba, Bb, Ca, Cb, W;
		A.zero(); Ba.zero(); Bb.zero(); Ca.zero(); Cb.zero(); W.zero();
		if ((m_cps > 0) || (m_rps > 0)) {
			// body A
			vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
			vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
			vec3d va = vat * alpha + vap * (1 - alpha);
			quatd qai = RBa.GetRotation() * RBa.m_qp.Inverse(); qai.MakeUnit();
			vec3d cai = qai.GetVector() * (2 * tan(qai.GetAngle() / 2));
			mat3d Ta = I + skew(cai) / 2 + dyad(cai) / 4;
			vec3d wa = RBa.m_wt * alpha + RBa.m_wp * (1 - alpha);

			// body b
			vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
			vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
			vec3d vb = vbt * alpha + vbp * (1 - alpha);
			quatd qbi = RBb.GetRotation() * RBb.m_qp.Inverse(); qbi.MakeUnit();
			vec3d cbi = qbi.GetVector() * (2 * tan(qbi.GetAngle() / 2));
			mat3d Tb = I + skew(cbi) / 2 + dyad(cbi) / 4;
			vec3d wb = RBb.m_wt * alpha + RBb.m_wp * (1 - alpha);

			m_F += (vb - va)*m_cps;

			vec3d w = wb - wa;

			// angular damping along all directions if rotation is prescribed
			mat3ds P = m_bq ? I : I - dyad(ea[0]);

			m_M += P*w*m_rps;

			A = I*(gamma/beta/dt);
			Ba = zathat * Ta.transpose() * (gamma / beta / dt) + skew(RBa.m_wt) * zathat;
			Bb = zbthat * Tb.transpose() * (gamma / beta / dt) + skew(RBb.m_wt) * zbthat;
			Ca = P * Ta.transpose() * (gamma / beta / dt);
			Cb = P * Tb.transpose() * (gamma / beta / dt);
			W = (mat3dd(ea[0] * w) + (ea[0] & w)) * eathat[0];
		}

		mat3d Fhat; Fhat.skew(m_F);
		mat3d Mhat; Mhat.skew(m_M);

		FEElementMatrix ke(12, 12);
		ke.zero();

		// row 1
		ke.set(0, 0, I * (m_eps) + A * (m_cps));
		ke.set(0, 3, zathat * (-m_eps) + Ba * (-m_cps));
		ke.set(0, 6, I * (-m_eps) + A * (-m_cps));
		ke.set(0, 9, zbthat * (m_eps) + Bb * (m_cps));

		// row 2
		ke.set(3, 0, zahat* (m_eps) + zahat * A * (m_cps));
		ke.set(3, 3, (zahat * zathat * m_eps + Fhat * zathat + Wba * m_ups) * (-1.0));
		ke.set(3, 6, zahat * (-m_eps) + zahat * A * (-m_cps));
		ke.set(3, 9, (zahat* zbthat* m_eps + Wab * m_ups));

		if (m_cps)
		{
			ke.add(3, 3, eathat[0] * (m_Mp)+(zahat * Ba) * (-m_cps));
			ke.add(3, 9, (zahat * Bb) * (m_cps));
		}
		if (m_rps)
		{
			ke.add(3, 3, (W - Ca) * (m_rps));
			ke.add(3, 9, Cb * (m_rps));
		}

		// row 3
		ke.set(6, 0, I* (-m_eps) + A * (-m_cps));
		ke.set(6, 3, zathat * (m_eps) + Ba * (m_cps));
		ke.set(6, 6, I * (m_eps) + A * (m_cps));
		ke.set(6, 9, zbthat * (-m_eps) + Bb * (-m_cps));

		// row 4
		ke.set(9, 0, zbhat* (-m_eps) + zbhat * A * (-m_cps));
		ke.set(9, 3, zbhat* zathat* m_eps + Wba * m_ups);
		ke.set(9, 6, zbhat* (m_eps) + zbhat * A * (m_cps));
		ke.set(9, 9, (zbhat * zbthat * m_eps - Fhat * zbthat + Wab * m_ups) * (-1.0));

		if (m_cps)
		{
			ke.add(9, 0, -eathat[0] * (m_Mp)+(zbhat * Ba) * (m_cps));
			ke.add(9, 9, (zbhat * Bb) * (-m_cps));
		}
		if (m_rps)
		{
			ke.add(9, 0, -(W - Ca) * (m_rps));
			ke.add(9, 9, -Cb * (m_rps));
		}

		ke *= alpha;

		vector<int> LM(12);
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
		FEElementMatrix ke; ke.resize(18, 18);
		ke.zero();

		mat3dd I(1.0);

		// translation constraint
		mat3da Fhat(m_F);
		ke.sub( 3,  3, -Fhat * zahat);
		ke.sub( 9,  9, Fhat* zbhat);
		ke.sub( 0, 12, -I);
		ke.sub( 3, 12, -zahat);
		ke.sub( 6, 12, I);
		ke.sub( 9, 12, zbhat);
		ke.sub(12,  0, -I);
		ke.sub(12,  3, zahat);
		ke.sub(12,  6, I);
		ke.sub(12,  9, -zbhat);

		// rotational constraint 1
		mat3da Mhat(m_M);
		ke.sub(3, 3, -Mhat*eahat[0]);
		ke.sub(9, 9,  Mhat*ebhat[0]);
		ke.sub(3, 15, -eahat[0]);
		ke.sub(9, 15,  ebhat[0]);
		ke.sub(15, 3,  eahat[0]);
		ke.sub(15, 9, -ebhat[0]);

		vector<int> LM;
		UnpackLM(LM);
		ke.SetIndices(LM);
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
bool FERigidRevoluteJoint::Augment(int naug, const FETimeInfo& tp)
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
    
    c = rb + zb - ra - za;
    
    normF0 = sqrt(m_L*m_L);
    
    // calculate trial multiplier
    Lm = m_L + c*m_eps;
    
    normF1 = sqrt(Lm*Lm);
    
    if (m_bq) {
		quatd q = (alpha*RBb.GetRotation() + (1 - alpha)*RBb.m_qp)*(alpha*RBa.GetRotation() + (1 - alpha)*RBa.m_qp).Inverse();
        quatd a(m_qp,ea[0]);
        quatd r = a*q.Inverse();
        r.MakeUnit();
        ksi = r.GetVector()*r.GetAngle();
    }
    else
        ksi = (ea[0] ^ eb[0])/2;
    
    normM0 = sqrt(m_U*m_U);
    
    // calculate trial multiplier
    Um = m_U + ksi*m_ups;
    
    normM1 = sqrt(Um*Um);
    
    // check convergence of constraints
    feLog(" rigid connector # %d (revolute joint)\n", m_nID+1);
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
void FERigidRevoluteJoint::Update()
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
    
    vec3d c = rb + zb - ra - za;
    m_F = m_L + c*m_eps;
    
    vec3d ksi;
    if (m_bq) {
		quatd q = (alpha*RBb.GetRotation() + (1 - alpha)*RBb.m_qp)*(alpha*RBa.GetRotation() + (1 - alpha)*RBa.m_qp).Inverse();
        quatd a(m_qp,ea[0]);
        quatd r = a*q.Inverse();
        r.MakeUnit();
        ksi = r.GetVector()*r.GetAngle();
    }
    else
        ksi = (ea[0] ^ eb[0])/2;
    m_M = m_U + ksi*m_ups + ea[0]*m_Mp;
    
    // add damping
    if (m_cps > 0) {
        // body A
        vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
        vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
        vec3d va = vat*alpha + vap*(1-alpha);
        
        // body b
        vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
        vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
        vec3d vb = vbt*alpha + vbp*(1-alpha);
        
        m_F += (vb - va)*m_cps;
    }
    if (m_rps > 0) {
        // body A
        vec3d wa = RBa.m_wt*alpha + RBa.m_wp*(1-alpha);
        
        // body b
        vec3d wb = RBb.m_wt*alpha + RBb.m_wp*(1-alpha);
        
        mat3ds P = m_bq ? mat3dd(1) : mat3dd(1) - dyad(ea[0]);
        
        m_M += P*(wb - wa)*m_rps;
    }
    
}

//-----------------------------------------------------------------------------
void FERigidRevoluteJoint::Reset()
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
vec3d FERigidRevoluteJoint::RelativeTranslation(const bool global)
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
vec3d FERigidRevoluteJoint::RelativeRotation(const bool global)
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
int FERigidRevoluteJoint::InitEquations(int neq)
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
void FERigidRevoluteJoint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> lm;
	UnpackLM(lm);

	// add it to the pile
	M.build_add(lm);
}

void FERigidRevoluteJoint::UnpackLM(vector<int>& lm)
{
	// add the dofs of rigid body A
	lm.resize(21);
	lm[0] = m_rbA->m_LM[0];
	lm[1] = m_rbA->m_LM[1];
	lm[2] = m_rbA->m_LM[2];
	lm[3] = m_rbA->m_LM[3];
	lm[4] = m_rbA->m_LM[4];
	lm[5] = m_rbA->m_LM[5];

	// add the dofs of rigid body B
	lm[ 6] = m_rbB->m_LM[0];
	lm[ 7] = m_rbB->m_LM[1];
	lm[ 8] = m_rbB->m_LM[2];
	lm[ 9] = m_rbB->m_LM[3];
	lm[10] = m_rbB->m_LM[4];
	lm[11] = m_rbB->m_LM[5];

	// add the LM equations
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		for (int i = 0; i < m_LM.size(); ++i) lm[12 + i] = m_LM[i];
	}
}

void FERigidRevoluteJoint::PrepStep()
{
	m_Fp = m_F;
	m_Up = m_M;
}

void FERigidRevoluteJoint::Update(const std::vector<double>& Ui, const std::vector<double>& ui)
{
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		m_F.x = m_Fp.x + Ui[m_LM[0]] + ui[m_LM[0]];
		m_F.y = m_Fp.y + Ui[m_LM[1]] + ui[m_LM[1]];
		m_F.z = m_Fp.z + Ui[m_LM[2]] + ui[m_LM[2]];

		m_M.x = m_Up.x + Ui[m_LM[3]] + ui[m_LM[3]];
		m_M.y = m_Up.y + Ui[m_LM[4]] + ui[m_LM[4]];
		m_M.z = m_Up.z + Ui[m_LM[5]] + ui[m_LM[5]];
	}
}

void FERigidRevoluteJoint::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui)
{
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		Ui[m_LM[0]] += ui[m_LM[0]];
		Ui[m_LM[1]] += ui[m_LM[1]];
		Ui[m_LM[2]] += ui[m_LM[2]];
		Ui[m_LM[3]] += ui[m_LM[3]];
		Ui[m_LM[4]] += ui[m_LM[4]];
		Ui[m_LM[5]] += ui[m_LM[5]];
	}
}
