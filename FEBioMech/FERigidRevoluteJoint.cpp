/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidRevoluteJoint, FERigidConnector);
	ADD_PARAMETER(m_blaugon, "laugon");
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
	m_blaugon = true; // on by default for backward compatibility
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
    
    FEModel& fem = *GetFEModel();
    
    // initialize joint basis
    m_e0[0].unit();
    m_e0[2] = m_e0[0] ^ m_e0[1]; m_e0[2].unit();
    m_e0[1] = m_e0[2] ^ m_e0[0];
    
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
void FERigidRevoluteJoint::Serialize(DumpStream& ar)
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
    
    fa[0] = m_F.x;
    fa[1] = m_F.y;
    fa[2] = m_F.z;
    
    fa[3] = za.y*m_F.z-za.z*m_F.y + m_M.x;
    fa[4] = za.z*m_F.x-za.x*m_F.z + m_M.y;
    fa[5] = za.x*m_F.y-za.y*m_F.x + m_M.z;
    
    fb[0] = -m_F.x;
    fb[1] = -m_F.y;
    fb[2] = -m_F.z;
    
    fb[3] = -zb.y*m_F.z+zb.z*m_F.y - m_M.x;
    fb[4] = -zb.z*m_F.x+zb.x*m_F.z - m_M.y;
    fb[5] = -zb.x*m_F.y+zb.y*m_F.x - m_M.z;
    
    for (int i=0; i<6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
    for (int i=0; i<6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];
    
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
    
    vector<int> LM(12);
    FEElementMatrix ke(12,12);
    ke.zero();
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

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
    mat3d zahat; zahat.skew(za);
    mat3d zathat; zathat.skew(zat);
    
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
    mat3d zbhat; zbhat.skew(zb);
    mat3d zbthat; zbthat.skew(zbt);
    
    vec3d c = rb + zb - ra - za;
    m_F = m_L + c*m_eps;
    mat3dd I(1);
    
    vec3d ksi;
    quatd q, a, r;
    if (m_bq) {
		q = (alpha*RBb.GetRotation() + (1 - alpha)*RBb.m_qp)*(alpha*RBa.GetRotation() + (1 - alpha)*RBa.m_qp).Inverse();
        a = quatd(m_qp,ea[0]);
        r = a*q.Inverse();
        r.MakeUnit();
        ksi = r.GetVector()*r.GetAngle();
    }
    else
        ksi = (ea[0] ^ eb[0])/2;
    m_M = m_U + ksi*m_ups + ea[0]*m_Mp;
    
    mat3d eahat[3], ebhat[3], eathat[3], ebthat[3];
    for (int j=0; j<3; ++j) {
        eahat[j] = skew(ea[j]);
        ebhat[j] = skew(eb[j]);
        eathat[j] = skew(eat[j]);
        ebthat[j] = skew(ebt[j]);
    }
    mat3d K, Wba, Wab;
    Wba = (ebhat[0]*eathat[0])/2;
    Wab = (eahat[0]*ebthat[0])/2;
    if (m_bq) {
		quatd qa = RBa.GetRotation()*(alpha*RBa.GetRotation() + (1 - alpha)*RBa.m_qp).Inverse();
		quatd qb = RBb.GetRotation()*(alpha*RBb.GetRotation() + (1 - alpha)*RBb.m_qp).Inverse();
        qa.MakeUnit();
        qb.MakeUnit();
        mat3d Qa = qa.RotationMatrix();
        mat3d Qb = qb.RotationMatrix();
        mat3d A = a.RotationMatrix();
        mat3d R = r.RotationMatrix();
        mat3dd I(1);
        Wba = A*(I*Qa.trace()-Qa)/2;
        Wab = R*(I*Qb.trace()-Qb)/2;
    }
    
    // add damping
    mat3ds A;
    mat3d Ba, Bb, Ca, Cb, W;
    A.zero(); Ba.zero(); Bb.zero(); Ca.zero(); Cb.zero(); W.zero();
    if ((m_cps > 0) || (m_rps > 0)) {
        // body A
        vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
        vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
        vec3d va = vat*alpha + vap*(1-alpha);
		quatd qai = RBa.GetRotation()*RBa.m_qp.Inverse(); qai.MakeUnit();
        vec3d cai = qai.GetVector()*(2*tan(qai.GetAngle()/2));
        mat3d Ta = I + skew(cai)/2 + dyad(cai)/4;
        vec3d wa = RBa.m_wt*alpha + RBa.m_wp*(1-alpha);
        
        // body b
        vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
        vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
        vec3d vb = vbt*alpha + vbp*(1-alpha);
		quatd qbi = RBb.GetRotation()*RBb.m_qp.Inverse(); qbi.MakeUnit();
        vec3d cbi = qbi.GetVector()*(2*tan(qbi.GetAngle()/2));
        mat3d Tb = I + skew(cbi)/2 + dyad(cbi)/4;
        vec3d wb = RBb.m_wt*alpha + RBb.m_wp*(1-alpha);
        
        m_F += (vb - va)*m_cps;
        
        vec3d w = wb - wa;
        
        // angular damping along all directions if rotation is prescribed
        mat3ds P = m_bq ? I : I - dyad(ea[0]);
        
        m_M += P*w*m_rps;

        A = I*(gamma/beta/dt);
        Ba = zathat*Ta.transpose()*(gamma/beta/dt) + skew(RBa.m_wt)*zathat;
        Bb = zbthat*Tb.transpose()*(gamma/beta/dt) + skew(RBb.m_wt)*zbthat;
        Ca = P*Ta.transpose()*(gamma/beta/dt);
        Cb = P*Tb.transpose()*(gamma/beta/dt);
        W = (mat3dd(ea[0]*w) + (ea[0] & w))*eathat[0];
    }
    
    // (1,1)
    K = I*(alpha*m_eps) + A*(alpha*m_cps);
    ke[0][0] = K[0][0]; ke[0][1] = K[0][1]; ke[0][2] = K[0][2];
    ke[1][0] = K[1][0]; ke[1][1] = K[1][1]; ke[1][2] = K[1][2];
    ke[2][0] = K[2][0]; ke[2][1] = K[2][1]; ke[2][2] = K[2][2];
    
    // (1,2)
    K = zathat*(-m_eps*alpha) + Ba*(-alpha*m_cps);
    ke[0][3] = K[0][0]; ke[0][4] = K[0][1]; ke[0][5] = K[0][2];
    ke[1][3] = K[1][0]; ke[1][4] = K[1][1]; ke[1][5] = K[1][2];
    ke[2][3] = K[2][0]; ke[2][4] = K[2][1]; ke[2][5] = K[2][2];
    
    // (1,3)
    K = I*(-alpha*m_eps) + A*(-alpha*m_cps);
    ke[0][6] = K[0][0]; ke[0][7] = K[0][1]; ke[0][8] = K[0][2];
    ke[1][6] = K[1][0]; ke[1][7] = K[1][1]; ke[1][8] = K[1][2];
    ke[2][6] = K[2][0]; ke[2][7] = K[2][1]; ke[2][8] = K[2][2];
    
    // (1,4)
    K = zbthat*(alpha*m_eps) + Bb*(alpha*m_cps);
    ke[0][9] = K[0][0]; ke[0][10] = K[0][1]; ke[0][11] = K[0][2];
    ke[1][9] = K[1][0]; ke[1][10] = K[1][1]; ke[1][11] = K[1][2];
    ke[2][9] = K[2][0]; ke[2][10] = K[2][1]; ke[2][11] = K[2][2];
    
    // (2,1)
    K = zahat*(alpha*m_eps) + zahat*A*(alpha*m_cps);
    ke[3][0] = K[0][0]; ke[3][1] = K[0][1]; ke[3][2] = K[0][2];
    ke[4][0] = K[1][0]; ke[4][1] = K[1][1]; ke[4][2] = K[1][2];
    ke[5][0] = K[2][0]; ke[5][1] = K[2][1]; ke[5][2] = K[2][2];
    
    // (2,2)
    K = (zahat*zathat*m_eps + Wba*m_ups)*(-alpha)
    + eathat[0]*(m_Mp*alpha) + (zahat*Ba)*(-alpha*m_cps)
    + (W - Ca)*(alpha*m_rps);
    ke[3][3] = K[0][0]; ke[3][4] = K[0][1]; ke[3][5] = K[0][2];
    ke[4][3] = K[1][0]; ke[4][4] = K[1][1]; ke[4][5] = K[1][2];
    ke[5][3] = K[2][0]; ke[5][4] = K[2][1]; ke[5][5] = K[2][2];
    
    // (2,3)
    K = zahat*(-alpha*m_eps) + zahat*A*(-alpha*m_cps);
    ke[3][6] = K[0][0]; ke[3][7] = K[0][1]; ke[3][8] = K[0][2];
    ke[4][6] = K[1][0]; ke[4][7] = K[1][1]; ke[4][8] = K[1][2];
    ke[5][6] = K[2][0]; ke[5][7] = K[2][1]; ke[5][8] = K[2][2];
    
    // (2,4)
    K = (zahat*zbthat*m_eps + Wab*m_ups)*alpha + (zahat*Bb)*(alpha*m_cps)
    + Cb*(alpha*m_rps);
    ke[3][9] = K[0][0]; ke[3][10] = K[0][1]; ke[3][11] = K[0][2];
    ke[4][9] = K[1][0]; ke[4][10] = K[1][1]; ke[4][11] = K[1][2];
    ke[5][9] = K[2][0]; ke[5][10] = K[2][1]; ke[5][11] = K[2][2];
    
    
    // (3,1)
    K = I*(-alpha*m_eps) + A*(-alpha*m_cps);
    ke[6][0] = K[0][0]; ke[6][1] = K[0][1]; ke[6][2] = K[0][2];
    ke[7][0] = K[1][0]; ke[7][1] = K[1][1]; ke[7][2] = K[1][2];
    ke[8][0] = K[2][0]; ke[8][1] = K[2][1]; ke[8][2] = K[2][2];
    
    // (3,2)
    K = zathat*(m_eps*alpha) + Ba*(alpha*m_cps);
    ke[6][3] = K[0][0]; ke[6][4] = K[0][1]; ke[6][5] = K[0][2];
    ke[7][3] = K[1][0]; ke[7][4] = K[1][1]; ke[7][5] = K[1][2];
    ke[8][3] = K[2][0]; ke[8][4] = K[2][1]; ke[8][5] = K[2][2];
    
    // (3,3)
    K = I*(alpha*m_eps) + A*(alpha*m_cps);
    ke[6][6] = K[0][0]; ke[6][7] = K[0][1]; ke[6][8] = K[0][2];
    ke[7][6] = K[1][0]; ke[7][7] = K[1][1]; ke[7][8] = K[1][2];
    ke[8][6] = K[2][0]; ke[8][7] = K[2][1]; ke[8][8] = K[2][2];
    
    // (3,4)
    K = zbthat*(-alpha*m_eps) + Bb*(-alpha*m_cps);
    ke[6][9] = K[0][0]; ke[6][10] = K[0][1]; ke[6][11] = K[0][2];
    ke[7][9] = K[1][0]; ke[7][10] = K[1][1]; ke[7][11] = K[1][2];
    ke[8][9] = K[2][0]; ke[8][10] = K[2][1]; ke[8][11] = K[2][2];
    

    // (4,1)
    K = zbhat*(-alpha*m_eps) + zbhat*A*(-alpha*m_cps);
    ke[9 ][0] = K[0][0]; ke[ 9][1] = K[0][1]; ke[ 9][2] = K[0][2];
    ke[10][0] = K[1][0]; ke[10][1] = K[1][1]; ke[10][2] = K[1][2];
    ke[11][0] = K[2][0]; ke[11][1] = K[2][1]; ke[11][2] = K[2][2];
    
    // (4,2)
    K = (zbhat*zathat*m_eps + Wba*m_ups)*alpha
    - eathat[0]*(m_Mp*alpha) + (zbhat*Ba)*(alpha*m_cps)
    - (W - Ca)*(alpha*m_rps);
    ke[9 ][3] = K[0][0]; ke[ 9][4] = K[0][1]; ke[ 9][5] = K[0][2];
    ke[10][3] = K[1][0]; ke[10][4] = K[1][1]; ke[10][5] = K[1][2];
    ke[11][3] = K[2][0]; ke[11][4] = K[2][1]; ke[11][5] = K[2][2];
    
    // (4,3)
    K = zbhat*(alpha*m_eps) + zbhat*A*(alpha*m_cps);
    ke[9 ][6] = K[0][0]; ke[ 9][7] = K[0][1]; ke[ 9][8] = K[0][2];
    ke[10][6] = K[1][0]; ke[10][7] = K[1][1]; ke[10][8] = K[1][2];
    ke[11][6] = K[2][0]; ke[11][7] = K[2][1]; ke[11][8] = K[2][2];
    
    // (4,4)
    K = (zbhat*zbthat*m_eps + Wab*m_ups)*(-alpha) + (zbhat*Bb)*(-alpha*m_cps)
    - Cb*(alpha*m_rps);
    ke[9 ][9] = K[0][0]; ke[ 9][10] = K[0][1]; ke[ 9][11] = K[0][2];
    ke[10][9] = K[1][0]; ke[10][10] = K[1][1]; ke[10][11] = K[1][2];
    ke[11][9] = K[2][0]; ke[11][10] = K[2][1]; ke[11][11] = K[2][2];
    
    for (int j=0; j<6; ++j)
    {
        LM[j  ] = RBa.m_LM[j];
        LM[j+6] = RBb.m_LM[j];
    }
    
	ke.SetIndices(LM);
	LS.Assemble(ke);
}

//-----------------------------------------------------------------------------
bool FERigidRevoluteJoint::Augment(int naug, const FETimeInfo& tp)
{
	if (m_blaugon == false) return true;

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
    vec3d ra, rb;
    vec3d za, zb;
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	FETimeInfo& tp = GetFEModel()->GetTime();
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
