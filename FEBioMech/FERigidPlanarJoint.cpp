//
//  FERigidPlanarJoint.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/28/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FERigidPlanarJoint.h"
#include "FECore/FERigidSystem.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidPlanarJoint, FERigidConnector);
ADD_PARAMETER(m_atol, FE_PARAM_DOUBLE, "tolerance"     );
ADD_PARAMETER(m_gtol, FE_PARAM_DOUBLE, "gaptol"        );
ADD_PARAMETER(m_qtol, FE_PARAM_DOUBLE, "angtol"        );
ADD_PARAMETER(m_eps , FE_PARAM_DOUBLE, "force_penalty" );
ADD_PARAMETER(m_ups , FE_PARAM_DOUBLE, "moment_penalty");
ADD_PARAMETER(m_q0  , FE_PARAM_VEC3D , "joint_origin"  );
ADD_PARAMETER(m_e0[0], FE_PARAM_VEC3D, "rotation_axis" );
ADD_PARAMETER(m_e0[1], FE_PARAM_VEC3D, "translation_axis_1");
ADD_PARAMETER(m_naugmin,FE_PARAM_INT , "minaug"        );
ADD_PARAMETER(m_naugmax,FE_PARAM_INT , "maxaug"        );
ADD_PARAMETER(m_bqx , FE_PARAM_BOOL  , "prescribed_rotation");
ADD_PARAMETER(m_qpx , FE_PARAM_DOUBLE, "rotation"      );
ADD_PARAMETER(m_bdy , FE_PARAM_BOOL  , "prescribed_translation_1");
ADD_PARAMETER(m_dpy , FE_PARAM_DOUBLE, "translation_1" );
ADD_PARAMETER(m_bdz , FE_PARAM_BOOL  , "prescribed_translation_2");
ADD_PARAMETER(m_dpz , FE_PARAM_DOUBLE, "translation_2" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidPlanarJoint::FERigidPlanarJoint(FEModel* pfem) : FERigidConnector(pfem)
{
    m_nID = m_ncount++;
    m_binit = false;
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
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidPlanarJoint::Init()
{
    if (m_binit) return true;
    
    // initialize joint basis
    m_e0[0].unit();
    m_e0[2] = m_e0[0] ^ m_e0[1]; m_e0[2].unit();
    m_e0[1] = m_e0[2] ^ m_e0[0]; m_e0[1].unit();
    
    // reset force
    m_F = vec3d(0,0,0); m_L = vec3d(0,0,0);
    m_M = vec3d(0,0,0); m_U = vec3d(0,0,0);
    
    FEModel& fem = *GetFEModel();
    
    // When the rigid joint is read in, the ID's correspond to the rigid materials.
    // Now we want to make the ID's refer to the rigid body ID's
    
    FEMaterial* pm = fem.GetMaterial(m_nRBa-1);
    if (pm->IsRigid() == false)
    {
        felog.printbox("FATAL ERROR", "Rigid connector %d (planar joint) does not connect two rigid bodies\n", m_nID+1);
        return false;
    }
    m_nRBa = pm->GetRigidBodyID();
    
    pm = fem.GetMaterial(m_nRBb-1);
    if (pm->IsRigid() == false)
    {
        felog.printbox("FATAL ERROR", "Rigid connector %d (planar joint) does not connect two rigid bodies\n", m_nID+1);
        return false;
    }
    m_nRBb = pm->GetRigidBodyID();
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
    m_qa0 = m_q0 - RBa.m_r0;
    m_qb0 = m_q0 - RBb.m_r0;
    
    m_ea0[0] = m_e0[0]; m_ea0[1] = m_e0[1]; m_ea0[2] = m_e0[2];
    m_eb0[0] = m_e0[0]; m_eb0[1] = m_e0[1]; m_eb0[2] = m_e0[2];
    
    m_binit = true;
    
    return true;
}

//-----------------------------------------------------------------------------
void FERigidPlanarJoint::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
    if (ar.IsSaving())
    {
		ar << m_binit;
        ar << m_qa0 << m_qb0;
        ar << m_L << m_U;
        ar << m_e0[0] << m_e0[1] << m_e0[2];
        ar << m_ea0[0] << m_ea0[1] << m_ea0[2];
        ar << m_eb0[0] << m_eb0[1] << m_eb0[2];
    }
    else
    {
		ar >> m_binit;
        ar >> m_qa0 >> m_qb0;
        ar >> m_L >> m_U;
        ar >> m_e0[0] >> m_e0[1] >> m_e0[2];
        ar >> m_ea0[0] >> m_ea0[1] >> m_ea0[2];
        ar >> m_eb0[0] >> m_eb0[1] >> m_eb0[2];
    }
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidPlanarJoint::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<double> fa(6);
    vector<double> fb(6);
    
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
    double alpha = tp.alpha;
    
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
    
    RBa.m_Fr += vec3d(fa[0],fa[1],fa[2]);
    RBa.m_Mr += vec3d(fa[3],fa[4],fa[5]);
    RBb.m_Fr += vec3d(fb[0],fb[1],fb[2]);
    RBb.m_Mr += vec3d(fb[3],fb[4],fb[5]);
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidPlanarJoint::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    double alpha = tp.alpha;
    
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
    int j;
    
    vector<int> LM(12);
    matrix ke(12,12);
    ke.zero();
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
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
    
    mat3d eahat[3], ebhat[3], eathat[3], ebthat[3];
    for (j=0; j<3; ++j) {
        eahat[j] = skew(ea[j]);
        ebhat[j] = skew(eb[j]);
        eathat[j] = skew(eat[j]);
        ebthat[j] = skew(ebt[j]);
    }
    
    mat3ds P;
    vec3d p;
    mat3d Q, Wba, Wab;
    vec3d d = rb + zb - ra - za;
    if (m_bdy && m_bdz) {
        P = mat3dd(1);
        p = ea[1]*m_dpy + ea[2]*m_dpz;
        Q = mat3dd(0);
    }
    else if (m_bdy) {
        P = mat3dd(1) - dyad(ea[2]);
        p = ea[1]*m_dpy;
        Q = ((ea[2] & d) + mat3dd(1)*(ea[2]*d))*eathat[2];
    }
    else if (m_bdz) {
        P = mat3dd(1) - dyad(ea[1]);
        p = ea[2]*m_dpz;
        Q = ((ea[1] & d) + mat3dd(1)*(ea[1]*d))*eathat[1];
    }
    else {
        P = dyad(ea[0]);
        p = vec3d(0, 0, 0);
        Q = ((ea[0] & d) + mat3dd(1)*(ea[0]*d))*eathat[0]*(-1);
    }
    vec3d c = P*d - p;
    m_F = m_L + c*m_eps;
    
    vec3d ksi;
    if (m_bqx) {
		quatd q = (alpha*RBb.GetRotation() + (1 - alpha)*RBb.m_qp)*(alpha*RBa.GetRotation() + (1 - alpha)*RBa.m_qp).Inverse();
        quatd a(m_qpx,ea[0]);
        quatd r = a*q.Inverse();
        r.MakeUnit();
        ksi = r.GetVector()*r.GetAngle();
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
    else
    {
        ksi = (ea[0] ^ eb[0])/2;
        Wba = (ebhat[0]*eathat[0])/2;
        Wab = (eahat[0]*ebthat[0])/2;
    }
    m_M = m_U + ksi*m_ups;
    
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
    K = (zahat*(P*zathat+Q)*m_eps + Wba*m_ups)*(-alpha);
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
    K = zbhat*P*(alpha*m_eps);
    ke[9 ][6] = K[0][0]; ke[ 9][7] = K[0][1]; ke[ 9][8] = K[0][2];
    ke[10][6] = K[1][0]; ke[10][7] = K[1][1]; ke[10][8] = K[1][2];
    ke[11][6] = K[2][0]; ke[11][7] = K[2][1]; ke[11][8] = K[2][2];
    
    // (4,4)
    K = (zbhat*P*zbthat*m_eps + Wab*m_ups)*(-alpha);
    ke[9 ][9] = K[0][0]; ke[ 9][10] = K[0][1]; ke[ 9][11] = K[0][2];
    ke[10][9] = K[1][0]; ke[10][10] = K[1][1]; ke[10][11] = K[1][2];
    ke[11][9] = K[2][0]; ke[11][10] = K[2][1]; ke[11][11] = K[2][2];
    
    for (j=0; j<6; ++j)
    {
        LM[j  ] = RBa.m_LM[j];
        LM[j+6] = RBb.m_LM[j];
    }
    
    psolver->AssembleStiffness(LM, ke);
}

//-----------------------------------------------------------------------------
bool FERigidPlanarJoint::Augment(int naug, const FETimeInfo& tp)
{
    vec3d ra, rb, qa, qb, c, ksi, Lm;
    vec3d za, zb;
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    double normF0, normF1;
    vec3d Um;
    double normM0, normM1;
    bool bconv = true;
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
    double alpha = tp.alpha;
    
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
    felog.printf(" rigid connector # %d (planar joint)\n", m_nID+1);
    felog.printf("                  CURRENT        REQUIRED\n");
    double pctn = 0;
    double gap = c.norm();
    double qap = ksi.norm();
    if (fabs(normF1) > 1e-10) pctn = fabs((normF1 - normF0)/normF1);
    if (m_atol) felog.printf("    force : %15le %15le\n", pctn, m_atol);
    else        felog.printf("    force : %15le        ***\n", pctn);
    if (m_gtol) felog.printf("    gap   : %15le %15le\n", gap, m_gtol);
    else        felog.printf("    gap   : %15le        ***\n", gap);
    double qctn = 0;
    if (fabs(normM1) > 1e-10) qctn = fabs((normM1 - normM0)/normM1);
    if (m_atol) felog.printf("    moment: %15le %15le\n", qctn, m_atol);
    else        felog.printf("    moment: %15le        ***\n", qctn);
    if (m_qtol) felog.printf("    angle : %15le %15le\n", qap, m_qtol);
    else        felog.printf("    angle : %15le        ***\n", qap);
    
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
    
    return bconv;
}

//-----------------------------------------------------------------------------
void FERigidPlanarJoint::Update(const FETimeInfo& tp)
{
    vec3d ra, rb;
    vec3d za, zb;
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
    double alpha = tp.alpha;
    
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
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
    m_qa0 = m_q0 - RBa.m_r0;
    m_qb0 = m_q0 - RBb.m_r0;
    
    m_ea0[0] = m_e0[0]; m_ea0[1] = m_e0[1]; m_ea0[2] = m_e0[2];
    m_eb0[0] = m_e0[0]; m_eb0[1] = m_e0[1]; m_eb0[2] = m_e0[2];
}
