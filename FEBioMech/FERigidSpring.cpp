//
//  FERigidSpring.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/11/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FERigidSpring.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidSpring, FERigidConnector);
    ADD_PARAMETER(m_k   , FE_PARAM_DOUBLE, "k"          );
    ADD_PARAMETER(m_a0  , FE_PARAM_VEC3D , "insertion_a");
    ADD_PARAMETER(m_b0  , FE_PARAM_VEC3D , "insertion_b");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidSpring::FERigidSpring(FEModel* pfem) : FERigidConnector(pfem)
{
    m_nID = m_ncount++;
    m_binit = false;
    m_L0 = 0;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidSpring::Init()
{
    if (m_binit) return true;
    
    // reset force
    m_F = vec3d(0,0,0);
    
    FEModel& fem = *GetFEModel();
    
    // When the rigid spring is read in, the ID's correspond to the rigid materials.
    // Now we want to make the ID's refer to the rigid body ID's
    
    FEMaterial* pm = fem.GetMaterial(m_nRBa-1);
    if (pm->IsRigid() == false)
    {
        felog.printbox("FATAL ERROR", "Rigid connector %d (spring) does not connect two rigid bodies\n", m_nID+1);
        return false;
    }
    m_nRBa = pm->GetRigidBodyID();
    
    pm = fem.GetMaterial(m_nRBb-1);
    if (pm->IsRigid() == false)
    {
        felog.printbox("FATAL ERROR", "Rigid connector %d (spring) does not connect two rigid bodies\n", m_nID);
        return false;
    }
    m_nRBb = pm->GetRigidBodyID();

	FERigidSystem& rs = *fem.GetRigidSystem();
    FERigidBody& ra = *rs.Object(m_nRBa);
    FERigidBody& rb = *rs.Object(m_nRBb);

    // get initial length of spring
    m_L0 = (m_b0 - m_a0).norm();
    // set spring insertions relative to rigid body center of mass
    m_qa0 = m_a0 - ra.m_r0;
    m_qb0 = m_b0 - rb.m_r0;
    
    m_binit = true;
    
    return true;
}

//-----------------------------------------------------------------------------
//! create a shallow copy
void FERigidSpring::ShallowCopy(DumpStream& dmp, bool bsave)
{
    if (bsave)
    {
        dmp << m_qa0 << m_qb0;
        dmp << m_F << m_L0;
    }
    else
    {
        dmp >> m_qa0 >> m_qb0;
        dmp >> m_F >> m_L0;
    }
}

//-----------------------------------------------------------------------------
void FERigidSpring::Serialize(DumpFile& ar)
{
	FERigidConnector::Serialize(ar);
    if (ar.IsSaving())
    {
		ar << m_binit;
        ar << m_qa0 << m_qb0;
        ar << m_L0 << m_k;
    }
    else
    {
		ar >> m_binit;
        ar >> m_qa0 >> m_qb0;
        ar >> m_L0 >> m_k;
    }
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidSpring::Residual(FEGlobalVector& R, const FETimePoint& tp)
{
    vector<double> fa(6);
    vector<double> fb(6);
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);

	double alpha = tp.alpha;

    // body A
    vec3d ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
    
    // body b
    vec3d rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
    
    vec3d c = rb + zb - ra - za;
    double L = c.norm();
    m_F = (L > 0) ? c*((1 - m_L0/L)*m_k) : c*m_k;
    
    fa[0] = m_F.x;
    fa[1] = m_F.y;
    fa[2] = m_F.z;
    
    fa[3] = za.y*m_F.z-za.z*m_F.y;
    fa[4] = za.z*m_F.x-za.x*m_F.z;
    fa[5] = za.x*m_F.y-za.y*m_F.x;
    
    fb[0] = -m_F.x;
    fb[1] = -m_F.y;
    fb[2] = -m_F.z;
    
    fb[3] = -zb.y*m_F.z+zb.z*m_F.y;
    fb[4] = -zb.z*m_F.x+zb.x*m_F.z;
    fb[5] = -zb.x*m_F.y+zb.y*m_F.x;
    
    for (int i=0; i<6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
    for (int i=0; i<6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidSpring::StiffnessMatrix(FESolver* psolver, const FETimePoint& tp)
{
	double alpha = tp.alpha;
    
    int j;
    
    vector<int> LM(12);
    matrix ke(12,12);
    ke.zero();
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
    // body A
    vec3d ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
    mat3d zahat; zahat.skew(za);
    mat3d zathat; zathat.skew(zat);
    
    // body b
    vec3d rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
    mat3d zbhat; zbhat.skew(zb);
    mat3d zbthat; zbthat.skew(zbt);
    
    vec3d c = rb + zb - ra - za;
    double L = c.norm();
    vec3d n = c/L;
    m_F = (L > 0) ? c*((1 - m_L0/L)*m_k) : c*m_k;
    mat3ds P;
    mat3dd I(1);
    P = (L > 0) ? I*(1 - m_L0/L) + dyad(n)*(m_L0/L) : I;
    
    mat3d K;
    
    // (1,1)
    K = P*(alpha*m_k);
    ke[0][0] = K[0][0]; ke[0][1] = K[0][1]; ke[0][2] = K[0][2];
    ke[1][0] = K[1][0]; ke[1][1] = K[1][1]; ke[1][2] = K[1][2];
    ke[2][0] = K[2][0]; ke[2][1] = K[2][1]; ke[2][2] = K[2][2];
    
    // (1,2)
    K = (P*zathat)*(-alpha*m_k);
    ke[0][3] = K[0][0]; ke[0][4] = K[0][1]; ke[0][5] = K[0][2];
    ke[1][3] = K[1][0]; ke[1][4] = K[1][1]; ke[1][5] = K[1][2];
    ke[2][3] = K[2][0]; ke[2][4] = K[2][1]; ke[2][5] = K[2][2];
    
    // (1,3)
    K = P*(-alpha*m_k);
    ke[0][6] = K[0][0]; ke[0][7] = K[0][1]; ke[0][8] = K[0][2];
    ke[1][6] = K[1][0]; ke[1][7] = K[1][1]; ke[1][8] = K[1][2];
    ke[2][6] = K[2][0]; ke[2][7] = K[2][1]; ke[2][8] = K[2][2];
    
    // (1,4)
    K = P*zbthat*(alpha*m_k);
    ke[0][9] = K[0][0]; ke[0][10] = K[0][1]; ke[0][11] = K[0][2];
    ke[1][9] = K[1][0]; ke[1][10] = K[1][1]; ke[1][11] = K[1][2];
    ke[2][9] = K[2][0]; ke[2][10] = K[2][1]; ke[2][11] = K[2][2];
    
    // (2,1)
    K = zahat*P*(alpha*m_k);
    ke[3][0] = K[0][0]; ke[3][1] = K[0][1]; ke[3][2] = K[0][2];
    ke[4][0] = K[1][0]; ke[4][1] = K[1][1]; ke[4][2] = K[1][2];
    ke[5][0] = K[2][0]; ke[5][1] = K[2][1]; ke[5][2] = K[2][2];
    
    // (2,2)
    K = (zahat*P*zathat)*(-alpha*m_k);
    ke[3][3] = K[0][0]; ke[3][4] = K[0][1]; ke[3][5] = K[0][2];
    ke[4][3] = K[1][0]; ke[4][4] = K[1][1]; ke[4][5] = K[1][2];
    ke[5][3] = K[2][0]; ke[5][4] = K[2][1]; ke[5][5] = K[2][2];
    
    // (2,3)
    K = zahat*P*(-alpha*m_k);
    ke[3][6] = K[0][0]; ke[3][7] = K[0][1]; ke[3][8] = K[0][2];
    ke[4][6] = K[1][0]; ke[4][7] = K[1][1]; ke[4][8] = K[1][2];
    ke[5][6] = K[2][0]; ke[5][7] = K[2][1]; ke[5][8] = K[2][2];
    
    // (2,4)
    K = (zahat*P*zbthat)*(alpha*m_k);
    ke[3][9] = K[0][0]; ke[3][10] = K[0][1]; ke[3][11] = K[0][2];
    ke[4][9] = K[1][0]; ke[4][10] = K[1][1]; ke[4][11] = K[1][2];
    ke[5][9] = K[2][0]; ke[5][10] = K[2][1]; ke[5][11] = K[2][2];
    
    
    // (3,1)
    K = P*(-alpha*m_k);
    ke[6][0] = K[0][0]; ke[6][1] = K[0][1]; ke[6][2] = K[0][2];
    ke[7][0] = K[1][0]; ke[7][1] = K[1][1]; ke[7][2] = K[1][2];
    ke[8][0] = K[2][0]; ke[8][1] = K[2][1]; ke[8][2] = K[2][2];
    
    // (3,2)
    K = (P*zathat)*(alpha*m_k);
    ke[6][3] = K[0][0]; ke[6][4] = K[0][1]; ke[6][5] = K[0][2];
    ke[7][3] = K[1][0]; ke[7][4] = K[1][1]; ke[7][5] = K[1][2];
    ke[8][3] = K[2][0]; ke[8][4] = K[2][1]; ke[8][5] = K[2][2];
    
    // (3,3)
    K = P*(alpha*m_k);
    ke[6][6] = K[0][0]; ke[6][7] = K[0][1]; ke[6][8] = K[0][2];
    ke[7][6] = K[1][0]; ke[7][7] = K[1][1]; ke[7][8] = K[1][2];
    ke[8][6] = K[2][0]; ke[8][7] = K[2][1]; ke[8][8] = K[2][2];
    
    // (3,4)
    K = P*zbthat*(-alpha*m_k);
    ke[6][9] = K[0][0]; ke[6][10] = K[0][1]; ke[6][11] = K[0][2];
    ke[7][9] = K[1][0]; ke[7][10] = K[1][1]; ke[7][11] = K[1][2];
    ke[8][9] = K[2][0]; ke[8][10] = K[2][1]; ke[8][11] = K[2][2];
    
    
    // (4,1)
    K = zbhat*P*(-alpha*m_k);
    ke[9 ][0] = K[0][0]; ke[ 9][1] = K[0][1]; ke[ 9][2] = K[0][2];
    ke[10][0] = K[1][0]; ke[10][1] = K[1][1]; ke[10][2] = K[1][2];
    ke[11][0] = K[2][0]; ke[11][1] = K[2][1]; ke[11][2] = K[2][2];
    
    // (4,2)
    K = (zbhat*P*zathat)*(alpha*m_k);
    ke[9 ][3] = K[0][0]; ke[ 9][4] = K[0][1]; ke[ 9][5] = K[0][2];
    ke[10][3] = K[1][0]; ke[10][4] = K[1][1]; ke[10][5] = K[1][2];
    ke[11][3] = K[2][0]; ke[11][4] = K[2][1]; ke[11][5] = K[2][2];
    
    // (4,3)
    K = zbhat*P*(alpha*m_k);
    ke[9 ][6] = K[0][0]; ke[ 9][7] = K[0][1]; ke[ 9][8] = K[0][2];
    ke[10][6] = K[1][0]; ke[10][7] = K[1][1]; ke[10][8] = K[1][2];
    ke[11][6] = K[2][0]; ke[11][7] = K[2][1]; ke[11][8] = K[2][2];
    
    // (4,4)
    K = (zbhat*P*zbthat)*(-alpha*m_k);
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
bool FERigidSpring::Augment(int naug, const FETimePoint& tp)
{
    return true;
}

//-----------------------------------------------------------------------------
void FERigidSpring::Update(const FETimePoint& tp)
{
    vec3d ra, rb, c;
    vec3d za, zb;
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);

	double alpha = tp.alpha;

    ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*alpha + zap*(1-alpha);
    
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*alpha + zbp*(1-alpha);
    
    c = rb + zb - ra - za;
    double L = c.norm();
    m_F = (L > 0) ? c*((1 - m_L0/L)*m_k) : c*m_k;
}

//-----------------------------------------------------------------------------
void FERigidSpring::Reset()
{
    m_F = vec3d(0,0,0);
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
    m_qa0 = m_a0 - RBa.m_r0;
    m_qb0 = m_b0 - RBb.m_r0;
}
