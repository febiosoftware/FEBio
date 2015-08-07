//
//  FERigidDamper.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/12/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FERigidDamper.h"
#include "stdafx.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidDamper, FERigidConnector);
ADD_PARAMETER(m_c   , FE_PARAM_DOUBLE, "c"          );
ADD_PARAMETER(m_a0  , FE_PARAM_VEC3D , "insertion_a");
ADD_PARAMETER(m_b0  , FE_PARAM_VEC3D , "insertion_b");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidDamper::FERigidDamper(FEModel* pfem) : FERigidConnector(pfem)
{
    m_nID = m_ncount++;
    m_binit = false;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidDamper::Init()
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
        felog.printbox("FATAL ERROR", "Rigid connector %d (damper) does not connect two rigid bodies\n", m_nID+1);
        return false;
    }
    m_nRBa = pm->GetRigidBodyID();
    
    pm = fem.GetMaterial(m_nRBb-1);
    if (pm->IsRigid() == false)
    {
        felog.printbox("FATAL ERROR", "Rigid connector %d (damper) does not connect two rigid bodies\n", m_nID+1);
        return false;
    }
    m_nRBb = pm->GetRigidBodyID();
    
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rigid.Object(m_nRBa);
    FERigidBody& RBb = *rigid.Object(m_nRBb);
    
    // set spring insertions relative to rigid body center of mass
    m_qa0 = m_a0 - RBa.m_r0;
    m_qb0 = m_b0 - RBb.m_r0;
    
    m_binit = true;
    
    return true;
}

//-----------------------------------------------------------------------------
//! create a shallow copy
void FERigidDamper::ShallowCopy(DumpStream& dmp, bool bsave)
{
    if (bsave)
    {
        dmp << m_qa0 << m_qb0;
        dmp << m_F;
    }
    else
    {
        dmp >> m_qa0 >> m_qb0;
        dmp >> m_F;
    }
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidDamper::Residual(FEGlobalVector& R, const FETimePoint& tp)
{
    vector<double> fa(6);
    vector<double> fb(6);
    
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rigid.Object(m_nRBa);
    FERigidBody& RBb = *rigid.Object(m_nRBb);

	double alpha = tp.alpha;
    
    // body A
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
    vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
    vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
    vec3d va = vat*alpha + vap*(1-alpha);
    
    // body b
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
    vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
    vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
    vec3d vb = vbt*alpha + vbp*(1-alpha);
    
    m_F = (vb - va)*m_c;
    
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
void FERigidDamper::StiffnessMatrix(FESolver* psolver, const FETimePoint& tp)
{
	double alpha = tp.alpha;
	double beta  = tp.beta;
	double gamma = tp.gamma;
    
    // get time increment
    FEModel& fem = *GetFEModel();
    double dt = fem.GetCurrentStep()->m_dt;
    
    int j;
    
    vector<int> LM(12);
    matrix ke(12,12);
    ke.zero();
    
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rigid.Object(m_nRBa);
    FERigidBody& RBb = *rigid.Object(m_nRBb);
    
    mat3dd I(1);
    
    // body A
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
    mat3d zahat; zahat.skew(za);
    mat3d zathat; zathat.skew(zat);
    vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
    vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
    vec3d va = vat*alpha + vap*(1-alpha);
    quatd qai = RBa.m_qt*RBa.m_qp.Inverse(); qai.MakeUnit();
    vec3d cai = qai.GetVector()*(2*tan(qai.GetAngle()/2));
    mat3d Ta = I + skew(cai)/2 + dyad(cai)/4;
    
    // body b
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
    mat3d zbhat; zbhat.skew(zb);
    mat3d zbthat; zbthat.skew(zbt);
    vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
    vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
    vec3d vb = vbt*alpha + vbp*(1-alpha);
    quatd qbi = RBb.m_qt*RBb.m_qp.Inverse(); qbi.MakeUnit();
    vec3d cbi = qbi.GetVector()*(2*tan(qbi.GetAngle()/2));
    mat3d Tb = I + skew(cbi)/2 + dyad(cbi)/4;
    
    m_F = (vb - va)*m_c;
    
    mat3ds A = I*(gamma/beta/dt);
    mat3d Ba = zathat*Ta.transpose()*(gamma/beta/dt) + skew(RBa.m_wt)*zathat;
    mat3d Bb = zbthat*Tb.transpose()*(gamma/beta/dt) + skew(RBb.m_wt)*zbthat;
    
    mat3d K;
    
    // (1,1)
    K = A*(alpha*m_c);
    ke[0][0] = K[0][0]; ke[0][1] = K[0][1]; ke[0][2] = K[0][2];
    ke[1][0] = K[1][0]; ke[1][1] = K[1][1]; ke[1][2] = K[1][2];
    ke[2][0] = K[2][0]; ke[2][1] = K[2][1]; ke[2][2] = K[2][2];
    
    // (1,2)
    K = Ba*(-alpha*m_c);
    ke[0][3] = K[0][0]; ke[0][4] = K[0][1]; ke[0][5] = K[0][2];
    ke[1][3] = K[1][0]; ke[1][4] = K[1][1]; ke[1][5] = K[1][2];
    ke[2][3] = K[2][0]; ke[2][4] = K[2][1]; ke[2][5] = K[2][2];
    
    // (1,3)
    K = A*(-alpha*m_c);
    ke[0][6] = K[0][0]; ke[0][7] = K[0][1]; ke[0][8] = K[0][2];
    ke[1][6] = K[1][0]; ke[1][7] = K[1][1]; ke[1][8] = K[1][2];
    ke[2][6] = K[2][0]; ke[2][7] = K[2][1]; ke[2][8] = K[2][2];
    
    // (1,4)
    K = Bb*(alpha*m_c);
    ke[0][9] = K[0][0]; ke[0][10] = K[0][1]; ke[0][11] = K[0][2];
    ke[1][9] = K[1][0]; ke[1][10] = K[1][1]; ke[1][11] = K[1][2];
    ke[2][9] = K[2][0]; ke[2][10] = K[2][1]; ke[2][11] = K[2][2];
    
    // (2,1)
    K = zahat*A*(alpha*m_c);
    ke[3][0] = K[0][0]; ke[3][1] = K[0][1]; ke[3][2] = K[0][2];
    ke[4][0] = K[1][0]; ke[4][1] = K[1][1]; ke[4][2] = K[1][2];
    ke[5][0] = K[2][0]; ke[5][1] = K[2][1]; ke[5][2] = K[2][2];
    
    // (2,2)
    K = (zahat*Ba)*(-alpha*m_c);
    ke[3][3] = K[0][0]; ke[3][4] = K[0][1]; ke[3][5] = K[0][2];
    ke[4][3] = K[1][0]; ke[4][4] = K[1][1]; ke[4][5] = K[1][2];
    ke[5][3] = K[2][0]; ke[5][4] = K[2][1]; ke[5][5] = K[2][2];
    
    // (2,3)
    K = zahat*A*(-alpha*m_c);
    ke[3][6] = K[0][0]; ke[3][7] = K[0][1]; ke[3][8] = K[0][2];
    ke[4][6] = K[1][0]; ke[4][7] = K[1][1]; ke[4][8] = K[1][2];
    ke[5][6] = K[2][0]; ke[5][7] = K[2][1]; ke[5][8] = K[2][2];
    
    // (2,4)
    K = (zahat*Bb)*(alpha*m_c);
    ke[3][9] = K[0][0]; ke[3][10] = K[0][1]; ke[3][11] = K[0][2];
    ke[4][9] = K[1][0]; ke[4][10] = K[1][1]; ke[4][11] = K[1][2];
    ke[5][9] = K[2][0]; ke[5][10] = K[2][1]; ke[5][11] = K[2][2];
    
    
    // (3,1)
    K = A*(-alpha*m_c);
    ke[6][0] = K[0][0]; ke[6][1] = K[0][1]; ke[6][2] = K[0][2];
    ke[7][0] = K[1][0]; ke[7][1] = K[1][1]; ke[7][2] = K[1][2];
    ke[8][0] = K[2][0]; ke[8][1] = K[2][1]; ke[8][2] = K[2][2];
    
    // (3,2)
    K = Ba*(alpha*m_c);
    ke[6][3] = K[0][0]; ke[6][4] = K[0][1]; ke[6][5] = K[0][2];
    ke[7][3] = K[1][0]; ke[7][4] = K[1][1]; ke[7][5] = K[1][2];
    ke[8][3] = K[2][0]; ke[8][4] = K[2][1]; ke[8][5] = K[2][2];
    
    // (3,3)
    K = A*(alpha*m_c);
    ke[6][6] = K[0][0]; ke[6][7] = K[0][1]; ke[6][8] = K[0][2];
    ke[7][6] = K[1][0]; ke[7][7] = K[1][1]; ke[7][8] = K[1][2];
    ke[8][6] = K[2][0]; ke[8][7] = K[2][1]; ke[8][8] = K[2][2];
    
    // (3,4)
    K = Bb*(-alpha*m_c);
    ke[6][9] = K[0][0]; ke[6][10] = K[0][1]; ke[6][11] = K[0][2];
    ke[7][9] = K[1][0]; ke[7][10] = K[1][1]; ke[7][11] = K[1][2];
    ke[8][9] = K[2][0]; ke[8][10] = K[2][1]; ke[8][11] = K[2][2];
    
    
    // (4,1)
    K = zbhat*A*(-alpha*m_c);
    ke[9 ][0] = K[0][0]; ke[ 9][1] = K[0][1]; ke[ 9][2] = K[0][2];
    ke[10][0] = K[1][0]; ke[10][1] = K[1][1]; ke[10][2] = K[1][2];
    ke[11][0] = K[2][0]; ke[11][1] = K[2][1]; ke[11][2] = K[2][2];
    
    // (4,2)
    K = (zbhat*Ba)*(alpha*m_c);
    ke[9 ][3] = K[0][0]; ke[ 9][4] = K[0][1]; ke[ 9][5] = K[0][2];
    ke[10][3] = K[1][0]; ke[10][4] = K[1][1]; ke[10][5] = K[1][2];
    ke[11][3] = K[2][0]; ke[11][4] = K[2][1]; ke[11][5] = K[2][2];
    
    // (4,3)
    K = zbhat*A*(alpha*m_c);
    ke[9 ][6] = K[0][0]; ke[ 9][7] = K[0][1]; ke[ 9][8] = K[0][2];
    ke[10][6] = K[1][0]; ke[10][7] = K[1][1]; ke[10][8] = K[1][2];
    ke[11][6] = K[2][0]; ke[11][7] = K[2][1]; ke[11][8] = K[2][2];
    
    // (4,4)
    K = (zbhat*Bb)*(-alpha*m_c);
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
bool FERigidDamper::Augment(int naug, const FETimePoint& tp)
{
    return true;
}

//-----------------------------------------------------------------------------
void FERigidDamper::Serialize(DumpFile& ar)
{
    if (ar.IsSaving())
    {
        ar << m_nID;
        ar << m_nRBa << m_nRBb;
        ar << m_qa0 << m_qb0;
        ar << m_F << m_c;
    }
    else
    {
        ar >> m_nID;
        ar >> m_nRBa >> m_nRBb;
        ar >> m_qa0 >> m_qb0;
        ar >> m_F >> m_c;
    }
}

//-----------------------------------------------------------------------------
void FERigidDamper::Update(const FETimePoint& tp)
{
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rigid.Object(m_nRBa);
    FERigidBody& RBb = *rigid.Object(m_nRBb);

	double alpha = tp.alpha;

    // body A
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
    vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
    vec3d va = vat*alpha + vap*(1-alpha);
    
    // body b
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
    vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
    vec3d vb = vbt*alpha + vbp*(1-alpha);
    
    m_F = (vb - va)*m_c;
}

//-----------------------------------------------------------------------------
void FERigidDamper::Reset()
{
    m_F = vec3d(0,0,0);
    
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rigid.Object(m_nRBa);
    FERigidBody& RBb = *rigid.Object(m_nRBb);
    
    m_qa0 = m_a0 - RBa.m_r0;
    m_qb0 = m_b0 - RBb.m_r0;
}
