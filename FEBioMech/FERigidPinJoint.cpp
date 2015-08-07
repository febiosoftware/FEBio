//
//  FERigidPinJoint.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/25/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FERigidPinJoint.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidPinJoint, FENLConstraint);
    ADD_PARAMETER(m_atol, FE_PARAM_DOUBLE, "tolerance"     );
    ADD_PARAMETER(m_gtol, FE_PARAM_DOUBLE, "gaptol"        );
    ADD_PARAMETER(m_qtol, FE_PARAM_DOUBLE, "angtol"        );
    ADD_PARAMETER(m_eps , FE_PARAM_DOUBLE, "force_penalty" );
    ADD_PARAMETER(m_ups , FE_PARAM_DOUBLE, "moment_penalty");
    ADD_PARAMETER(m_nRBa, FE_PARAM_INT   , "body_a"        );
    ADD_PARAMETER(m_nRBb, FE_PARAM_INT   , "body_b"        );
    ADD_PARAMETER(m_q0  , FE_PARAM_VEC3D , "joint"  );
    ADD_PARAMETER(m_n0  , FE_PARAM_VEC3D , "axis"          );
    ADD_PARAMETER(m_naugmin,FE_PARAM_INT , "minaug"        );
    ADD_PARAMETER(m_naugmax,FE_PARAM_INT , "maxaug"        );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidPinJoint::FERigidPinJoint(FEModel* pfem) : FENLConstraint(pfem)
{
	static int count = 1;
	m_nID = count++;
	m_binit = false;
    m_atol = 0;
    m_gtol = 0;
    m_qtol = 0;
	m_naugmin = 0;
	m_naugmax = 10;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidPinJoint::Init()
{
	if (m_binit) return true;
    
	// reset force
	m_F = vec3d(0,0,0); m_L = vec3d(0,0,0);
    m_M = vec3d(0,0,0); m_U = vec3d(0,0,0);
    
	FEModel& fem = *GetFEModel();
    
	// When the rigid joint is read in, the ID's correspond to the rigid materials.
	// Now we want to make the ID's refer to the rigid body ID's
    
	FEMaterial* pm = fem.GetMaterial(m_nRBa);
	if (pm->IsRigid() == false)
	{
		felog.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", m_nID);
		return false;
	}
	m_nRBa = pm->GetRigidBodyID();
    
	pm = fem.GetMaterial(m_nRBb);
	if (pm->IsRigid() == false)
	{
		felog.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", m_nID);
		return false;
	}
	m_nRBb = pm->GetRigidBodyID();
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
	m_qa0 = m_q0 - RBa.m_r0;
	m_qb0 = m_q0 - RBb.m_r0;
    
    m_na0 = m_n0;
    m_nb0 = m_n0;
    
	m_binit = true;
    
	return true;
}

//-----------------------------------------------------------------------------
//! create a shallow copy
void FERigidPinJoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_q0 << m_qa0 << m_qb0;
		dmp << m_F << m_L;
		dmp << m_n0 << m_na0 << m_nb0;
		dmp << m_M << m_U;
	}
	else
	{
		dmp >> m_q0 >> m_qa0 >> m_qb0;
		dmp >> m_F >> m_L;
		dmp >> m_n0 >> m_na0 >> m_nb0;
		dmp >> m_M >> m_U;
	}
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidPinJoint::Residual(FEGlobalVector& R, const FETimePoint& tp)
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
	vec3d nat = m_na0; RBa.m_qt.RotateVector(nat);
	vec3d nap = m_na0; RBa.m_qp.RotateVector(nap);
    vec3d na = nat*alpha + nap*(1-alpha);
    
	// body b
    vec3d rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
	vec3d nbt = m_nb0; RBb.m_qt.RotateVector(nbt);
	vec3d nbp = m_nb0; RBb.m_qp.RotateVector(nbp);
    vec3d nb = nbt*alpha + nbp*(1-alpha);
    
	vec3d c = rb + zb - ra - za;
	m_F = m_L + c*m_eps;
    
    vec3d q = na ^ nb;
    m_M = m_U + q*m_ups;
    
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
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidPinJoint::StiffnessMatrix(FESolver* psolver, const FETimePoint& tp)
{
	double alpha = tp.alpha;
    
	int j, k;
    
	vector<int> LM(12);
	matrix ke(12,12);
	ke.zero();
	vec3d a, at, ap;
    vec3d na, nat, nap, nb, nbt, nbp;
    
	double y1[3][3], y2[3][3], y1h[3][3], y2h[3][3];
    double y11[3][3], y12[3][3], y21[3][3], y22[3][3];
    
	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
	// body A
    at = m_qa0; RBa.m_qt.RotateVector(at);
    ap = m_qa0; RBa.m_qp.RotateVector(ap);
    a = at*alpha + ap*(1-alpha);
    
	y1h[0][0] =    0; y1h[0][1] =  a.z; y1h[0][2] = -a.y;
	y1h[1][0] = -a.z; y1h[1][1] =    0; y1h[1][2] =  a.x;
	y1h[2][0] =  a.y; y1h[2][1] = -a.x; y1h[2][2] =    0;
    
	y1[0][0] =     0; y1[0][1] =  at.z; y1[0][2] = -at.y;
	y1[1][0] = -at.z; y1[1][1] =     0; y1[1][2] =  at.x;
	y1[2][0] =  at.y; y1[2][1] = -at.x; y1[2][2] =     0;
    
	// body b
    at = m_qb0; RBb.m_qt.RotateVector(at);
    ap = m_qb0; RBb.m_qp.RotateVector(ap);
    a = at*alpha + ap*(1-alpha);
    
	y2h[0][0] =    0; y2h[0][1] =  a.z; y2h[0][2] = -a.y;
	y2h[1][0] = -a.z; y2h[1][1] =    0; y2h[1][2] =  a.x;
	y2h[2][0] =  a.y; y2h[2][1] = -a.x; y2h[2][2] =    0;
    
	y2[0][0] =     0; y2[0][1] =  at.z; y2[0][2] = -at.y;
	y2[1][0] = -at.z; y2[1][1] =     0; y2[1][2] =  at.x;
	y2[2][0] =  at.y; y2[2][1] = -at.x; y2[2][2] =     0;
    
	nat = m_na0; RBa.m_qt.RotateVector(nat);
	nap = m_na0; RBa.m_qp.RotateVector(nap);
    na = nat*alpha + nap*(1-alpha);
	nbt = m_nb0; RBb.m_qt.RotateVector(nbt);
	nbp = m_nb0; RBb.m_qp.RotateVector(nbp);
    nb = nbt*alpha + nbp*(1-alpha);
    mat3d Nab = ((na*nb)*mat3dd(1) - (na & nb))*(m_ups*alpha);
    mat3d Nba = Nab.transpose();
    
	for (j=0; j<3; ++j)
	{
		y11[j][0] = y1h[0][j]*y1[0][0]+y1h[1][j]*y1[1][0]+y1h[2][j]*y1[2][0];
		y11[j][1] = y1h[0][j]*y1[0][1]+y1h[1][j]*y1[1][1]+y1h[2][j]*y1[2][1];
		y11[j][2] = y1h[0][j]*y1[0][2]+y1h[1][j]*y1[1][2]+y1h[2][j]*y1[2][2];
        
		y12[j][0] = y1h[0][j]*y2[0][0]+y1h[1][j]*y2[1][0]+y1h[2][j]*y2[2][0];
		y12[j][1] = y1h[0][j]*y2[0][1]+y1h[1][j]*y2[1][1]+y1h[2][j]*y2[2][1];
		y12[j][2] = y1h[0][j]*y2[0][2]+y1h[1][j]*y2[1][2]+y1h[2][j]*y2[2][2];
        
		y21[j][0] = y2h[0][j]*y1[0][0]+y2h[1][j]*y1[1][0]+y2h[2][j]*y1[2][0];
		y21[j][1] = y2h[0][j]*y1[0][1]+y2h[1][j]*y1[1][1]+y2h[2][j]*y1[2][1];
		y21[j][2] = y2h[0][j]*y1[0][2]+y2h[1][j]*y1[1][2]+y2h[2][j]*y1[2][2];
        
		y22[j][0] = y2h[0][j]*y2[0][0]+y2h[1][j]*y2[1][0]+y2h[2][j]*y2[2][0];
		y22[j][1] = y2h[0][j]*y2[0][1]+y2h[1][j]*y2[1][1]+y2h[2][j]*y2[2][1];
		y22[j][2] = y2h[0][j]*y2[0][2]+y2h[1][j]*y2[1][2]+y2h[2][j]*y2[2][2];
	}
    
	ke[0][0] = ke[1][1] = ke[2][2] =  1;
	ke[0][6] = ke[1][7] = ke[2][8] = -1;
	ke[6][0] = ke[7][1] = ke[8][2] = -1;
	ke[6][6] = ke[7][7] = ke[8][8] =  1;
    
	ke[0][3] = y1[0][0]; ke[0][4] = y1[0][1]; ke[0][5] = y1[0][2];
	ke[1][3] = y1[1][0]; ke[1][4] = y1[1][1]; ke[1][5] = y1[1][2];
	ke[2][3] = y1[2][0]; ke[2][4] = y1[2][1]; ke[2][5] = y1[2][2];
    
	ke[0][9] = -y2[0][0]; ke[0][10] = -y2[0][1]; ke[0][11] = -y2[0][2];
	ke[1][9] = -y2[1][0]; ke[1][10] = -y2[1][1]; ke[1][11] = -y2[1][2];
	ke[2][9] = -y2[2][0]; ke[2][10] = -y2[2][1]; ke[2][11] = -y2[2][2];
    
	ke[3][0] = y1h[0][0]; ke[3][1] = y1h[1][0]; ke[3][2] = y1h[2][0];
	ke[4][0] = y1h[0][1]; ke[4][1] = y1h[1][1]; ke[4][2] = y1h[2][1];
	ke[5][0] = y1h[0][2]; ke[5][1] = y1h[1][2]; ke[5][2] = y1h[2][2];
    
	ke[3][3] = y11[0][0]; ke[3][4] = y11[0][1]; ke[3][5] = y11[0][2];
	ke[4][3] = y11[1][0]; ke[4][4] = y11[1][1]; ke[4][5] = y11[1][2];
	ke[5][3] = y11[2][0]; ke[5][4] = y11[2][1]; ke[5][5] = y11[2][2];
    
	ke[3][6] = -y1h[0][0]; ke[3][7] = -y1h[1][0]; ke[3][8] = -y1h[2][0];
	ke[4][6] = -y1h[0][1]; ke[4][7] = -y1h[1][1]; ke[4][8] = -y1h[2][1];
	ke[5][6] = -y1h[0][2]; ke[5][7] = -y1h[1][2]; ke[5][8] = -y1h[2][2];
    
	ke[3][9] = -y12[0][0]; ke[3][10] = -y12[0][1]; ke[3][11] = -y12[0][2];
	ke[4][9] = -y12[1][0]; ke[4][10] = -y12[1][1]; ke[4][11] = -y12[1][2];
	ke[5][9] = -y12[2][0]; ke[5][10] = -y12[2][1]; ke[5][11] = -y12[2][2];
    
	ke[6][3] = -y1[0][0]; ke[6][4] = -y1[0][1]; ke[6][5] = -y1[0][2];
	ke[7][3] = -y1[1][0]; ke[7][4] = -y1[1][1]; ke[7][5] = -y1[1][2];
	ke[8][3] = -y1[2][0]; ke[8][4] = -y1[2][1]; ke[8][5] = -y1[2][2];
    
	ke[6][9] = y2[0][0]; ke[6][10] = y2[0][1]; ke[6][11] = y2[0][2];
	ke[7][9] = y2[1][0]; ke[7][10] = y2[1][1]; ke[7][11] = y2[1][2];
	ke[8][9] = y2[2][0]; ke[8][10] = y2[2][1]; ke[8][11] = y2[2][2];
    
	ke[ 9][0] = -y2h[0][0]; ke[ 9][1] = -y2h[1][0]; ke[ 9][2] = -y2h[2][0];
	ke[10][0] = -y2h[0][1]; ke[10][1] = -y2h[1][1]; ke[10][2] = -y2h[2][1];
	ke[11][0] = -y2h[0][2]; ke[11][1] = -y2h[1][2]; ke[11][2] = -y2h[2][2];
    
	ke[ 9][3] = -y21[0][0]; ke[ 9][4] = -y21[0][1]; ke[ 9][5] = -y21[0][2];
	ke[10][3] = -y21[1][0]; ke[10][4] = -y21[1][1]; ke[10][5] = -y21[1][2];
	ke[11][3] = -y21[2][0]; ke[11][4] = -y21[2][1]; ke[11][5] = -y21[2][2];
    
	ke[ 9][6] = y2h[0][0]; ke[ 9][7] = y2h[1][0]; ke[ 9][8] = y2h[2][0];
	ke[10][6] = y2h[0][1]; ke[10][7] = y2h[1][1]; ke[10][8] = y2h[2][1];
	ke[11][6] = y2h[0][2]; ke[11][7] = y2h[1][2]; ke[11][8] = y2h[2][2];
    
	ke[ 9][9] = y22[0][0]; ke[ 9][10] = y22[0][1]; ke[ 9][11] = y22[0][2];
	ke[10][9] = y22[1][0]; ke[10][10] = y22[1][1]; ke[10][11] = y22[1][2];
	ke[11][9] = y22[2][0]; ke[11][10] = y22[2][1]; ke[11][11] = y22[2][2];
    
	for (j=0; j<12; ++j)
		for (k=0; k<12; ++k)
		{
			ke[j][k] *= m_eps*alpha;
		}
    
	ke[3][3] += Nab(0,0); ke[3][4] += Nab(0,1); ke[3][5] += Nab(0,2);
	ke[4][3] += Nab(1,0); ke[4][4] += Nab(1,1); ke[4][5] += Nab(1,2);
	ke[5][3] += Nab(2,0); ke[5][4] += Nab(2,1); ke[5][5] += Nab(2,2);
    
	ke[3][9] += -Nba(0,0); ke[3][10] += -Nba(0,1); ke[3][11] += -Nba(0,2);
	ke[4][9] += -Nba(1,0); ke[4][10] += -Nba(1,1); ke[4][11] += -Nba(1,2);
	ke[5][9] += -Nba(2,0); ke[5][10] += -Nba(2,1); ke[5][11] += -Nba(2,2);
    
	ke[ 9][3] += -Nab(0,0); ke[ 9][4] += -Nab(1,0); ke[ 9][5] += -Nab(2,0);
	ke[10][3] += -Nab(0,1); ke[10][4] += -Nab(1,1); ke[10][5] += -Nab(2,1);
	ke[11][3] += -Nab(0,2); ke[11][4] += -Nab(1,2); ke[11][5] += -Nab(2,2);
    
	ke[ 9][9] += Nba(0,0); ke[ 9][10] += Nba(0,1); ke[ 9][11] += Nba(0,2);
	ke[10][9] += Nba(1,0); ke[10][10] += Nba(1,1); ke[10][11] += Nba(1,2);
	ke[11][9] += Nba(2,0); ke[11][10] += Nba(2,1); ke[11][11] += Nba(2,2);
    
	for (j=0; j<6; ++j)
	{
		LM[j  ] = RBa.m_LM[j];
		LM[j+6] = RBb.m_LM[j];
	}
    
	psolver->AssembleStiffness(LM, ke);
}

//-----------------------------------------------------------------------------
bool FERigidPinJoint::Augment(int naug, const FETimePoint& tp)
{
	vec3d ra, rb, qa, qb, c,  Lm;
    vec3d za, zb;
	double normF0, normF1;
	vec3d na, nat, nap, nb, nbt, nbp;
    vec3d q,  Um;
	double normM0, normM1;
	bool bconv = true;

	double alpha = tp.alpha;

	FERigidSystem& rs = *GetFEModel()->GetRigidSystem();
    FERigidBody& RBa = *rs.Object(m_nRBa);
    FERigidBody& RBb = *rs.Object(m_nRBb);
    
    ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*alpha + zap*(1-alpha);
    
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*alpha + zbp*(1-alpha);
    
	c = rb + zb - ra - za;
    
	normF0 = sqrt(m_L*m_L);
    
	// calculate trial multiplier
	Lm = m_L + c*m_eps;
    
	normF1 = sqrt(Lm*Lm);
    
	nat = m_na0; RBa.m_qt.RotateVector(nat);
	nap = m_na0; RBa.m_qp.RotateVector(nap);
    na = nat*alpha + nap*(1-alpha);
	nbt = m_nb0; RBb.m_qt.RotateVector(nbt);
	nbp = m_nb0; RBb.m_qp.RotateVector(nbp);
    nb = nbt*alpha + nbp*(1-alpha);
    
    q = na ^ nb;
    
	normM0 = sqrt(m_U*m_U);
    
	// calculate trial multiplier
	Um = m_U + q*m_ups;
    
	normM1 = sqrt(Um*Um);
    
	// check convergence of constraints
	felog.printf(" rigid joint # %d\n", m_nID);
	felog.printf("                  CURRENT        REQUIRED\n");
	double pctn = 0;
    double gap = c.norm();
    double qap = q.norm();
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
void FERigidPinJoint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_nRBa << m_nRBb;
		ar << m_q0 << m_qa0 << m_qb0;
		ar << m_F << m_L << m_eps << m_atol;
		ar << m_n0 << m_na0 << m_nb0;
		ar << m_M << m_U << m_ups;
	}
	else
	{
		ar >> m_nID;
		ar >> m_nRBa >> m_nRBb;
		ar >> m_q0 >> m_qa0 >> m_qb0;
		ar >> m_F >> m_L >> m_eps >> m_atol;
		ar >> m_n0 >> m_na0 >> m_nb0;
		ar >> m_M >> m_U >> m_ups;
	}
}

//-----------------------------------------------------------------------------
void FERigidPinJoint::Update(const FETimePoint& tp)
{
	vec3d ra, rb, c;
    vec3d za, zb;
	vec3d na, nat, nap, nb, nbt, nbp;
    
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
    
	m_F = m_L + c*m_eps;
    
	nat = m_na0; RBa.m_qt.RotateVector(nat);
	nap = m_na0; RBa.m_qp.RotateVector(nap);
    na = nat*alpha + nap*(1-alpha);
	nbt = m_nb0; RBb.m_qt.RotateVector(nbt);
	nbp = m_nb0; RBb.m_qp.RotateVector(nbp);
    nb = nbt*alpha + nbp*(1-alpha);
    
    vec3d q = na ^ nb;
    m_M = m_U + q*m_ups;
    
}

//-----------------------------------------------------------------------------
void FERigidPinJoint::Reset()
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
    
    m_na0 = m_n0;
    m_nb0 = m_n0;
}
