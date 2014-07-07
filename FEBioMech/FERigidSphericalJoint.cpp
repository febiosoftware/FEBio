// FERigidSphericalJoint.cpp: implementation of the FERigidSphericalJoint class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidSphericalJoint.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FESolidSolver2.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidSphericalJoint, FENLConstraint);
ADD_PARAMETER(m_atol, FE_PARAM_DOUBLE, "tolerance");
ADD_PARAMETER(m_gtol, FE_PARAM_DOUBLE, "gaptol"   );
ADD_PARAMETER(m_eps , FE_PARAM_DOUBLE, "penalty"  );
ADD_PARAMETER(m_nRBa, FE_PARAM_INT   , "body_a"   );
ADD_PARAMETER(m_nRBb, FE_PARAM_INT   , "body_b"   );
ADD_PARAMETER(m_q0  , FE_PARAM_VEC3D , "joint"    );
ADD_PARAMETER(m_naugmin,FE_PARAM_INT, "minaug"    );
ADD_PARAMETER(m_naugmax,FE_PARAM_INT, "maxaug"    );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidSphericalJoint::FERigidSphericalJoint(FEModel* pfem) : FENLConstraint(pfem)
{
	static int count = 1;
	m_nID = count++;
	m_binit = false;
    m_atol = 0;
    m_gtol = 0;
	m_naugmin = 0;
	m_naugmax = 10;
    m_alpha = 0.5;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidSphericalJoint::Init()
{
	if (m_binit) return true;
    
	// reset force
	m_F = vec3d(0,0,0);
    
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
    
	FERigidBody& ra = dynamic_cast<FERigidBody&>(*fem.Object(m_nRBa));
	FERigidBody& rb = dynamic_cast<FERigidBody&>(*fem.Object(m_nRBb));
    
	m_qa0 = m_q0 - ra.m_r0;
	m_qb0 = m_q0 - rb.m_r0;
    
	m_binit = true;
    
	return true;
}

//-----------------------------------------------------------------------------
//! create a shallow copy
void FERigidSphericalJoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_q0 << m_qa0 << m_qb0;
		dmp << m_F << m_L;
	}
	else
	{
		dmp >> m_q0 >> m_qa0 >> m_qb0;
		dmp >> m_F >> m_L;
	}
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidSphericalJoint::Residual(FEGlobalVector& R)
{
	vector<double> fa(6);
	vector<double> fb(6);
    
	FERigidBody& RBa = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBa));
	FERigidBody& RBb = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBb));
    
	// body A
    vec3d ra = RBa.m_rt*m_alpha + RBa.m_rp*(1-m_alpha);
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*m_alpha + zap*(1-m_alpha);
    
	// body b
    vec3d rb = RBb.m_rt*m_alpha + RBb.m_rp*(1-m_alpha);
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*m_alpha + zbp*(1-m_alpha);
    
	vec3d c = rb + zb - ra - za;
	m_F = m_L + c*m_eps;
    
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
void FERigidSphericalJoint::StiffnessMatrix(FESolver* psolver)
{
    // get m_alpha from solver
    m_alpha = dynamic_cast<FESolidSolver2*>(psolver)->m_alpha;
    
	int j, k;
    
	vector<int> LM(12);
	matrix ke(12,12);
	ke.zero();
    vec3d a;
    
	double y1[3][3], y2[3][3], y1h[3][3], y2h[3][3];
    double y11[3][3], y12[3][3], y21[3][3], y22[3][3];
    
	FERigidBody& RBa = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBa));
	FERigidBody& RBb = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBb));
    
	// body A
    vec3d at = m_qa0; RBa.m_qt.RotateVector(at);
    vec3d ap = m_qa0; RBa.m_qp.RotateVector(ap);
    a = at*m_alpha + ap*(1-m_alpha);
    
	y1h[0][0] =    0; y1h[0][1] =  a.z; y1h[0][2] = -a.y;
	y1h[1][0] = -a.z; y1h[1][1] =    0; y1h[1][2] =  a.x;
	y1h[2][0] =  a.y; y1h[2][1] = -a.x; y1h[2][2] =    0;
    
	y1[0][0] =     0; y1[0][1] =  at.z; y1[0][2] = -at.y;
	y1[1][0] = -at.z; y1[1][1] =     0; y1[1][2] =  at.x;
	y1[2][0] =  at.y; y1[2][1] = -at.x; y1[2][2] =     0;
    
	// body b
    at = m_qb0; RBb.m_qt.RotateVector(at);
    ap = m_qb0; RBb.m_qp.RotateVector(ap);
    a = at*m_alpha + ap*(1-m_alpha);
    
	y2h[0][0] =    0; y2h[0][1] =  a.z; y2h[0][2] = -a.y;
	y2h[1][0] = -a.z; y2h[1][1] =    0; y2h[1][2] =  a.x;
	y2h[2][0] =  a.y; y2h[2][1] = -a.x; y2h[2][2] =    0;
    
	y2[0][0] =     0; y2[0][1] =  at.z; y2[0][2] = -at.y;
	y2[1][0] = -at.z; y2[1][1] =     0; y2[1][2] =  at.x;
	y2[2][0] =  at.y; y2[2][1] = -at.x; y2[2][2] =     0;
    
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
			ke[j][k] *= m_eps*m_alpha;
		}
    
	for (j=0; j<6; ++j)
	{
		LM[j  ] = RBa.m_LM[j];
		LM[j+6] = RBb.m_LM[j];
	}
    
	psolver->AssembleStiffness(LM, ke);
}

//-----------------------------------------------------------------------------
bool FERigidSphericalJoint::Augment(int naug)
{
	vec3d ra, rb, c,  Lm;
    vec3d za, zb;
	double normF0, normF1;
	bool bconv = true;
    
	FERigidBody& RBa = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBa));
	FERigidBody& RBb = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBb));
    
    ra = RBa.m_rt*m_alpha + RBa.m_rp*(1-m_alpha);
    rb = RBb.m_rt*m_alpha + RBb.m_rp*(1-m_alpha);
    
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*m_alpha + zap*(1-m_alpha);

    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*m_alpha + zbp*(1-m_alpha);
    
	c = rb + zb - ra - za;
    
	normF0 = m_L.norm();
    
	// calculate trial multiplier
	Lm = m_L + c*m_eps;
    
	normF1 = Lm.norm();
    
	// check convergence of constraints
	felog.printf(" rigid joint # %d\n", m_nID);
	felog.printf("                  CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normF1) > 1e-10) pctn = fabs((normF1 - normF0)/normF1);
    double gap = c.norm();
    if (m_atol) felog.printf("    force : %15le %15le\n", pctn, m_atol);
    else        felog.printf("    force : %15le        ***\n", pctn);
    if (m_gtol) felog.printf("    gap   : %15le %15le\n", gap, m_gtol);
    else        felog.printf("    gap   : %15le        ***\n", gap);
    
    if ((m_atol > 0) && (pctn >= m_atol)) bconv = false;
    if ((m_gtol > 0) && (gap  >= m_gtol)) bconv = false;
	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;
    
    // update multiplier
	if (!bconv) m_L = Lm;
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FERigidSphericalJoint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_nRBa << m_nRBb;
		ar << m_q0 << m_qa0 << m_qb0;
		ar << m_F << m_L << m_eps << m_atol;
	}
	else
	{
		ar >> m_nID;
		ar >> m_nRBa >> m_nRBb;
		ar >> m_q0 >> m_qa0 >> m_qb0;
		ar >> m_F >> m_L >> m_eps >> m_atol;
	}
}

//-----------------------------------------------------------------------------
void FERigidSphericalJoint::Update()
{
	vec3d ra, rb, c;
    vec3d za, zb;
    
	FERigidBody& RBa = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBa));
	FERigidBody& RBb = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBb));
    
    ra = RBa.m_rt*m_alpha + RBa.m_rp*(1-m_alpha);
    rb = RBb.m_rt*m_alpha + RBb.m_rp*(1-m_alpha);
    
    vec3d zat = m_qa0; RBa.m_qt.RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*m_alpha + zap*(1-m_alpha);
    
    vec3d zbt = m_qb0; RBb.m_qt.RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*m_alpha + zbp*(1-m_alpha);
    
	c = rb + zb - ra - za;
    
	m_F = m_L + c*m_eps;
}

//-----------------------------------------------------------------------------
void FERigidSphericalJoint::Reset()
{
	m_F = vec3d(0,0,0);
	m_L = vec3d(0,0,0);
    
	FERigidBody& RBa = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBa));
	FERigidBody& RBb = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(m_nRBb));
    
	m_qa0 = m_q0 - RBa.m_r0;
	m_qb0 = m_q0 - RBb.m_r0;
}
