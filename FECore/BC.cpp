#include "stdafx.h"
#include "BC.h"
#include "FEModel.h"
#include "FESolver.h"
#include "FERigidBody.h"

//-----------------------------------------------------------------------------
FERigidBodyForce::FERigidBodyForce(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
{
	m_bfollow = false;
}

//-----------------------------------------------------------------------------
void FERigidBodyForce::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << GetID() << IsActive();
		ar << ntype << bc << id << lc << sf;
	}
	else
	{
		int nid;
		bool bactive;
		ar >> nid >> bactive;
		ar >> ntype >> bc >> id >> lc >> sf;
		SetID(nid);
		if (bactive) Activate(); else Deactivate();
	}
}

//-----------------------------------------------------------------------------
double FERigidBodyForce::Value()
{
	FEModel& fem = *GetFEModel();
	if (lc >= 0)
		return fem.GetLoadCurve(lc)->Value()*sf;
	else
		return sf;
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidBodyForce::Residual(FEGlobalVector& R, FETimePoint& tp)
{
	FEModel& fem = *GetFEModel();
	FERigidBody& rb = static_cast<FERigidBody&>(*fem.Object(id));

	// setup the force vector
	vec3d f(0,0,0);
	if      (bc == 0) f.x = Value();
	else if (bc == 1) f.y = Value();
	else if (bc == 2) f.z = Value();
	
	// for follower forces apply the rigid body rotation
	if (m_bfollow) f = rb.m_qt*f;

	// add to the residual
	int n;
	n = rb.m_LM[0]; if (n >= 0) R[n] += f.x;
	n = rb.m_LM[1]; if (n >= 0) R[n] += f.y;
	n = rb.m_LM[2]; if (n >= 0) R[n] += f.z;
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
void FERigidBodyForce::StiffnessMatrix(FESolver* psolver, FETimePoint& tp)
{
	// I think for follower forces, I need to contribute to the stiffness matrix, but I'm not sure yet what.
}

//-----------------------------------------------------------------------------
double FERigidBodyDisplacement::Value()
{
	FEModel& fem = *GetFEModel();
	if (lc < 0) return 0;
	else return sf*fem.GetLoadCurve(lc)->Value() + ref;
}

//=============================================================================
FERigidAxialForce::FERigidAxialForce(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
{
	m_ida = m_idb = -1;
	m_ra0 = m_rb0 = vec3d(0,0,0);
	m_s = 0.0;
	m_lc = -1;
}

//-----------------------------------------------------------------------------
void FERigidAxialForce::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << GetID() << IsActive();
	}
	else
	{
		int nid;
		bool bactive;
		ar >> nid >> bactive;
		SetID(nid);
		if (bactive) Activate(); else Deactivate();
	}
}

//-----------------------------------------------------------------------------
double FERigidAxialForce::Value()
{
	FEModel& fem = *GetFEModel();
	if (m_lc > 0)
		return fem.GetLoadCurve(m_lc - 1)->Value()*m_s;
	else
		return m_s;
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidAxialForce::Residual(FEGlobalVector& R, FETimePoint& tp)
{
	FEModel& fem = *GetFEModel();
	FERigidBody& bodyA = static_cast<FERigidBody&>(*fem.Object(m_ida));
	FERigidBody& bodyB = static_cast<FERigidBody&>(*fem.Object(m_idb));

	// get the attachment position in global coordinates for body A
	vec3d da0 = m_ra0 - bodyA.m_r0;
	vec3d da = bodyA.m_qt*da0;
	vec3d a = da + bodyA.m_rt;

	// get the attachment position in global coordinates for body B
	vec3d db0 = m_rb0 - bodyB.m_r0;
	vec3d db = bodyB.m_qt*db0;
	vec3d b = db + bodyB.m_rt;

	// get the unit axial vector
	vec3d N = b - a; N.unit();

	// calculate the force value
	double f = Value();

	// get the axial force and torques
	vec3d F = N*f;
	vec3d Ma = F^da;
	vec3d Mb = F^db;

	// apply force and torque to body A
	int n;
	n = bodyA.m_LM[0]; if (n >= 0) R[n] += F.x;
	n = bodyA.m_LM[1]; if (n >= 0) R[n] += F.y;
	n = bodyA.m_LM[2]; if (n >= 0) R[n] += F.z;
	n = bodyA.m_LM[3]; if (n >= 0) R[n] += Ma.x;
	n = bodyA.m_LM[4]; if (n >= 0) R[n] += Ma.y;
	n = bodyA.m_LM[5]; if (n >= 0) R[n] += Ma.z;

	// apply force and torque to body B
	n = bodyB.m_LM[0]; if (n >= 0) R[n] -= F.x;
	n = bodyB.m_LM[1]; if (n >= 0) R[n] -= F.y;
	n = bodyB.m_LM[2]; if (n >= 0) R[n] -= F.z;
	n = bodyB.m_LM[3]; if (n >= 0) R[n] -= Mb.x;
	n = bodyB.m_LM[4]; if (n >= 0) R[n] -= Mb.y;
	n = bodyB.m_LM[5]; if (n >= 0) R[n] -= Mb.z;
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
//! TODO: Only the stiffness contribution in the were the axial forces are applied
//!       to the center of mass has been implemented. 
void FERigidAxialForce::StiffnessMatrix(FESolver* psolver, const FETimePoint& tp)
{
	// Get the rigid bodies
	FEModel& fem = *GetFEModel();
	FERigidBody& bodyA = static_cast<FERigidBody&>(*fem.Object(m_ida));
	FERigidBody& bodyB = static_cast<FERigidBody&>(*fem.Object(m_idb));

	// get the attachment position in global coordinates for body A
	vec3d da0 = m_ra0 - bodyA.m_r0;
	vec3d da = bodyA.m_qt*da0;
	vec3d pa = da + bodyA.m_rt;

	// get the attachment position in global coordinates for body B
	vec3d db0 = m_rb0 - bodyB.m_r0;
	vec3d db = bodyB.m_qt*db0;
	vec3d pb = db + bodyB.m_rt;

	// setup the axial unit vector
	vec3d N = pb - pa; 
	double L = N.unit();
	double Linv = 1.0 / L;

	// calculate the force value
	double f = -Value()*Linv;

	// build the stiffness matrix
	double k[3][3];
	k[0][0] = f*(-1.0 + N.x*N.x); k[0][1] = f*N.x*N.y; k[0][2] = f*N.x*N.z;
	k[1][1] = f*(-1.0 + N.y*N.y); k[1][0] = f*N.y*N.x; k[1][2] = f*N.y*N.z;
	k[2][2] = f*(-1.0 + N.z*N.z); k[2][0] = f*N.z*N.x; k[2][1] = f*N.z*N.y;

	matrix K(6,6);
	K[0][0] = k[0][0]; K[0][1] = k[0][1]; K[0][2] = k[0][2]; K[0][3] = -k[0][0]; K[0][4] = -k[0][1]; K[0][5] = -k[0][2];
	K[1][0] = k[1][0]; K[1][1] = k[1][1]; K[1][2] = k[1][2]; K[1][3] = -k[1][0]; K[1][4] = -k[1][1]; K[1][5] = -k[1][2];
	K[2][0] = k[2][0]; K[2][1] = k[2][1]; K[2][2] = k[2][2]; K[2][3] = -k[2][0]; K[2][4] = -k[2][1]; K[2][5] = -k[2][2];

	K[3][0] = -k[0][0]; K[3][1] = -k[0][1]; K[3][2] = -k[0][2]; K[3][3] = k[0][0]; K[3][4] = k[0][1]; K[3][5] = k[0][2];
	K[4][0] = -k[1][0]; K[4][1] = -k[1][1]; K[4][2] = -k[1][2]; K[4][3] = k[1][0]; K[4][4] = k[1][1]; K[4][5] = k[1][2];
	K[5][0] = -k[2][0]; K[5][1] = -k[2][1]; K[5][2] = -k[2][2]; K[5][3] = k[2][0]; K[5][4] = k[2][1]; K[5][5] = k[2][2];

	// get the equation numbers
	vector<int> lm(6);
	lm[0] = bodyA.m_LM[0];
	lm[1] = bodyA.m_LM[1];
	lm[2] = bodyA.m_LM[2];
	lm[3] = bodyB.m_LM[0];
	lm[4] = bodyB.m_LM[1];
	lm[5] = bodyB.m_LM[2];

	// assemble into global matrix
	psolver->AssembleStiffness(lm, K);
}
