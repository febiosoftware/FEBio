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
void FERigidAxialForce::StiffnessMatrix(FESolver* psolver, FETimePoint& tp)
{
	// I think for follower forces, I need to contribute to the stiffness matrix, but I'm not sure yet what.
}
