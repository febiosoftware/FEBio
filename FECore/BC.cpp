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
