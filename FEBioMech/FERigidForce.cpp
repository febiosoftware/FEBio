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
#include "FERigidForce.h"
#include "FERigidBody.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEMaterial.h"
#include "FECore/FELoadCurve.h"
#include "FEMechModel.h"
#include "FERigidMaterial.h"
#include <FECore/FELinearSystem.h>
#include <FECore/log.h>

//=============================================================================
BEGIN_FECORE_CLASS(FERigidAxialForce, FERigidLoad);
	ADD_PARAMETER(m_ida      , "rbA"     )->setEnums("$(rigid_materials)");
	ADD_PARAMETER(m_idb      , "rbB"     )->setEnums("$(rigid_materials)");
	ADD_PARAMETER(m_ra0      , "ra"      );
	ADD_PARAMETER(m_rb0      , "rb"      );
	ADD_PARAMETER(m_s        , "force"   );
	ADD_PARAMETER(m_brelative, "relative");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidAxialForce::FERigidAxialForce(FEModel* pfem) : FERigidLoad(pfem)
{
	m_ida = m_idb = -1;
	m_ra0 = m_rb0 = vec3d(0,0,0);
	m_s = 0.0;
	m_brelative = false;
}

//-----------------------------------------------------------------------------
//! do some sanity checks
bool FERigidAxialForce::Init()
{
	// At this point the rigid ID's are still associated with the materials.
	// We want to associate them with the rigid objects instead.
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_ida-1));
	if (pm == 0) return false;
	m_ida = pm->GetRigidBodyID(); if (m_ida < 0) return false;
	pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_idb-1));
	if (pm == 0) return false;
	m_idb = pm->GetRigidBodyID(); if (m_idb < 0) return false;

	// get the actual rigid bodies
	FERigidBody& bodyA = *fem.GetRigidBody(m_ida);
	FERigidBody& bodyB = *fem.GetRigidBody(m_idb);

	// get the attachment position in global coordinates for body A
	vec3d da0 = (m_brelative ? m_ra0 : m_ra0 - bodyA.m_r0);
	vec3d da = bodyA.GetRotation()*da0;
	vec3d a = da + bodyA.m_rt;

	// get the attachment position in global coordinates for body B
	vec3d db0 = (m_brelative ? m_rb0 : m_rb0 - bodyB.m_r0);
	vec3d db = bodyB.GetRotation()*db0;
	vec3d b = db + bodyB.m_rt;

	// get the unit axial vector
	// and make sure it is of finite length
	vec3d N = b - a; 
	double L = N.unit();
	if (L < 1e-17) return false;

	// all is well in the world
	return true;
}

//-----------------------------------------------------------------------------
void FERigidAxialForce::Serialize(DumpStream& ar)
{
	FEModelLoad::Serialize(ar);
	ar & m_ida & m_idb;
	ar & m_ra0 & m_rb0;
	ar & m_s & m_brelative;
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidAxialForce::LoadVector(FEGlobalVector& R)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& bodyA = *fem.GetRigidBody(m_ida);
	FERigidBody& bodyB = *fem.GetRigidBody(m_idb);

	// get the attachment position in global coordinates for body A
	vec3d da0 = (m_brelative ? m_ra0 : m_ra0 - bodyA.m_r0);
	vec3d da = bodyA.GetRotation()*da0;
	vec3d a = da + bodyA.m_rt;

	// get the attachment position in global coordinates for body B
	vec3d db0 = (m_brelative ? m_rb0 : m_rb0 - bodyB.m_r0);
	vec3d db = bodyB.GetRotation()*db0;
	vec3d b = db + bodyB.m_rt;

	// get the unit axial vector
	vec3d N = b - a; N.unit();

	// calculate the force value
	double f = m_s;

	// get the axial force and torques
	vec3d F = N*f;
	vec3d Ma = da^F;
	vec3d Mb = db^F;

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
    
    bodyA.m_Fr += F;
    bodyA.m_Mr += Ma;
    bodyB.m_Fr -= F;
    bodyB.m_Mr -= Mb;
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
//! TODO: Only the stiffness contribution in the were the axial forces are applied
//!       to the center of mass has been implemented. 
void FERigidAxialForce::StiffnessMatrix(FELinearSystem& LS)
{
	// Get the rigid bodies
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& bodyA = *fem.GetRigidBody(m_ida);
	FERigidBody& bodyB = *fem.GetRigidBody(m_idb);

	// get the attachment position in global coordinates for body A
	vec3d da0 = (m_brelative ? m_ra0 : m_ra0 - bodyA.m_r0);
	vec3d da = bodyA.GetRotation()*da0;
	vec3d pa = da + bodyA.m_rt;

	// get the attachment position in global coordinates for body B
	vec3d db0 = (m_brelative ? m_rb0 : m_rb0 - bodyB.m_r0);
	vec3d db = bodyB.GetRotation()*db0;
	vec3d pb = db + bodyB.m_rt;

	// setup the axial unit vector
	vec3d N = pb - pa; 
	double L = N.unit();

	// calculate the force value
	double f = -m_s / L;

	// build the stiffness matrix components
	mat3ds M = mat3dd(1.0) - dyad(N);
	mat3da A(-da);
	mat3da B(-db);
	mat3da S(-N);

	mat3d MA = M*A;
	mat3d MB = M*B;
	mat3d AM = A*M;

	mat3d AMA = A*MA;
	mat3d BMB = B*MB;
	mat3d AMB = A*MB;

	mat3d SA = S*A;
	mat3d SB = S*B;

	// put it all together
	FEElementMatrix K; K.resize(12, 12); K.zero();
	K.sub(0,0, M);
	K.sub(0,3, MA);
	K.add(0,6, M);
	K.add(0,9, MB);
	K.add(3,3, SA*L + AMA);
	K.add(3,6, AM);
	K.sub(3,9, AMB);
	K.sub(6,6,M);
	K.sub(6,9, MB);
	K.add(9,9, SB*(-L) + BMB);

	// since this is a symmetric matrix, fill the bottom triangular part
	K.copy_ut();

	// and multiply by f
	K *= f;

	// get the equation numbers
	vector<int> lm(12);
	lm[ 0] = bodyA.m_LM[0];
	lm[ 1] = bodyA.m_LM[1];
	lm[ 2] = bodyA.m_LM[2];
	lm[ 3] = bodyA.m_LM[3];
	lm[ 4] = bodyA.m_LM[4];
	lm[ 5] = bodyA.m_LM[5];

	lm[ 6] = bodyB.m_LM[0];
	lm[ 7] = bodyB.m_LM[1];
	lm[ 8] = bodyB.m_LM[2];
	lm[ 9] = bodyB.m_LM[3];
	lm[10] = bodyB.m_LM[4];
	lm[11] = bodyB.m_LM[5];

	// assemble into global matrix
	K.SetIndices(lm);
	LS.Assemble(K);
}

//=============================================================================

BEGIN_FECORE_CLASS(FERigidBodyForce, FERigidLoad)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
	ADD_PARAMETER(m_dof, "dof", 0, "Rx\0Ry\0Rz");
	ADD_PARAMETER(m_force, "value")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE)->setUnits(UNIT_FORCE);
	ADD_PARAMETER(m_ntype, "load_type")->setEnums("LOAD\0FOLLOW\0TARGET");
	ADD_PARAMETER(m_brelative, "relative");
END_FECORE_CLASS();

FERigidBodyForce::FERigidBodyForce(FEModel* pfem) : FERigidLoad(pfem)
{
	m_rigidMat = -1;
	m_ntype = FORCE_LOAD;
	m_rid = -1;
	m_brelative = false;
	m_force0 = 0.0;
	m_force = 0.0;
	m_dof = 0;
}

void FERigidBodyForce::SetRigidMaterialID(int nid) { m_rigidMat = nid; }

void FERigidBodyForce::SetDOF(int bc) { m_dof = bc; }

void FERigidBodyForce::SetLoadType(int loadType) { m_ntype = loadType; }

void FERigidBodyForce::SetForce(double f) { m_force = f; }

//-----------------------------------------------------------------------------
bool FERigidBodyForce::Init()
{
	// At this point the rigid ID's are still associated with the materials.
	// We want to associate them with the rigid objects instead.
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	if (pm == 0) return false;

	// check the dof value
	if ((m_dof < 0) || (m_dof >= 3)) return false;

	return FEModelLoad::Init();
}

//-----------------------------------------------------------------------------
void FERigidBodyForce::Activate()
{
	FEModelLoad::Activate();

	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));

	m_rid = pm->GetRigidBodyID(); assert(m_rid >= 0);

	FERigidBody& rb = *fem.GetRigidBody(m_rid);
	if (m_ntype == FORCE_TARGET)
	{
		switch (m_dof)
		{
		case 0: m_force0 = rb.m_Fr.x; break;
		case 1: m_force0 = rb.m_Fr.y; break;
		case 2: m_force0 = rb.m_Fr.z; break;
		}
	}
	if ((m_ntype == FORCE_LOAD) && m_brelative)
	{
		switch (m_dof)
		{
		case 0: m_force0 = rb.m_Fr.x; break;
		case 1: m_force0 = rb.m_Fr.y; break;
		case 2: m_force0 = rb.m_Fr.z; break;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyForce::Serialize(DumpStream& ar)
{
	FEModelLoad::Serialize(ar);
	ar & m_ntype & m_dof & m_rigidMat & m_force & m_rid;
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidBodyForce::LoadVector(FEGlobalVector& R)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(m_rid);
	double t = CurrentTime();

	if (m_ntype == FORCE_FOLLOW)
	{
		// setup the force vector
		vec3d f(0,0,0);
		if      (m_dof == 0) f.x = m_force;
		else if (m_dof == 1) f.y = m_force;
		else if (m_dof == 2) f.z = m_force;
	
		// apply the rigid body rotation
		f = rb.GetRotation()*f;

		// add to the residual
		int n;
		n = rb.m_LM[0]; if (n >= 0) R[n] += f.x;
		n = rb.m_LM[1]; if (n >= 0) R[n] += f.y;
		n = rb.m_LM[2]; if (n >= 0) R[n] += f.z;
	}
	else if (m_ntype == FORCE_LOAD)
	{
		int I = rb.m_LM[m_dof];
		if (I >= 0)
		{
			R[I] += m_force + m_force0;
		}
	}
	else if (m_ntype == FORCE_TARGET)
	{
		// get the current analysis step
		FEAnalysis* pstep = fem.GetCurrentStep();

		double t0 = pstep->m_tstart;
		double t1 = pstep->m_tend;
		double w = (t - t0) / (t1 - t0);
		assert((w >= -0.0000001) && (w <= 1.0000001));
		double f0 = m_force0, f1 = m_force;

		double f = f0 * (1.0 - w) + f1 * w;

		int I = rb.m_LM[m_dof];
		R[I] += f;
	}
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
void FERigidBodyForce::StiffnessMatrix(FELinearSystem& LS)
{
	// I think for follower forces I need to contribute to the stiffness matrix, but I'm not sure yet what.
}


//=============================================================================

BEGIN_FECORE_CLASS(FERigidBodyMoment, FERigidLoad)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
	ADD_PARAMETER(m_dof, "dof", 0, "Ru\0Rv\0Rw\0");
	ADD_PARAMETER(m_value, "value")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE)->setUnits(UNIT_MOMENT);
	ADD_PARAMETER(m_brelative, "relative");
	ADD_PARAMETER(m_ntype, "load_type")->setEnums("LOAD\0FOLLOW\0TARGET");
END_FECORE_CLASS();

FERigidBodyMoment::FERigidBodyMoment(FEModel* pfem) : FERigidLoad(pfem)
{
	m_rigidMat = -1;
	m_rid = -1;
	m_brelative = false;
	m_value0 = 0.0;
	m_value = 0.0;
	m_dof = 0;
	m_ntype = MOMENT_LOAD;
}

void FERigidBodyMoment::SetRigidMaterialID(int nid) { m_rigidMat = nid; }

void FERigidBodyMoment::SetDOF(int bc) { m_dof = bc; }

void FERigidBodyMoment::SetValue(double f) { m_value = f; }

//-----------------------------------------------------------------------------
bool FERigidBodyMoment::Init()
{
	// At this point the rigid ID's are still associated with the materials.
	// We want to associate them with the rigid objects instead.
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	if (pm == 0) return false;

	// follower moments are not supported yet.
	if (m_ntype == MOMENT_FOLLOW)
	{
		feLogError("Follower moments are not supported.");
		return false;
	}

	return FEModelLoad::Init();
}

//-----------------------------------------------------------------------------
void FERigidBodyMoment::Activate()
{
	FEModelLoad::Activate();

	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));

	m_rid = pm->GetRigidBodyID(); assert(m_rid >= 0);

	FERigidBody& rb = *fem.GetRigidBody(m_rid);
	if (m_ntype == MOMENT_TARGET)
	{
		switch (m_dof)
		{
		case 0: m_value0 = rb.m_Mr.x; break;
		case 1: m_value0 = rb.m_Mr.y; break;
		case 2: m_value0 = rb.m_Mr.z; break;
		}
	}
	if ((m_ntype == MOMENT_LOAD) && m_brelative)
	{
		switch (m_dof)
		{
		case 0: m_value0 = rb.m_Mr.x; break;
		case 1: m_value0 = rb.m_Mr.y; break;
		case 2: m_value0 = rb.m_Mr.z; break;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyMoment::Serialize(DumpStream& ar)
{
	FEModelLoad::Serialize(ar);
	ar & m_value0 & m_rid;
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidBodyMoment::LoadVector(FEGlobalVector& R)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(m_rid);
	double t = CurrentTime();

	int I = rb.m_LM[m_dof + 3];
	if (I >= 0)
	{
		if (m_ntype == MOMENT_LOAD)
		{
			R[I] += m_value + m_value0;
		}
		else if (m_ntype == MOMENT_TARGET)
		{
			// get the current analysis step
			FEAnalysis* pstep = fem.GetCurrentStep();

			double t0 = pstep->m_tstart;
			double t1 = pstep->m_tend;
			double w = (t - t0) / (t1 - t0);
			assert((w >= -0.0000001) && (w <= 1.0000001));
			double M0 = m_value0, M1 = m_value;

			double M = M0 * (1.0 - w) + M1 * w;
			R[I] += M;
		}
	}
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
void FERigidBodyMoment::StiffnessMatrix(FELinearSystem& LS)
{

}
