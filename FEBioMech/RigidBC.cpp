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
#include "RigidBC.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include "FERigidBody.h"
#include <FECore/FEMaterial.h>
#include <FECore/FELoadCurve.h>
#include "FEMechModel.h"
#include "FERigidMaterial.h"
#include <FECore/DumpStream.h>
#include <FECore/log.h>

REGISTER_SUPER_CLASS(FERigidBC, FERIGIDBC_ID);

//=============================================================================
BEGIN_FECORE_CLASS(FERigidBC, FEStepComponent)
	// NOTE: This parameter is hidden, since FEBio Studio implements its own mechanism
	//       for assigning the rigid material ID.
	ADD_PARAMETER(m_rigidMat, "rb")->SetFlags(FE_PARAM_HIDDEN);
END_FECORE_CLASS();

FERigidBC::FERigidBC(FEModel* fem) : FEStepComponent(fem)
{
	m_rigidMat = -1;
	m_rb = -1;
}

bool FERigidBC::Init()
{
	// Make sure the rigid material ID is valid
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	if (pm == nullptr) return false;

	return FEStepComponent::Init();
}

void FERigidBC::Activate()
{
	// Get the Rigidbody ID
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1)); assert(pm);
	m_rb = pm->GetRigidBodyID(); assert(m_rb >= 0);

	FEStepComponent::Activate();
}

void FERigidBC::Serialize(DumpStream& ar)
{
	ar & m_rb;
	FEStepComponent::Serialize(ar);
}

FERigidBody& FERigidBC::GetRigidBody()
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& RB = *fem.GetRigidBody(m_rb);
	return RB;
}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidBodyFixedBC, FERigidBC)
	ADD_PARAMETER(m_dofs, "dofs", 0, "Rx\0Ry\0Rz\0Ru\0Rv\0Rw\0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidBodyFixedBC::FERigidBodyFixedBC(FEModel* pfem) : FERigidBC(pfem)
{
	m_binit = false;
}

//-----------------------------------------------------------------------------
bool FERigidBodyFixedBC::Init()
{
	// make sure we have a valid dof
	// TODO: Can I test this automatically during validation?
	for (int i = 0; i < m_dofs.size(); ++i)
	{
		int dof_i = m_dofs[i];
		if ((dof_i < 0) || (dof_i >= 6))
		{
			feLogError("Invalid value for dofs of fixed rigid constraint %s", GetName().c_str());
			return false;
		}
	}

	m_binit = true;
	return FERigidBC::Init();
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Activate()
{
	FERigidBC::Activate();

	if (m_binit == false) Init();
	if (m_binit)
	{
		FERigidBody& RB = GetRigidBody();

		// we only fix the open dofs. If a user accidentally applied a fixed and prescribed
		// rigid degree of freedom, then we make sure the prescribed takes precedence.
		for (int i = 0; i < m_dofs.size(); ++i)
		{
			int dof_i = m_dofs[i];
			if (RB.m_BC[dof_i] == DOF_OPEN) RB.m_BC[dof_i] = DOF_FIXED;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Deactivate()
{
	if (m_binit)
	{
		FERigidBody& RB = GetRigidBody();

		// Since fixed rigid dofs can be overwritten by prescribed dofs, 
		// we have to make sure that this dof is actually a fixed dof.
		for (int i = 0; i < m_dofs.size(); ++i)
		{
			int dof_i = m_dofs[i];
			if (RB.m_BC[dof_i] == DOF_FIXED) RB.m_BC[dof_i] = DOF_OPEN;
		}
	}

	FERigidBC::Deactivate();
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Serialize(DumpStream& ar)
{
	FERigidBC::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dofs & m_binit;
}

//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FERigidBodyDisplacement, FERigidBC)
	ADD_PARAMETER(m_dof, "dof", 0, "Rx\0Ry\0Rz\0Ru\0Rv\0Rw\0");
	ADD_PARAMETER(m_val, "value");
	ADD_PARAMETER(m_brel, "relative");
END_FECORE_CLASS();

FERigidBodyDisplacement::FERigidBodyDisplacement(FEModel* pfem) : FERigidBC(pfem)
{
	m_dof = -1;
	m_val = 0.0;
	m_ref= 0.0; 
	m_brel = false; 
	m_binit = false;
}

//-----------------------------------------------------------------------------
bool FERigidBodyDisplacement::Init()
{
	// make sure we have a valid dof
	if ((m_dof < 0)||(m_dof >=6)) return false;
	m_binit = true;
	return FERigidBC::Init();
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Activate()
{
	// don't forget to call the base class
	FERigidBC::Activate();

	if (m_binit == false) Init();

	// get the rigid body
	FERigidBody& RB = GetRigidBody();

	// set some stuff
	RB.m_pDC[m_dof] = this;

	// mark the dof as prescribed
	RB.m_BC[m_dof] = DOF_PRESCRIBED;

	// set the relative offset
	m_ref = 0.0;
	if (m_brel)
	{
		quatd Q = RB.GetRotation();
		vec3d q = Q.GetRotationVector();
		switch (m_dof)
		{
		case 0: m_ref = RB.m_rt.x - RB.m_r0.x; break;
		case 1: m_ref = RB.m_rt.y - RB.m_r0.y; break;
		case 2: m_ref = RB.m_rt.z - RB.m_r0.z; break;
		case 3: m_ref = q.x; break;
		case 4: m_ref = q.y; break;
		case 5: m_ref = q.z; break;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Deactivate()
{
	FERigidBC::Deactivate();

	// get the rigid body
	// Since Deactivate is called before Init (for multi-step analysis; in the FEBio input)
	// we have to make sure the data is initialized
	if (m_binit)
	{
		FERigidBody& RB = GetRigidBody();

		// turn off the prescribed displacement
		RB.m_pDC[m_dof] = 0;
		RB.m_BC[m_dof] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Serialize(DumpStream& ar)
{
	FERigidBC::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof & m_ref & m_binit & m_brel;
}

//-----------------------------------------------------------------------------
double FERigidBodyDisplacement::Value()
{
	return m_val + m_ref;
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::InitTimeStep()
{
	FERigidBody& RB = GetRigidBody();
	int I = GetBC();
	RB.m_dul[I] = Value() - RB.m_Ut[I];
}

//=============================================================================
FERigidIC::FERigidIC(FEModel* fem) : FERigidBC(fem)
{

}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidBodyVelocity, FERigidBC)
	ADD_PARAMETER(m_vel, "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidBodyVelocity::FERigidBodyVelocity(FEModel* pfem) : FERigidIC(pfem) 
{
	m_vel = vec3d(0, 0, 0);
}

//-----------------------------------------------------------------------------
void FERigidBodyVelocity::Activate()
{
	FERigidIC::Activate();
	FERigidBody& RB = GetRigidBody();
	RB.m_vp = RB.m_vt = m_vel;
}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidBodyAngularVelocity, FERigidBC)
	ADD_PARAMETER(m_w, "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidBodyAngularVelocity::FERigidBodyAngularVelocity(FEModel* pfem) : FERigidIC(pfem)
{
	m_w = vec3d(0, 0, 0);
}

//-----------------------------------------------------------------------------
void FERigidBodyAngularVelocity::Activate()
{
	FERigidIC::Activate();
	FERigidBody& RB = GetRigidBody();
	RB.m_wp = RB.m_wt = m_w;
}
