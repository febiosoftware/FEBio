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

//=============================================================================
FERigidBC::FERigidBC(FEModel* fem) : FEBoundaryCondition(fem)
{
	m_rigidMat = -1;
	m_rb = -1;

	m_binit = false;
}

bool FERigidBC::Init()
{
	// Make sure the rigid material ID is valid
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	if (pm == nullptr) return false;

	m_rb = pm->GetRigidBodyID(); assert(m_rb >= 0);

	m_binit = true;

	return FEBoundaryCondition::Init();
}

void FERigidBC::CopyFrom(FEBoundaryCondition* pbc)
{
	FERigidBC* rbc = dynamic_cast<FERigidBC*>(pbc);
	GetParameterList() = rbc->GetParameterList();
	m_rb = rbc->m_rb;
}

void FERigidBC::Serialize(DumpStream& ar)
{
	ar & m_rb & m_binit;
	FEBoundaryCondition::Serialize(ar);
}

FERigidBody& FERigidBC::GetRigidBody()
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& RB = *fem.GetRigidBody(m_rb);
	return RB;
}

//=============================================================================
FERigidFixedBC::FERigidFixedBC(FEModel* pfem) : FERigidBC(pfem)
{
}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidFixedBCNew, FERigidFixedBC)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
	ADD_PARAMETER(m_dof[0], "Rx_dof")->setLongName("X-displacement");
	ADD_PARAMETER(m_dof[1], "Ry_dof")->setLongName("Y-displacement");
	ADD_PARAMETER(m_dof[2], "Rz_dof")->setLongName("Z-displacement");
	ADD_PARAMETER(m_dof[3], "Ru_dof")->setLongName("X-rotation");
	ADD_PARAMETER(m_dof[4], "Rv_dof")->setLongName("Y-rotation");
	ADD_PARAMETER(m_dof[5], "Rw_dof")->setLongName("Z-rotation");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidFixedBCNew::FERigidFixedBCNew(FEModel* pfem) : FERigidFixedBC(pfem)
{
	for (int i = 0; i < 6; ++i) m_dof[i] = false;
	m_binit = false;
}

//-----------------------------------------------------------------------------
void FERigidFixedBCNew::Activate()
{
	FERigidFixedBC::Activate();

	if (m_binit == false) Init();
	if (m_binit)
	{
		FERigidBody& RB = GetRigidBody();

		// we only fix the open dofs. If a user accidentally applied a fixed and prescribed
		// rigid degree of freedom, then we make sure the prescribed takes precedence.
		for (int i = 0; i < 6; ++i)
		{
			if (m_dof[i] && (RB.m_BC[i] == DOF_OPEN)) RB.m_BC[i] = DOF_FIXED;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidFixedBCNew::Deactivate()
{
	if (m_binit)
	{
		FERigidBody& RB = GetRigidBody();

		// Since fixed rigid dofs can be overwritten by prescribed dofs, 
		// we have to make sure that this dof is actually a fixed dof.
		for (int i = 0; i < 6; ++i)
		{
			if (m_dof[i] && (RB.m_BC[i] == DOF_FIXED)) RB.m_BC[i] = DOF_OPEN;
		}
	}

	FERigidBC::Deactivate();
}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidFixedBCOld, FERigidFixedBC)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
	ADD_PARAMETER(m_dofs, "dofs", 0, "Rx\0Ry\0Rz\0Ru\0Rv\0Rw\0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidFixedBCOld::FERigidFixedBCOld(FEModel* pfem) : FERigidFixedBC(pfem)
{

}

//-----------------------------------------------------------------------------
bool FERigidFixedBCOld::Init()
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

	return FERigidFixedBC::Init();
}

//-----------------------------------------------------------------------------
void FERigidFixedBCOld::Activate()
{
	FERigidFixedBC::Activate();

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
void FERigidFixedBCOld::Deactivate()
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
}

//=============================================================================
FERigidPrescribedBC::FERigidPrescribedBC(FEModel* pfem) : FERigidBC(pfem)
{
	m_dof = -1;
	m_val = 0.0;
	m_ref= 0.0; 
	m_brel = false; 
	m_binit = false;
}

//-----------------------------------------------------------------------------
bool FERigidPrescribedBC::Init()
{
	// make sure we have a valid dof
	if ((m_dof < 0)||(m_dof >=6)) return false;
	m_binit = true;
	return FERigidBC::Init();
}

//-----------------------------------------------------------------------------
void FERigidPrescribedBC::Activate()
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
void FERigidPrescribedBC::Deactivate()
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
void FERigidPrescribedBC::Serialize(DumpStream& ar)
{
	FERigidBC::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof & m_ref & m_binit & m_brel;
}

//-----------------------------------------------------------------------------
double FERigidPrescribedBC::Value()
{
	return m_val + m_ref;
}

//-----------------------------------------------------------------------------
void FERigidPrescribedBC::InitTimeStep()
{
	FERigidBody& RB = GetRigidBody();
	int I = GetBC();
	RB.m_dul[I] = Value() - RB.m_Ut[I];
}

//===============================================================================
BEGIN_FECORE_CLASS(FERigidDisplacement, FERigidPrescribedBC)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
	ADD_PARAMETER(m_bc, "dof", 0, "$(dof_list:displacement)");
	ADD_PARAMETER(m_val, "value")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE)->setUnits(UNIT_LENGTH);
	ADD_PARAMETER(m_brel, "relative");
END_FECORE_CLASS();

FERigidDisplacement::FERigidDisplacement(FEModel* fem) : FERigidPrescribedBC(fem)
{
	m_bc = -1;
}

bool FERigidDisplacement::Init()
{
	int dofX = GetDOFIndex("x");
	int dofY = GetDOFIndex("y");
	int dofZ = GetDOFIndex("z");
	if (m_bc == dofX) m_dof = 0;
	if (m_bc == dofY) m_dof = 1;
	if (m_bc == dofZ) m_dof = 2;

	return FERigidPrescribedBC::Init();
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidRotation, FERigidPrescribedBC)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
	ADD_PARAMETER(m_bc, "dof", 0, "$(dof_list:rigid rotation)");
	ADD_PARAMETER(m_val, "value")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE)->setUnits(UNIT_RADIAN);
	ADD_PARAMETER(m_brel, "relative");
END_FECORE_CLASS();

FERigidRotation::FERigidRotation(FEModel* fem) : FERigidPrescribedBC(fem)
{
	m_bc = -1;
}

bool FERigidRotation::Init()
{
	int dofRu = GetDOFIndex("Ru");
	int dofRv = GetDOFIndex("Rv");
	int dofRw = GetDOFIndex("Rw");
	if (m_bc == dofRu) m_dof = 3;
	if (m_bc == dofRv) m_dof = 4;
	if (m_bc == dofRw) m_dof = 5;

	return FERigidPrescribedBC::Init();
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidPrescribedOld, FERigidPrescribedBC)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
	ADD_PARAMETER(m_dof, "dof", 0, "Rx\0Ry\0Rz\0Ru\0Rv\0Rw\0");
	ADD_PARAMETER(m_val, "value");
	ADD_PARAMETER(m_brel, "relative");
END_FECORE_CLASS();

FERigidPrescribedOld::FERigidPrescribedOld(FEModel* fem) : FERigidPrescribedBC(fem) {}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidIC, FEInitialCondition)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
END_FECORE_CLASS();

FERigidIC::FERigidIC(FEModel* fem) : FEInitialCondition(fem)
{
	m_rigidMat = -1;
	m_rb = -1;
}

bool FERigidIC::Init()
{
	// Make sure the rigid material ID is valid
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	if (pm == nullptr) return false;

	return FEInitialCondition::Init();
}

void FERigidIC::Activate()
{
	// Get the Rigidbody ID
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1)); assert(pm);
	m_rb = pm->GetRigidBodyID(); assert(m_rb >= 0);

	FEInitialCondition::Activate();
}

void FERigidIC::Serialize(DumpStream& ar)
{
	ar& m_rb;
	FEInitialCondition::Serialize(ar);
}

FERigidBody& FERigidIC::GetRigidBody()
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& RB = *fem.GetRigidBody(m_rb);
	return RB;
}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidBodyVelocity, FERigidIC)
	ADD_PARAMETER(m_vel, "value")->setUnits(UNIT_VELOCITY);
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
BEGIN_FECORE_CLASS(FERigidBodyAngularVelocity, FERigidIC)
	ADD_PARAMETER(m_w, "value")->setUnits(UNIT_ANGULAR_VELOCITY);
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
