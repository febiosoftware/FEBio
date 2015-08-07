#include "stdafx.h"
#include "FESolidAnalysis.h"
#include "FERigidMaterial.h"
#include "FEContactInterface.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FESolidAnalysis::FESolidAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_SOLID) 
{

}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FESolidAnalysis::Activate()
{
	// initialize base class data
	FEAnalysis::Activate();

	// reset nodal ID's
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		// fix all degrees of freedom
		FENode& node = mesh.Node(i);
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = DOF_FIXED;

		// open degrees of freedom
		if (node.m_bexclude == false)
		{
			// turn on displacement dofs (for non-rigid nodes)
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = DOF_OPEN;
				node.m_ID[DOF_Y] = DOF_OPEN;
				node.m_ID[DOF_Z] = DOF_OPEN;
			}

			// turn on rotation dofs for non-rigid shell nodes
			if (node.m_bshell)
			{
				node.m_ID[DOF_U] = DOF_OPEN;
				node.m_ID[DOF_V] = DOF_OPEN;
				node.m_ID[DOF_W] = DOF_OPEN;
			}
		}
	}

	// apply fixed dofs
	for (int i=0; i<m_fem.FixedBCs(); ++i)
	{
		FEFixedBC& bc = *m_fem.FixedBC(i);
		bc.Activate();
	}

	// apply prescribed dofs
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		if (DC.IsActive())
		{
			FENode& node = m_fem.GetMesh().Node(DC.node);
			node.m_ID[DC.bc] = DOF_PRESCRIBED;
		}
	}

	// initialize equations
	// ----->
	if (m_psolver->InitEquations() == false) return false;

	// initialize linear constraints
	// Must be done after equations are initialized
	if (InitLinearConstraints() == false) return false;
	// ----->

	// evaluate the relative prescribed values
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		int nid = DC.node;
		int bc  = DC.bc;
		bool br = DC.br;

		FENode& node = m_fem.GetMesh().Node(nid); 

		if (DC.IsActive())
		{
			switch (bc)
			{
			case DOF_X: DC.r = br ? node.m_rt.x - node.m_r0.x : 0; break;
			case DOF_Y: DC.r = br ? node.m_rt.y - node.m_r0.y : 0; break;
			case DOF_Z: DC.r = br ? node.m_rt.z - node.m_r0.z : 0; break;
			case DOF_U: DC.r = br ? node.m_Dt.x - node.m_D0.x : 0; break;
			case DOF_V: DC.r = br ? node.m_Dt.y - node.m_D0.y : 0; break;
			case DOF_W: DC.r = br ? node.m_Dt.z - node.m_D0.z : 0; break;
			}
		}
	}

	// activate the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il) il->Activate();
	}

	if (m_nanalysis == FE_DYNAMIC)
	{
		FERigidSystem& rigid = *m_fem.GetRigidSystem();

		// set the initial velocities of all rigid nodes
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& n = mesh.Node(i);
			if (n.m_rid >= 0)
			{
				FERigidBody& rb = *rigid.Object(n.m_rid);
				vec3d V = rb.m_vt;
				vec3d W = rb.m_wt;
				vec3d r = n.m_rt - rb.m_rt;

				vec3d v = V + (W ^ r); 
				n.m_v0 = n.m_vp = n.m_vt = v;

				vec3d a = (W ^ V)*2.0 + (W ^ (W ^ r));
				n.m_ap = n.m_at = a;
			}
		}
	}

	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		felog.printbox("FATAL ERROR","Failed to initialize solver.\nAborting run.\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FESolidAnalysis::Deactivate()
{
	FEAnalysis::Deactivate();

	// store the current rigid body reaction forces
	FERigidSystem& rigid = *m_fem.GetRigidSystem();
	for (int i=0; i<rigid.Objects(); ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		RB.m_Fp = RB.m_Fr;
		RB.m_Mp = RB.m_Mr;
	}
}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FEExplicitSolidAnalysis::Activate()
{
	// initialize base class data
	FEAnalysis::Activate();

	// reset nodal ID's
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = DOF_FIXED;

		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = 0;
				node.m_ID[DOF_Y] = 0;
				node.m_ID[DOF_Z] = 0;
			}

			if (node.m_bshell)
			{
				node.m_ID[DOF_U] = 0;
				node.m_ID[DOF_V] = 0;
				node.m_ID[DOF_W] = 0;
			}
		}
	}

	// apply fixed dofs
	for (int i=0; i<m_fem.FixedBCs(); ++i)
	{
		FEFixedBC& bc = *m_fem.FixedBC(i);
		bc.Activate();
	}

	// apply prescribed dofs
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		if (DC.IsActive())
		{
			FENode& node = m_fem.GetMesh().Node(DC.node);
			node.m_ID[DC.bc] = DOF_PRESCRIBED;
		}
	}

	// initialize equations
	// ----->
	if (m_psolver->InitEquations() == false) return false;

	// initialize linear constraints
	// Must be done after equations are initialized
	if (InitLinearConstraints() == false) return false;
	// ----->

	// Now we adjust the equation numbers of prescribed dofs according to the above rule
	// Make sure that a prescribed dof has not been fixed
	// TODO: maybe this can be moved to the FESolver::InitEquations function
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		int nid = DC.node;
		int bc  = DC.bc;
		bool br = DC.br;

		FENode& node = m_fem.GetMesh().Node(nid); 

		if (DC.IsActive())
		{
			switch (bc)
			{
			case DOF_X: DC.r = br ? node.m_rt.x - node.m_r0.x : 0; break;
			case DOF_Y: DC.r = br ? node.m_rt.y - node.m_r0.y : 0; break;
			case DOF_Z: DC.r = br ? node.m_rt.z - node.m_r0.z : 0; break;
			case DOF_U: DC.r = br ? node.m_Dt.x - node.m_D0.x : 0; break;
			case DOF_V: DC.r = br ? node.m_Dt.y - node.m_D0.y : 0; break;
			case DOF_W: DC.r = br ? node.m_Dt.z - node.m_D0.z : 0; break;
			}
		}
	}

	// activate the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il) il->Activate();
	}

	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		felog.printbox("FATAL ERROR","Failed to initialize solver.\nAborting run.\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FELinearSolidAnalysis::Activate()
{
	// initialize base class data
	FEAnalysis::Activate();

	// reset nodal ID's
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = -1;

		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = 0;
				node.m_ID[DOF_Y] = 0;
				node.m_ID[DOF_Z] = 0;
			}

			if (node.m_bshell)
			{
				node.m_ID[DOF_U] = 0;
				node.m_ID[DOF_V] = 0;
				node.m_ID[DOF_W] = 0;
			}
		}
	}

	// apply fixed dofs
	for (int i=0; i<m_fem.FixedBCs(); ++i)
	{
		FEFixedBC& bc = *m_fem.FixedBC(i);
		bc.Activate();
	}

	// apply prescribed dofs
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		if (DC.IsActive())
		{
			FENode& node = m_fem.GetMesh().Node(DC.node);
			node.m_ID[DC.bc] = DOF_PRESCRIBED;
		}
	}

	// initialize equations
	// ----->
	if (m_psolver->InitEquations() == false) return false;

	// initialize linear constraints
	// Must be done after equations are initialized
	if (InitLinearConstraints() == false) return false;
	// ----->

	// Now we adjust the equation numbers of prescribed dofs according to the above rule
	// Make sure that a prescribed dof has not been fixed
	// TODO: maybe this can be moved to the FESolver::InitEquations function
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		int nid = DC.node;
		int bc  = DC.bc;
		bool br = DC.br;

		FENode& node = m_fem.GetMesh().Node(nid); 

		if (DC.IsActive())
		{
			switch (bc)
			{
			case DOF_X: DC.r = br ? node.m_rt.x - node.m_r0.x : 0; break;
			case DOF_Y: DC.r = br ? node.m_rt.y - node.m_r0.y : 0; break;
			case DOF_Z: DC.r = br ? node.m_rt.z - node.m_r0.z : 0; break;
			case DOF_U: DC.r = br ? node.m_Dt.x - node.m_D0.x : 0; break;
			case DOF_V: DC.r = br ? node.m_Dt.y - node.m_D0.y : 0; break;
			case DOF_W: DC.r = br ? node.m_Dt.z - node.m_D0.z : 0; break;
			}
		}
	}

	// activate the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il) il->Activate();
	}

	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		felog.printbox("FATAL ERROR","Failed to initialize solver.\nAborting run.\n");
		return false;
	}

	return true;
}
