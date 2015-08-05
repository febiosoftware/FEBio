#include "stdafx.h"
#include "FESolidAnalysis.h"
#include "FERigidMaterial.h"
#include "FEContactInterface.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

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

	// clear the active rigid body BC's
	// (don't overwrite prescribed displacements)
	int NRB = m_fem.Objects();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (int j=0; j<6; ++j) if (RB.m_LM[j] != DOF_PRESCRIBED) RB.m_LM[j] = DOF_OPEN;
	}

	// set the fixed rigid body BC's
	// (don't overwrite prescribed displacements)
	for (int i=0;i<(int) m_fem.m_RBC.size(); ++i)
	{
		FERigidBodyFixedBC* pbc = m_fem.m_RBC[i];
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(pbc->id));
		if (pbc->IsActive() && (RB.m_LM[pbc->bc] != DOF_PRESCRIBED)) RB.m_LM[pbc->bc] = DOF_FIXED;
	}

	// reset nodal ID's
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		// fix all degrees of freedom
		FENode& node = mesh.Node(i);
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = -1;

		// open degrees of freedom
		if (node.m_bexclude == false)
		{
			// turn on displacement dofs (for non-rigid nodes)
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = 0;
				node.m_ID[DOF_Y] = 0;
				node.m_ID[DOF_Z] = 0;
			}

			// turn on rotation dofs for non-rigid shell nodes
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

	// override prescribed displacements for rigid nodes
	bool bdisp = false;
	int nbc = m_fem.PrescribedBCs();
	for (int i=0; i<nbc; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive())
		{
			// if the node is not free we don't use this prescribed displacement
			// note that we don't do this for prescribed pressures and concentrations
			if ((dc.bc != DOF_P) && (dc.bc < DOF_C))
			{
				FENode& node = m_fem.GetMesh().Node(dc.node);
				if (node.m_rid >= 0) 
				{
					dc.Deactivate();
					bdisp = true;
				}
			}
		}
	}

	if (bdisp) felog.printbox("WARNING", "Rigid degrees of freedom cannot be prescribed.");

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
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		int nid = DC.node;
		int bc  = DC.bc;
		bool br = DC.br;

		FENode& node = m_fem.GetMesh().Node(nid); 

		if (DC.IsActive())
		{
			int n;
			switch (bc)
			{
			case DOF_X: // x-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.x - node.m_r0.x : 0;
				break;
			case DOF_Y: // y-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.y - node.m_r0.y : 0;
				break;
			case DOF_Z: // z-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.z - node.m_r0.z : 0;
				break;
			case DOF_U: // x-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.x - node.m_D0.x : 0;
				break;
			case DOF_V: // y-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.y - node.m_D0.y : 0;
				break;
			case DOF_W: // z-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.z - node.m_D0.z : 0;
				break;
			case DOF_P: // prescribed pressure
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_pt - node.m_p0 : 0;
				break;
			case DOF_T: // precribed temperature
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = 0;
				break;
/*			case DOF_C: // precribed concentration
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_ct[0] - node.m_c0[0] : 0;
				break;
			case DOF_C+1: // precribed concentration
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_ct[1] - node.m_c0[1] : 0;
				break;*/
				//--> TODO: change bc=20 to something else
			case 20: // y-z radial displacement
				n = node.m_ID[DOF_Y];
				node.m_ID[DOF_Y] = (n<0?n:-n-2);
				n = node.m_ID[DOF_Z];
				node.m_ID[DOF_Z] = (n<0?n:-n-2);
				DC.r = 0;
				break;
			default:	// all prescribed concentrations
				if ((bc >= DOF_C) && (bc < (int)node.m_ID.size())) {
					n = node.m_ID[bc];
					node.m_ID[bc] = (n<0?n:-n-2);
					int sid = bc - DOF_C;
					DC.r = br ? node.m_ct[sid] - node.m_c0[sid] : 0;
				}
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
		// set the initial velocities of all rigid nodes
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& n = mesh.Node(i);
			if (n.m_rid >= 0)
			{
				FERigidBody& rb = static_cast<FERigidBody&>(*m_fem.Object(n.m_rid));
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
	for (int i=0; i<m_fem.Objects(); ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
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

	// clear the active rigid body BC's
	int NRB = m_fem.Objects();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (int j=0; j<6; ++j) if (RB.m_LM[j] != DOF_PRESCRIBED) RB.m_LM[j] = DOF_OPEN;
	}

	// set the fixed rigid body BC's
	for (int i=0;i<(int) m_fem.m_RBC.size(); ++i)
	{
		FERigidBodyFixedBC* pbc = m_fem.m_RBC[i];
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(pbc->id));
		if (pbc->IsActive() && (RB.m_LM[pbc->bc] != DOF_PRESCRIBED)) RB.m_LM[pbc->bc] = DOF_FIXED;
	}

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

	// override prescribed displacements for rigid nodes
	bool bdisp = false;
	int nbc = m_fem.PrescribedBCs();
	for (int i=0; i<nbc; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive())
		{
			// if the node is not free we don't use this prescribed displacement
			// note that we don't do this for prescribed pressures and concentrations
			if ((dc.bc != DOF_P) && (dc.bc < DOF_C))
			{
				FENode& node = m_fem.GetMesh().Node(dc.node);
				if (node.m_rid >= 0) 
				{
					dc.Deactivate();
					bdisp = true;
				}
			}
		}
	}

	if (bdisp) felog.printbox("WARNING", "Rigid degrees of freedom cannot be prescribed.");

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
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		int nid = DC.node;
		int bc  = DC.bc;
		bool br = DC.br;

		FENode& node = m_fem.GetMesh().Node(nid); 

		if (DC.IsActive())
		{
			int n;
			switch (bc)
			{
			case DOF_X: // x-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.x - node.m_r0.x : 0;
				break;
			case DOF_Y: // y-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.y - node.m_r0.y : 0;
				break;
			case DOF_Z: // z-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.z - node.m_r0.z : 0;
				break;
			case DOF_U: // x-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.x - node.m_D0.x : 0;
				break;
			case DOF_V: // y-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.y - node.m_D0.y : 0;
				break;
			case DOF_W: // z-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.z - node.m_D0.z : 0;
				break;
			case DOF_P: // prescribed pressure
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_pt - node.m_p0 : 0;
				break;
			case DOF_T: // precribed temperature
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = 0;
				break;
/*			case DOF_C: // precribed concentration
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_ct[0] - node.m_c0[0] : 0;
				break;
			case DOF_C+1: // precribed concentration
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_ct[1] - node.m_c0[1] : 0;
				break;*/
				//--> TODO: change bc=20 to something else
			case 20: // y-z radial displacement
				n = node.m_ID[DOF_Y];
				node.m_ID[DOF_Y] = (n<0?n:-n-2);
				n = node.m_ID[DOF_Z];
				node.m_ID[DOF_Z] = (n<0?n:-n-2);
				DC.r = 0;
				break;
			default:	// all prescribed concentrations
				if ((bc >= DOF_C) && (bc < (int)node.m_ID.size())) {
					n = node.m_ID[bc];
					node.m_ID[bc] = (n<0?n:-n-2);
					int sid = bc - DOF_C;
					DC.r = br ? node.m_ct[sid] - node.m_c0[sid] : 0;
				}
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
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		int nid = DC.node;
		int bc  = DC.bc;
		bool br = DC.br;

		FENode& node = m_fem.GetMesh().Node(nid); 

		if (DC.IsActive())
		{
			int n;
			switch (bc)
			{
			case DOF_X: // x-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.x - node.m_r0.x : 0;
				break;
			case DOF_Y: // y-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.y - node.m_r0.y : 0;
				break;
			case DOF_Z: // z-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.z - node.m_r0.z : 0;
				break;
			case DOF_U: // x-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.x - node.m_D0.x : 0;
				break;
			case DOF_V: // y-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.y - node.m_D0.y : 0;
				break;
			case DOF_W: // z-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.z - node.m_D0.z : 0;
				break;
			case 20: // y-z radial displacement
				n = node.m_ID[DOF_Y];
				node.m_ID[DOF_Y] = (n<0?n:-n-2);
				n = node.m_ID[DOF_Z];
				node.m_ID[DOF_Z] = (n<0?n:-n-2);
				DC.r = 0;
				break;
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
