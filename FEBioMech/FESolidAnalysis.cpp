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
bool FESolidAnalysis::Init()
{
	// initialize base class data
	FEAnalysis::Init();

	// clear the rigid body BC's
	int NRB = m_fem.Objects();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (int j=0; j<6; ++j)
		{
			RB.m_pDC[j] = 0;
			if (RB.m_BC[j] != 2) RB.m_BC[j] = 0;
		}
	}

	// set the fixed rigid body BC's
	for (int i=0;i<(int) m_fem.m_RBC.size(); ++i)
	{
		FERigidBodyFixedBC* pbc = m_fem.m_RBC[i];
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(pbc->id));
		if (RB.IsActive() && pbc->IsActive()) RB.m_BC[pbc->bc] = -1;
	}

	// set the active rigid bodies BC's
	for (int i=0; i<(int) m_fem.m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *(m_fem.m_RDC[i]);
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(DC.id));
		if (RB.IsActive() && DC.IsActive())
		{
			// TODO: I commented this line out since we can have more rigid materials than rigid bodies
			//       I just probably find a better way to test this (or prevent this from happening)
//			assert(RB.m_BC[DC.bc] == 0);	// make sure we are not overriding a fixed bc
			RB.m_pDC[DC.bc] = &DC;
			RB.m_BC[DC.bc] = 1;
			DC.ref = 0.0;
			if (DC.brel)
			{
				switch (DC.bc)
				{
				case 0: DC.ref = RB.m_rt.x - RB.m_r0.x; break;
				case 1: DC.ref = RB.m_rt.y - RB.m_r0.y; break;
				case 2: DC.ref = RB.m_rt.z - RB.m_r0.z; break;
				}
			}
		}
	}

	// reset nodal ID's
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = -1;

		if (node.m_BC[DOF_X] != -1) node.m_ID[DOF_X] = 0;
		if (node.m_BC[DOF_Y] != -1) node.m_ID[DOF_Y] = 0;
		if (node.m_BC[DOF_Z] != -1) node.m_ID[DOF_Z] = 0;
		if (node.m_BC[DOF_U] != -1) node.m_ID[DOF_U] = 0;
		if (node.m_BC[DOF_V] != -1) node.m_ID[DOF_V] = 0;
		if (node.m_BC[DOF_W] != -1) node.m_ID[DOF_W] = 0;
	}

	// set the rigid nodes
	// Note that also the rotational degrees of freedom are fixed
	// for rigid nodes that do not belong to a non-rigid shell element.
	int nrn = m_fem.RigidNodes();
	for (int i=0; i<nrn; ++i)
	{
		FERigidNode& rn = *m_fem.RigidNode(i);
		if (rn.IsActive())
		{
			FENode& node = m_fem.GetMesh().Node(rn.nid);
			node.m_rid = rn.rid;

			// fix degrees of freedom
			node.m_ID[DOF_X] = -1;
			node.m_ID[DOF_Y] = -1;
			node.m_ID[DOF_Z] = -1;
			if (node.m_bshell == false)
			{
				node.m_ID[DOF_U] = -1;
				node.m_ID[DOF_V] = -1;
				node.m_ID[DOF_W] = -1;
			}
		}
		else 
		{
			FENode& node = m_fem.GetMesh().Node(rn.nid);
			node.m_rid = -1;
		}
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

	// Sometimes an (ignorant) user might have added a rigid body
	// that is not being used. Since this can cause problems we need
	// to find these rigid bodies.
	int nrb = m_fem.Objects();
	vector<int> mec; mec.assign(nrb, 0);
	for (int i=0; i<m_fem.GetMesh().Nodes(); ++i)
	{
		FENode& node = m_fem.GetMesh().Node(i);
		int n = node.m_rid;
		if (n >= 0) mec[n]++;
	}

	for (int i=0; i<nrb; ++i)
		if (mec[i] == 0)
		{
			FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
			felog.printbox("WARNING", "Rigid body %d is not being used.", RB.m_mat+1);
			RB.Activate(false);
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

	// modify the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il)
		{
			list<FELinearConstraint::SlaveDOF>::iterator is = il->slave.begin();
			for (int i=0; i<(int) il->slave.size(); ++i, ++is)
			{
				is->neq = m_fem.GetMesh().Node(is->node).m_ID[is->bc];
			}
		}
	}

	if (m_nanalysis == FE_DYNAMIC)
	{
		// Set the initial velocity for rigid bodies
		for (int i=0; i<(int)m_fem.m_RBV.size(); ++i)
		{
			FERigidBodyVelocity* pv = m_fem.m_RBV[i];
			if (pv->IsActive())
			{
				FERigidBody& rb = static_cast<FERigidBody&>(*m_fem.Object(i));
				rb.m_vp = rb.m_vt = pv->v;
			}
		}

		// Set the initial angular velocity for rigid bodies
		for (int i=0; i<(int)m_fem.m_RBW.size(); ++i)
		{
			FERigidBodyAngularVelocity* pw = m_fem.m_RBW[i];
			if (pw->IsActive())
			{
				FERigidBody& rb = static_cast<FERigidBody&>(*m_fem.Object(i));
				rb.m_wp = rb.m_wt = pw->w;
			}
		}

		// now set the initial velocities of all rigid nodes
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
void FESolidAnalysis::Finish()
{
	FEAnalysis::Finish();

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
bool FEExplicitSolidAnalysis::Init()
{
	// initialize base class data
	FEAnalysis::Init();

	// clear the active rigid body BC's
	int NRB = m_fem.Objects();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (int j=0; j<6; ++j)
		{
			RB.m_pDC[j] = 0;
			RB.m_BC[j] = 0;
		}
	}

	// set the fixed rigid body BC's
	for (int i=0;i<(int) m_fem.m_RBC.size(); ++i)
	{
		FERigidBodyFixedBC* pbc = m_fem.m_RBC[i];
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(pbc->id));
		if (RB.IsActive() && pbc->IsActive()) RB.m_BC[pbc->bc] = -1;
	}

	// set the active rigid bodies BC's
	for (int i=0; i<(int) m_fem.m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *(m_fem.m_RDC[i]);
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(DC.id));
		if (RB.IsActive() && DC.IsActive())
		{
			assert(RB.m_BC[DC.bc] == 0);	// make sure we are not overriding a fixed bc
			RB.m_pDC[DC.bc] = &DC;
			RB.m_BC[DC.bc] = 1;
			DC.ref = 0.0;
			if (DC.brel)
			{
				switch (DC.bc)
				{
				case 0: DC.ref = RB.m_rt.x - RB.m_r0.x; break;
				case 1: DC.ref = RB.m_rt.y - RB.m_r0.y; break;
				case 2: DC.ref = RB.m_rt.z - RB.m_r0.z; break;
				}
			}
		}
	}

	// reset nodal ID's
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = -1;

		if (node.m_BC[DOF_X] != -1) node.m_ID[DOF_X] = 0;
		if (node.m_BC[DOF_Y] != -1) node.m_ID[DOF_Y] = 0;
		if (node.m_BC[DOF_Z] != -1) node.m_ID[DOF_Z] = 0;
		if (node.m_BC[DOF_U] != -1) node.m_ID[DOF_U] = 0;
		if (node.m_BC[DOF_V] != -1) node.m_ID[DOF_V] = 0;
		if (node.m_BC[DOF_W] != -1) node.m_ID[DOF_W] = 0;
	}

	// set the rigid nodes
	// Note that also the rotational degrees of freedom are fixed
	// for rigid nodes that do not belong to a non-rigid shell element.
	int nrn = m_fem.RigidNodes();
	for (int i=0; i<nrn; ++i)
	{
		FERigidNode& rn = *m_fem.RigidNode(i);
		if (rn.IsActive())
		{
			FENode& node = m_fem.GetMesh().Node(rn.nid);
			node.m_rid = rn.rid;

			// fix degrees of freedom
			node.m_ID[DOF_X] = -1;
			node.m_ID[DOF_Y] = -1;
			node.m_ID[DOF_Z] = -1;
			if (node.m_bshell == false)
			{
				node.m_ID[DOF_U] = -1;
				node.m_ID[DOF_V] = -1;
				node.m_ID[DOF_W] = -1;
			}
		}
		else 
		{
			FENode& node = m_fem.GetMesh().Node(rn.nid);
			node.m_rid = -1;
		}
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

	// Sometimes an (ignorant) user might have added a rigid body
	// that is not being used. Since this can cause problems we need
	// to find these rigid bodies.
	int nrb = m_fem.Objects();
	vector<int> mec; mec.assign(nrb, 0);
	for (int i=0; i<m_fem.GetMesh().Nodes(); ++i)
	{
		FENode& node = m_fem.GetMesh().Node(i);
		int n = node.m_rid;
		if (n >= 0) mec[n]++;
	}

	for (int i=0; i<nrb; ++i)
		if (mec[i] == 0)
		{
			FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
			felog.printbox("WARNING", "Rigid body %d is not being used.", RB.m_mat+1);
			RB.Activate(false);
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

	// modify the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il)
		{
			list<FELinearConstraint::SlaveDOF>::iterator is = il->slave.begin();
			for (int i=0; i<(int) il->slave.size(); ++i, ++is)
			{
				is->neq = m_fem.GetMesh().Node(is->node).m_ID[is->bc];
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
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FELinearSolidAnalysis::Init()
{
	// initialize base class data
	FEAnalysis::Init();

	// reset nodal ID's
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = -1;

		if (node.m_BC[DOF_X] != -1) node.m_ID[DOF_X] = 0;
		if (node.m_BC[DOF_Y] != -1) node.m_ID[DOF_Y] = 0;
		if (node.m_BC[DOF_Z] != -1) node.m_ID[DOF_Z] = 0;
		if (node.m_BC[DOF_U] != -1) node.m_ID[DOF_U] = 0;
		if (node.m_BC[DOF_V] != -1) node.m_ID[DOF_V] = 0;
		if (node.m_BC[DOF_W] != -1) node.m_ID[DOF_W] = 0;
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

	// modify the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il)
		{
			list<FELinearConstraint::SlaveDOF>::iterator is = il->slave.begin();
			for (int i=0; i<(int) il->slave.size(); ++i, ++is)
			{
				is->neq = m_fem.GetMesh().Node(is->node).m_ID[is->bc];
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
