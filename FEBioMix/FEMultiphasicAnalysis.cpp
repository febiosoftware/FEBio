#include "stdafx.h"
#include "FEMultiphasicAnalysis.h"
#include "FEBioMech/FERigidMaterial.h"
#include "FEBioMech/FEContactInterface.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
void FEMultiphasicAnalysis::InitNodes()
{
    // get nodal DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
	// open all dofs we need
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int i=0; i<(int)node.m_ID.size(); ++i) node.m_ID[i] = -1;

		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = 0;
				node.m_ID[DOF_Y] = 0;
				node.m_ID[DOF_Z] = 0;
				node.m_ID[DOF_P] = 0;

				for (int k=0; k<MAX_CDOFS; ++k) {
					int dofc = DOF_C + k;
					node.m_ID[dofc] = 0;
				}
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
		mesh.Node(bc.m_node).m_ID[bc.m_dof] = -1;
	}

	// fix all mixture dofs that are not used that is, that are not part of a biphasic material.
	// This is done in three steps.
	// step 1. mark all pressure nodes
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		FEBiphasic*		  bm  = dynamic_cast<FEBiphasic*      >(dom.GetMaterial());
		FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(dom.GetMaterial());
		FETriphasic*      btm = dynamic_cast<FETriphasic*     >(dom.GetMaterial());
		FEMultiphasic*    bmm = dynamic_cast<FEMultiphasic*   >(dom.GetMaterial());
		if (bm || bsm || btm || bmm)
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (int j=0; j<N; ++j) 
					if (mesh.Node(n[j]).m_ID[DOF_P] == 0) mesh.Node(n[j]).m_ID[DOF_P] = 1;
			}
		}
		if (bsm)
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (int j=0; j<N; ++j) {
					int dofc = DOF_C + bsm->GetSolute()->GetSoluteID();
					if (mesh.Node(n[j]).m_ID[dofc] == 0) mesh.Node(n[j]).m_ID[dofc] = 1;
				}
			}
		}
		else if (btm)
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (int j=0; j<N; ++j) {
					for (int k=0; k<2; ++k) {
						int dofc = DOF_C + btm->m_pSolute[k]->GetSoluteID();
						if (mesh.Node(n[j]).m_ID[dofc] == 0) mesh.Node(n[j]).m_ID[dofc] = 1;
					}
				}
			}
		}
		else if (bmm)
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				int nsol = bmm->Solutes();
				for (int j=0; j<N; ++j) {
					for (int k=0; k<nsol; ++k) {
						int dofc = DOF_C + bmm->GetSolute(k)->GetSoluteID();
						if (mesh.Node(n[j]).m_ID[dofc] == 0) mesh.Node(n[j]).m_ID[dofc] = 1;
					}
				}
			}
		}	
	}
	
	// step 2. fix pressure dofs of all unmarked nodes
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		for (int i=0; i<dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			int N = el.Nodes();
			int* n = &el.m_node[0];
			for (int j=0; j<N; ++j) {
				if (mesh.Node(n[j]).m_ID[DOF_P] != 1) mesh.Node(n[j]).m_ID[DOF_P] = -1;
				for (int k=0; k<MAX_CDOFS; ++k) {
					int dofc = DOF_C + k;
					if (mesh.Node(n[j]).m_ID[dofc] != 1) mesh.Node(n[j]).m_ID[dofc] = -1;
				}
			}
		}
	}
	
	// step 3. free all marked dofs
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		for (int i=0; i<dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			int N = el.Nodes();
			int* n = &el.m_node[0];
			for (int j=0; j<N; ++j) {
				if (mesh.Node(n[j]).m_ID[DOF_P] == 1) mesh.Node(n[j]).m_ID[DOF_P] = 0;
				for (int k=0; k<MAX_CDOFS; ++k) {
					int dofc = DOF_C + k;
					if (mesh.Node(n[j]).m_ID[dofc] == 1) mesh.Node(n[j]).m_ID[dofc] = 0;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FEMultiphasicAnalysis::Activate()
{
	// initialize base class data
	FEAnalysis::Activate();

	// clear the active rigid body BC's
	int NRB = m_fem.Objects();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (int j=0; j<6; ++j)
		{
			RB.m_pDC[j] = 0;
			RB.m_LM[j] = 0;
		}
	}

	// set the fixed rigid body BC's
	for (int i=0;i<(int) m_fem.m_RBC.size(); ++i)
	{
		FERigidBodyFixedBC* pbc = m_fem.m_RBC[i];
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(pbc->id));
		if (pbc->IsActive()) RB.m_LM[pbc->bc] = -1;
	}

	// set the active rigid bodies BC's
	for (int i=0; i<(int) m_fem.m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *(m_fem.m_RDC[i]);
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(DC.id));
		if (DC.IsActive())
		{
			RB.m_pDC[DC.bc] = &DC;
			RB.m_LM[DC.bc] = 0;		// make sure this dof is free
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
	InitNodes();

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

	// modify the (aug lag) nonlinear constraints
	// TODO: I think this is already done in FEM::Init. Why do I need to do this again?
	int M = m_fem.NonlinearConstraints();
	for (int m=0; m<M; ++m) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(m);
		if (plc->IsActive()) plc->Init();
	}

	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		felog.printbox("FATAL ERROR","Failed to initialize solver.\nAborting run.\n");
		return false;
	}

	return true;
}
