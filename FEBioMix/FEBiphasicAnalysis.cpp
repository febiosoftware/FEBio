#include "stdafx.h"
#include "FEBiphasicAnalysis.h"
#include "FEBioMech/FERigidMaterial.h"
#include "FEBioMech/FEContactInterface.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include "FEBiphasic.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
void FEBiphasicAnalysis::InitNodes()
{
	// open all dofs we need
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int i=0; i<(int)node.m_ID.size(); ++i) node.m_ID[i] = -1; 

		// open the dofs for non-fixed nodes
		if (node.m_bexclude == false)
		{
			// Open displacement DOFs for non-rigid nodes
			// (rigid nodes are assigned the rigid body's equation numbers)
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = 0;
				node.m_ID[DOF_Y] = 0;
				node.m_ID[DOF_Z] = 0;
			}

			// Open rotational degrees of freedom for non-rigid shell nodes
			if (node.m_bshell)
			{
				node.m_ID[DOF_U] = 0;
				node.m_ID[DOF_V] = 0;
				node.m_ID[DOF_W] = 0;
			}

			// Open pressure DOF for all nodes
			// (These will be closed below for non-poro nodes)
			node.m_ID[DOF_P] = 0;
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

	// fix all mixture dofs that are not used that is, that are not part of a biphasic material.
	// This is done in two steps.
	// step 1. mark all pressure nodes
	vector<int> tag;
	tag.assign(mesh.Nodes(), 0);
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		FEBiphasic* bm  = dynamic_cast<FEBiphasic*>(dom.GetMaterial());
		if (bm)
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (int j=0; j<N; ++j) 
					if (mesh.Node(n[j]).m_ID[DOF_P] != DOF_FIXED) tag[ n[j] ] = 1;
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
				if (tag[ n[j] ] != 1) mesh.Node(n[j]).m_ID[DOF_P] = DOF_FIXED;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FEBiphasicAnalysis::Activate()
{
	// initialize base class data
	FEAnalysis::Activate();

	// reset nodal ID's
	InitNodes();

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
			switch (bc)
			{
			case DOF_X: DC.r = br ? node.m_rt.x - node.m_r0.x : 0; break;
			case DOF_Y: DC.r = br ? node.m_rt.y - node.m_r0.y : 0; break;
			case DOF_Z: DC.r = br ? node.m_rt.z - node.m_r0.z : 0; break;
			case DOF_U: DC.r = br ? node.m_Dt.x - node.m_D0.x : 0; break;
			case DOF_V: DC.r = br ? node.m_Dt.y - node.m_D0.y : 0; break;
			case DOF_W: DC.r = br ? node.m_Dt.z - node.m_D0.z : 0; break;
			case DOF_P: DC.r = br ? node.m_pt   - node.m_p0   : 0; break;
			}
		}
	}

	// modify the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il) il->Activate();
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
