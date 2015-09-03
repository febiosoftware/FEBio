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
#include "FECore/FEModel.h"

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
		for (int i=0; i<(int)node.m_ID.size(); ++i) node.m_ID[i] = DOF_FIXED;

		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = node.m_BC[DOF_X];
				node.m_ID[DOF_Y] = node.m_BC[DOF_Y];
				node.m_ID[DOF_Z] = node.m_BC[DOF_Z];
			}

			if (node.m_bshell)
			{
				node.m_ID[DOF_U] = node.m_BC[DOF_U];
				node.m_ID[DOF_V] = node.m_BC[DOF_V];
				node.m_ID[DOF_W] = node.m_BC[DOF_W];
			}

			node.m_ID[DOF_P] = node.m_BC[DOF_P];

			for (int k=0; k<MAX_CDOFS; ++k) {
				int dofc = DOF_C + k;
				node.m_ID[dofc] = node.m_BC[dofc];
			}
		}
	}

	// fix all mixture dofs that are not used that is, that are not part of a biphasic material.
	const int NN = mesh.Nodes();
	vector<int> tag;

	// do the pressure dofs first
	tag.assign(NN, 0);
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
					if (mesh.Node(n[j]).m_ID[DOF_P] != DOF_FIXED) tag[ n[j] ] = 1;
			}
		}
	}

	// fix pressure dofs of all unmarked nodes
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

	// now do the concentration dofs
	for (int k=0; k<MAX_CDOFS; ++k)
	{
		int dofc = DOF_C + k;
		tag.assign(NN, 0);
		for (int nd = 0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& dom = mesh.Domain(nd);

			// get the material
			FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(dom.GetMaterial());
			FETriphasic*      btm = dynamic_cast<FETriphasic*     >(dom.GetMaterial());
			FEMultiphasic*    bmm = dynamic_cast<FEMultiphasic*   >(dom.GetMaterial());

			// see if the material has this dof
			bool has_c = false;
			if (bsm && (dofc == DOF_C + bsm->GetSolute()->GetSoluteID())) has_c = true;
			if (btm)
			{
				if (dofc == DOF_C + btm->m_pSolute[0]->GetSoluteID()) has_c = true;
				if (dofc == DOF_C + btm->m_pSolute[1]->GetSoluteID()) has_c = true;
			}
			if (bmm)
			{
				int nsol = bmm->Solutes();
				for (int l=0; l<nsol; ++l)
					if (dofc == DOF_C + bmm->GetSolute(l)->GetSoluteID()) has_c = true;
			}

			// if so, mark all non-fixed dofs
			if (has_c)
			{
				for (int i=0; i<dom.Elements(); ++i)
				{
					FEElement& el = dom.ElementRef(i);
					int N = el.Nodes();
					int* n = &el.m_node[0];
					for (int j=0; j<N; ++j) {
						if (mesh.Node(n[j]).m_ID[dofc] != DOF_FIXED) tag[ n[j] ] = 1;
					}
				}
			}
		}
	
		// step 2. fix concentration dofs of all unmarked nodes
		for (int nd = 0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& dom = mesh.Domain(nd);
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (int j=0; j<N; ++j) {
					if (tag[ n[j] ] != 1) mesh.Node(n[j]).m_ID[dofc] = DOF_FIXED;
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

	// reset nodal ID's
	InitNodes();

	// initialize equations
	// ----->
	if (m_psolver->InitEquations() == false) return false;

	// initialize linear constraints
	// Must be done after equations are initialized
	if (InitLinearConstraints() == false) return false;
	// ----->

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
