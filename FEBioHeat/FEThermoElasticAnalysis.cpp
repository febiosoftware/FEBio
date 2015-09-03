#include "FEThermoElasticAnalysis.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
FEThermoElasticAnalysis::FEThermoElasticAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_THERMO_ELASTIC)
{
}

//-----------------------------------------------------------------------------
bool FEThermoElasticAnalysis::Activate()
{
	// initialize base class data
	FEAnalysis::Activate();

	// reset nodal ID's
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// fix all degrees of freedom
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = DOF_FIXED;

		// open the dofs for non-fixed nodes
		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = node.m_BC[DOF_X];
				node.m_ID[DOF_Y] = node.m_BC[DOF_Y];
				node.m_ID[DOF_Z] = node.m_BC[DOF_Z];
				node.m_ID[DOF_T] = node.m_BC[DOF_T];
			}

			if (node.m_bshell)
			{
				node.m_ID[DOF_U] = node.m_BC[DOF_U];
				node.m_ID[DOF_V] = node.m_BC[DOF_V];
				node.m_ID[DOF_W] = node.m_BC[DOF_W];
			}
		}
	}

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
