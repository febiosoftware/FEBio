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
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = -1;

		// open the dofs for non-fixed nodes
		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = 0;
				node.m_ID[DOF_Y] = 0;
				node.m_ID[DOF_Z] = 0;
				node.m_ID[DOF_T] = 0;
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
			int dof = DC.GetDOF();
			for (size_t j=0; j<DC.Items(); ++j)
			{
				FENode& node = mesh.Node(DC.NodeID(j));
				node.m_ID[dof] = DOF_PRESCRIBED;
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

	// Now we adjust the equation numbers of prescribed dofs according to the above rule
	// Make sure that a prescribed dof has not been fixed
	// TODO: maybe this can be moved to the FESolver::InitEquations function
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		if (DC.IsActive()) DC.Update();
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
