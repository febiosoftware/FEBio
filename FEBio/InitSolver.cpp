#include "stdafx.h"
#include "FESolver.h"
#include "fem.h"

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolver::Init
// Allocates and initializes the data structures used by the FESolver
//

bool FESolver::Init()
{
	// Now that we have determined the equation numbers we can continue
	// with creating the stiffness matrix. First we select the linear solver
	// The stiffness matrix is created in CreateStiffness
	// Note that if a particular solver was requested in the input file
	// then the solver might already be allocated. That's way we need to check it.
	if (m_psolver == 0)
	{
		switch (m_fem.m_nsolver)
		{
		case SKYLINE_SOLVER: m_psolver = new SkylineSolver(); break;
		case PSLDLT_SOLVER : m_psolver = new PSLDLTSolver (); break;
		case SUPERLU_SOLVER: m_psolver = new SuperLUSolver(); break;
		case SUPERLU_MT_SOLVER: m_psolver = new SuperLU_MT_Solver(); break;
		case PARDISO_SOLVER: m_psolver = new PardisoSolver(); break;
		case LU_SOLVER     : m_psolver = new LUSolver(); break;
		case WSMP_SOLVER   : m_psolver = new WSMPSolver(); break;
		case CG_ITERATIVE_SOLVER : m_psolver = new ConjGradIterSolver(); break;
		default:
			m_log.printbox("FATAL ERROR","Unknown solver type selected\n");
			return false;
		}
	}

	// create the stiffness matrix
	// we let the solver allocate the correct type of matrix format
	// Note that this does not construct the stiffness matrix. This
	// is done later in the StiffnessMatrix routine
	m_pK = new FEStiffnessMatrix(m_psolver->GetMatrix());
	if (m_pK == 0)
	{
		m_log.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
		return false;
	}


	// get nr of equations
	int neq = m_fem.m_neq;

	// allocate vectors
	m_Fn.create(neq); m_Fn.zero();
	m_Fd.create(neq); m_Fd.zero();
	m_Fr.create(neq); m_Fr.zero();
	m_ui.create(neq); m_ui.zero();
	m_Ui.create(neq); m_Ui.zero();
	m_Ut.create(neq); m_Ut.zero();
	m_R0.create(neq); m_R0.zero();
	m_R1.create(neq); m_R1.zero();

	// we need to fill the total displacement vector m_Ut
	// TODO: I need to find an easier way to do this
	FEMesh& mesh = m_fem.m_mesh;
	int i, n;
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// displacement dofs
		n = node.m_ID[0]; if (n >= 0) m_Ut[n] = node.m_rt.x - node.m_r0.x;
		n = node.m_ID[1]; if (n >= 0) m_Ut[n] = node.m_rt.y - node.m_r0.y;
		n = node.m_ID[2]; if (n >= 0) m_Ut[n] = node.m_rt.z - node.m_r0.z;

		// rotational dofs
		n = node.m_ID[3]; if (n >= 0) m_Ut[n] = node.m_Dt.x - node.m_D0.x;
		n = node.m_ID[4]; if (n >= 0) m_Ut[n] = node.m_Dt.y - node.m_D0.y;
		n = node.m_ID[5]; if (n >= 0) m_Ut[n] = node.m_Dt.z - node.m_D0.z;

		// pressure dofs
		n = node.m_ID[6]; if (n >= 0) m_Ut[n] = node.m_pt;
	}

	// allocate storage for BFGS update vectors
	m_V.Create(m_maxups, neq);
	m_W.Create(m_maxups, neq);

	m_D .create(neq);
	m_G .create(neq);
	m_H .create(neq);

	// set the create stiffness matrix flag
	m_breshape = true;

	// output some information about the direct linear solver
	m_log.printf(" LINEAR EQUATION SOLVER DATA\n");
	m_log.printf("===========================================================================\n");
	m_log.printf("\tSolver type ............................... : ");
	if (m_fem.m_nsolver == SKYLINE_SOLVER     ) m_log.printf("Skyline\n");
	if (m_fem.m_nsolver == PSLDLT_SOLVER      ) m_log.printf("PSLDLT\n");
	if (m_fem.m_nsolver == SUPERLU_SOLVER     ) m_log.printf("SuperLU\n");
	if (m_fem.m_nsolver == SUPERLU_MT_SOLVER  ) m_log.printf("SuperLU_MT\n");
	if (m_fem.m_nsolver == PARDISO_SOLVER  ) m_log.printf("Pardiso\n");
	if (m_fem.m_nsolver == LU_SOLVER          ) m_log.printf("LUSolver\n");
	if (m_fem.m_nsolver == CG_ITERATIVE_SOLVER) m_log.printf("Conjugate gradient\n");
	m_log.printf("\tNr of equations ........................... : %d\n", neq);
//	m_log.printf("\tNr of nonzero coefficients ................ : %d\n", m_pK->NonZeroes());
//	m_log.printf("\tin global stiffness matrix\n");
//	m_log.printf("\n\n");

	return true;
}
