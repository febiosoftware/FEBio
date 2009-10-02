// FESolidSolver.cpp: implementation of the FESolidSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESolver.h"
#include "fem.h"
#include "FECore/SkylineSolver.h"
#include "FECore/PSLDLTSolver.h"
#include "FECore/SuperLUSolver.h"
#include "FECore/SuperLU_MT_Solver.h"
#include "FECore/LUSolver.h"
#include "FECore/PardisoSolver.h"
#include "FECore/WSMPSolver.h"
#include "FECore/ConjGradIterSolver.h"

FESolver::FESolver(FEM& fem) : m_fem(fem), m_log(fem.m_log)
{
	// Stiffness matrix and linear solver are allocated in Init()
	m_pK = 0;
	m_psolver = 0;
}


FESolver::~FESolver()
{
	delete m_pK;		// clean up stiffnes matrix data
	delete m_psolver;	// clean up linear solver data
}

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
		case SKYLINE_SOLVER      : m_psolver = new SkylineSolver(); break;
		case PSLDLT_SOLVER       : m_psolver = new PSLDLTSolver (); break;
		case SUPERLU_SOLVER      : m_psolver = new SuperLUSolver(); break;
		case SUPERLU_MT_SOLVER   : m_psolver = new SuperLU_MT_Solver(); break;
		case PARDISO_SOLVER      : m_psolver = new PardisoSolver(); break;
		case LU_SOLVER           : m_psolver = new LUSolver(); break;
		case WSMP_SOLVER         : m_psolver = new WSMPSolver(); break;
		case CG_ITERATIVE_SOLVER : m_psolver = new ConjGradIterSolver(); break;
		default:
			m_log.printbox("FATAL ERROR","Unknown solver type selected\n");
			return false;
		}
	}

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the solver allocate the correct type of matrix format
	SparseMatrix* pS = m_psolver->GetMatrix(m_fem.m_bsymm? SPARSE_SYMMETRIC : SPARSE_UNSYMMETRIC);
	if (pS == 0)
	{
		m_log.printbox("FATAL ERROR", "The selected linear solver does not support the requested\n matrix format.\nPlease select a different linear solver.\n");
		return false;
	}


	// Create the stiffness matrix.
	// Note that this does not construct the stiffness matrix. This
	// is done later in the StiffnessMatrix routine.
	m_pK = new FEStiffnessMatrix(pS);
	if (m_pK == 0)
	{
		m_log.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
		return false;
	}

	return true;
}


///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolver::CreateStiffness
//  Creates the global stiffness matrix
//

bool FESolver::CreateStiffness(bool breset)
{
	// clean up the solver
	if (m_pK->NonZeroes()) m_psolver->Destroy(*m_pK);	// GAA

	// clean up the stiffness matrix
	m_pK->Clear();

	// create the stiffness matrix
	m_log.printf("===== reforming stiffness matrix:\n");
	if (m_pK->Create(m_fem, breset) == false) 
	{
		m_log.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
		return false;
	}
	else
	{
		// output some information about the direct linear solver
		int neq = m_fem.m_neq;
		int nnz = m_pK->NonZeroes();
		m_log.printf("\tNr of equations ........................... : %d\n", neq);
		m_log.printf("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
		m_log.printf("\n");
	}

	// Do the preprocessing of the solver
	m_SolverTime.start();
	{
		if (!m_psolver->PreProcess(*m_pK)) throw FatalError();
	}
	m_SolverTime.stop();

	// done!
	return true;
}
