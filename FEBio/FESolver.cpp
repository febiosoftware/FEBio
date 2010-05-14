// FESolidSolver.cpp: implementation of the FESolidSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESolver.h"
#include "fem.h"
#include "log.h"
#include "FECore/SkylineSolver.h"
#include "FECore/PSLDLTSolver.h"
#include "FECore/SuperLUSolver.h"
#include "FECore/SuperLU_MT_Solver.h"
#include "FECore/LUSolver.h"
#include "FECore/PardisoSolver.h"
#include "FECore/WSMPSolver.h"
#include "FECore/ConjGradIterSolver.h"

FESolver::FESolver(FEM& fem) : m_fem(fem)
{
	// Stiffness matrix and linear solver are allocated in Init()
	m_pK = 0;
	m_plinsolve = 0;

	m_Rmin = 1.0e-20;
}


FESolver::~FESolver()
{
	delete m_pK;		// clean up stiffnes matrix data
	delete m_plinsolve;	// clean up linear solver data
}

bool FESolver::Init()
{
	// get the logfile
	Logfile& log = GetLogfile();

	// Now that we have determined the equation numbers we can continue
	// with creating the stiffness matrix. First we select the linear solver
	// The stiffness matrix is created in CreateStiffness
	// Note that if a particular solver was requested in the input file
	// then the solver might already be allocated. That's way we need to check it.
	if (m_plinsolve == 0)
	{
		switch (m_fem.m_nsolver)
		{
		case SKYLINE_SOLVER      : m_plinsolve = new SkylineSolver(); break;
		case PSLDLT_SOLVER       : m_plinsolve = new PSLDLTSolver (); break;
		case SUPERLU_SOLVER      : m_plinsolve = new SuperLUSolver(); break;
		case SUPERLU_MT_SOLVER   : m_plinsolve = new SuperLU_MT_Solver(); break;
		case PARDISO_SOLVER      : m_plinsolve = new PardisoSolver(); break;
		case LU_SOLVER           : m_plinsolve = new LUSolver(); break;
		case WSMP_SOLVER         : m_plinsolve = new WSMPSolver(); break;
		case CG_ITERATIVE_SOLVER : m_plinsolve = new ConjGradIterSolver(); break;
		default:
			log.printbox("FATAL ERROR","Unknown solver type selected\n");
			return false;
		}
	}

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the solver allocate the correct type of matrix format
	SparseMatrix* pS = m_plinsolve->CreateSparseMatrix(m_fem.m_bsymm? SPARSE_SYMMETRIC : SPARSE_UNSYMMETRIC);
	if (pS == 0)
	{
		log.printbox("FATAL ERROR", "The selected linear solver does not support the requested\n matrix format.\nPlease select a different linear solver.\n");
		return false;
	}

	// clean up the stiffness matrix if we have one
	if (m_pK) delete m_pK; m_pK = 0;

	// Create the stiffness matrix.
	// Note that this does not construct the stiffness matrix. This
	// is done later in the StiffnessMatrix routine.
	m_pK = new FEStiffnessMatrix(pS);
	if (m_pK == 0)
	{
		log.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FESolver::Clean()
{
	if (m_plinsolve) m_plinsolve->Destroy();
}

//-----------------------------------------------------------------------------
// FUNCTION: FESolver::CreateStiffness
//  Creates the global stiffness matrix
//

bool FESolver::CreateStiffness(bool breset)
{
	// clean up the solver
	if (m_pK->NonZeroes()) m_plinsolve->Destroy();	// GAA

	// clean up the stiffness matrix
	m_pK->Clear();

	// get the logfile
	Logfile& log = GetLogfile();

	// create the stiffness matrix
	log.printf("===== reforming stiffness matrix:\n");
	if (m_pK->Create(m_fem, breset) == false) 
	{
		log.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
		return false;
	}
	else
	{
		// output some information about the direct linear solver
		int neq = m_fem.m_neq;
		int nnz = m_pK->NonZeroes();
		log.printf("\tNr of equations ........................... : %d\n", neq);
		log.printf("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
		log.printf("\n");
	}

	// Do the preprocessing of the solver
	m_SolverTime.start();
	{
		if (!m_plinsolve->PreProcess()) throw FatalError();
	}
	m_SolverTime.stop();

	// done!
	return true;
}
