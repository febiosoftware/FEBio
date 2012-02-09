// FESolidSolver.cpp: implementation of the FESolidSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESolver.h"
#include "fem.h"
#include "FEBioLib/log.h"
#include "NumCore/SkylineSolver.h"
#include "NumCore/LUSolver.h"
#include "NumCore/ConjGradIterSolver.h"
#include "FEBioLib/PSLDLTSolver.h"
#include "FEBioLib/SuperLUSolver.h"
#include "FEBioLib/SuperLU_MT_Solver.h"
#include "FEBioLib/PardisoSolver.h"
#include "FEBioLib/WSMPSolver.h"
#include "FEBioLib/RCICGSolver.h"

//-----------------------------------------------------------------------------
//! Constructor of FESolver base class
FESolver::FESolver(FEModel& fem) : FENLSolver(fem)
{
	// Stiffness matrix and linear solver are allocated in Init()
	m_pK = 0;
	m_plinsolve = 0;
	m_neq = 0;
	m_niter = 0;
}

//-----------------------------------------------------------------------------
//! Destructor.
FESolver::~FESolver()
{
	delete m_pK;		// clean up stiffnes matrix data
	delete m_plinsolve;	// clean up linear solver data
}

//-----------------------------------------------------------------------------
//! Initialization
// TODO: where is this function called?
bool FESolver::Init()
{
	FEM& fem = dynamic_cast<FEM&>(m_fem);

	// Now that we have determined the equation numbers we can continue
	// with creating the stiffness matrix. First we select the linear solver
	// The stiffness matrix is created in CreateStiffness
	// Note that if a particular solver was requested in the input file
	// then the solver might already be allocated. That's way we need to check it.
	if (m_plinsolve == 0)
	{
		switch (fem.m_nsolver)
		{
		case SKYLINE_SOLVER      : m_plinsolve = new SkylineSolver(); break;
		case PSLDLT_SOLVER       : m_plinsolve = new PSLDLTSolver (); break;
		case SUPERLU_SOLVER      : m_plinsolve = new SuperLUSolver(); break;
		case SUPERLU_MT_SOLVER   : m_plinsolve = new SuperLU_MT_Solver(); break;
		case PARDISO_SOLVER      : m_plinsolve = new PardisoSolver(); break;
		case LU_SOLVER           : m_plinsolve = new LUSolver(); break;
		case WSMP_SOLVER         : m_plinsolve = new WSMPSolver(); break;
		case CG_ITERATIVE_SOLVER : m_plinsolve = new ConjGradIterSolver(); break;
		case RCICG_SOLVER        : m_plinsolve = new RCICGSolver(); break;
		default:
			clog.printbox("FATAL ERROR","Unknown solver type selected\n");
			return false;
		}
	}

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the solver allocate the correct type of matrix format
	SparseMatrix* pS = m_plinsolve->CreateSparseMatrix(m_bsymm? SPARSE_SYMMETRIC : SPARSE_UNSYMMETRIC);
	if (pS == 0)
	{
		clog.printbox("FATAL ERROR", "The selected linear solver does not support the requested\n matrix format.\nPlease select a different linear solver.\n");
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
		clog.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Clean
// TODO: why can this not be done in destructor?
void FESolver::Clean()
{
	if (m_plinsolve) m_plinsolve->Destroy();
}

//-----------------------------------------------------------------------------
//!  Creates the global stiffness matrix
// TODO: Can we move this to the FEStiffnessMatrix::Create function?
bool FESolver::CreateStiffness(bool breset)
{
	// clean up the solver
	if (m_pK->NonZeroes()) m_plinsolve->Destroy();

	// clean up the stiffness matrix
	m_pK->Clear();

	// create the stiffness matrix
	clog.printf("===== reforming stiffness matrix:\n");
	if (m_pK->Create(this, m_neq, breset) == false) 
	{
		clog.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
		return false;
	}
	else
	{
		// output some information about the direct linear solver
		int neq = m_pK->Rows();
		int nnz = m_pK->NonZeroes();
		clog.printf("\tNr of equations ........................... : %d\n", neq);
		clog.printf("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
		clog.printf("\n");
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

//-----------------------------------------------------------------------------
//! TODO: This function is only used for rigid joints. I need to figure out if
//!       I can use the other assembly function.
void FESolver::AssembleStiffness(std::vector<int>& lm, matrix& ke)
{
	m_pK->Assemble(ke, lm);
}

//-----------------------------------------------------------------------------
//! TODO: This function is only used by the rigid joints. I need to figure out
//!       if I can use the ohter residual function
void FESolver::AssembleResidual(vector<int>& lm, vector<double>& fe, vector<double>& R)
{
	int n = lm.size();
	assert(fe.size() == n);
	for (int i=0; i<n; ++i) if (lm[i] >= 0) R[lm[i]] += fe[i];
}
