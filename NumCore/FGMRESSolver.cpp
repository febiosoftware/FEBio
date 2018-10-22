#include "stdafx.h"
#include "FGMRESSolver.h"
#include "CompactSymmMatrix.h"
#include "CompactUnSymmMatrix.h"

//-----------------------------------------------------------------------------
// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef MKL_ISS
#ifdef PARDISO
#undef PARDISO
#endif
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#endif // MKL_ISS

//-----------------------------------------------------------------------------
FGMRESSolver::FGMRESSolver(FEModel* fem) : IterativeLinearSolver(fem), m_pA(0)
{
	m_maxiter = 0; // use default min(N, 150)
	m_print_level = 0;
	m_doResidualTest = true;
	m_doZeroNormTest = true;
	m_tol = 0.0;
	m_nrestart = 0; // use default = maxiter

	m_doPreCond = false;
	m_P = 0; // we don't use a preconditioner for this solver
}

//-----------------------------------------------------------------------------
// set the preconditioner
void FGMRESSolver::SetPreconditioner(Preconditioner* P)
{
	if (m_P) delete m_P;
	m_P = P;
	m_doPreCond = (m_P != 0);
}

//-----------------------------------------------------------------------------
//! Set max nr of iterations
void FGMRESSolver::SetMaxIterations(int n)
{
	m_maxiter = n;
}

//-----------------------------------------------------------------------------
//! Set the nr of non-restarted iterations
void FGMRESSolver::SetNonRestartedIterations(int n)
{
	m_nrestart = n;
}

//-----------------------------------------------------------------------------
// Set the print level
void FGMRESSolver::SetPrintLevel(int n)
{
	m_print_level = n;
}

//-----------------------------------------------------------------------------
// set residual stopping test flag
void FGMRESSolver::DoResidualStoppingTest(bool b)
{
	m_doResidualTest = b;
}

//-----------------------------------------------------------------------------
// set zero norm stopping test flag
void FGMRESSolver::DoZeroNormStoppingTest(bool b)
{
	m_doZeroNormTest = b;
}

//-----------------------------------------------------------------------------
// set the convergence tolerance for the residual stopping test
void FGMRESSolver::SetResidualTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
//! This solver does not use a preconditioner
bool FGMRESSolver::HasPreconditioner() const 
{
	return m_P != 0; 
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRESSolver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	// Cleanup if necessary
	if (m_pA) delete m_pA; 
	m_pA = 0;

	// allocate new matrix
	switch(ntype)
	{
	case REAL_SYMMETRIC  : m_pA = new CompactSymmMatrix(0); break;
	case REAL_UNSYMMETRIC: m_pA = new CRSSparseMatrix(1); break;
	}

	// return the matrix (Can be null if matrix format not supported!)
	return m_pA;
#else
	return 0;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::SetSparseMatrix(SparseMatrix* pA)
{
	m_pA = pA;
	return (m_pA != 0);
}

//-----------------------------------------------------------------------------
//! Clean up
void FGMRESSolver::Destroy()
{
	m_tmp.clear();
	m_tmp.shrink_to_fit();
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::PreProcess() 
{
#ifdef MKL_ISS
	// number of equations
	MKL_INT N = m_pA->Rows();

	int M = (N < 150 ? N : 150); // this is the default value of ipar[14]

	if (m_nrestart > 0) M = m_nrestart;
	else if (m_maxiter > 0) M = m_maxiter;

	// allocate temp storage
	m_tmp.resize((N*(2 * M + 1) + (M*(M + 9)) / 2 + 1));

	// set the pre-conditioning flag
	m_doPreCond = (m_P != 0);
	return true; 
#else
	return false;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::BackSolve(double* x, double* b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// number of equations
	MKL_INT N = m_pA->Rows();

	// data allocation
	int M = (N < 150 ? N : 150); // this is the default value of ipar[4] and ipar[14]

	int nrestart = M;
	if (m_nrestart > 0) nrestart = m_nrestart;
	else if (m_maxiter > 0) nrestart = m_maxiter;

	int maxIter = M;
	if (m_maxiter > 0) maxIter = m_maxiter;

	// initialize the solver
	MKL_INT ipar[128] = { 0 };
	double dpar[128] = { 0.0 };
	MKL_INT ivar = N;
	MKL_INT RCI_request;
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, &m_tmp[0]);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	if (m_doPreCond && m_P)
	{
		if (m_P->Create(m_pA) == false)
		{
			MKL_Free_Buffers();
			return false;
		}
		m_doPreCond = false;
	}

	// Set the desired parameters:
	ipar[ 4] = maxIter;	                        // max number of iterations
	ipar[ 7] = 1;								// do the stopping test for maximal number of iterations
	ipar[ 8] = (m_doResidualTest ? 1 : 0);		// do residual stopping test
	ipar[ 9] = 0;								// do not request for the user defined stopping test
	ipar[10] = (m_P != 0 ? 1 : 0);				// do the pre-conditioned version of the FGMRES iterative solver
	ipar[11] = (m_doZeroNormTest ? 1 : 0);		// do the check of the norm of the next generated vector automatically
	ipar[14] = nrestart;	                    // number of non-restarted iterations
	if (m_tol > 0) dpar[0] = m_tol;				// set the relative tolerance

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, &m_tmp[0]);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// solve the problem
	bool bdone = false;
	bool bconverged = false;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, &m_tmp[0]);

		switch (RCI_request)
		{
		case 0: // solution converged. 
			bdone = true;
			bconverged = true;
			break;
		case 1:
			{
				// do matrix-vector multiplication
				m_pA->mult_vector(&m_tmp[ipar[21] - 1], &m_tmp[ipar[22] - 1]);

				if (m_print_level == 1)
				{
					fprintf(stderr, "%3d = %lg (%lg), %lg (%lg)\n", ipar[3], dpar[4], dpar[3], dpar[6], dpar[7]);
				}
			}
			break;
		case 3:	// do the pre-conditioning step
			{
				assert(m_P);
				m_P->mult_vector(&m_tmp[ipar[21] - 1], &m_tmp[ipar[22] - 1]);
			}
			break;
		case 4:
			break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// get the solution. 
	MKL_INT itercount;
	dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, &m_tmp[0], &itercount);
	if (m_print_level > 0)
	{
		fprintf(stderr, "%3d = %lg (%lg), %lg (%lg)\n", ipar[3], dpar[4], dpar[3], dpar[6], dpar[7]);
	}

	MKL_Free_Buffers();
	return bconverged;

#else
	return false;
#endif // MKL_ISS
}

//! convenience function for solving linear system Ax = b
bool FGMRESSolver::Solve(SparseMatrix* A, vector<double>& x, vector<double>& b)
{
	SetSparseMatrix(A);
	if (PreProcess() == false) return false;
	if (Factor() == false) return false;
	return BackSolve(&x[0], &b[0]);
}
