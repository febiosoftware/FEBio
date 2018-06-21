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
FGMRESSolver::FGMRESSolver() : m_pA(0)
{
	m_maxiter = 0; // use default min(N, 150)
	m_print_level = 0;
	m_doResidualTest = true;
	m_doZeroNormTest = true;
	m_tol = 0.0;
	m_nrestart = 0; // use default = maxiter
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
void FGMRESSolver::SetSparseMatrix(SparseMatrix* pA)
{
	m_pA = pA;
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

	return true; 
#else
	return false;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::BackSolve(vector<double>& x, vector<double>& b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// number of equations
	MKL_INT N = m_pA->Rows();

	// data allocation
	MKL_INT ipar[128];
	double dpar[128];
	MKL_INT RCI_request;
	int M = (N < 150 ? N : 150); // this is the default value of ipar[4] and ipar[14]

	int nrestart = M;
	if (m_nrestart > 0) nrestart = m_nrestart;
	else if (m_maxiter > 0) nrestart = m_maxiter;

	int maxIter = M;
	if (m_maxiter > 0) maxIter = m_maxiter;

	double* ptmp = &m_tmp[0];
	MKL_INT ivar = N;

	// initialize the solver
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// Set the desired parameters:
	ipar[ 4] = maxIter;	                        // max number of iterations
	ipar[ 8] = (m_doResidualTest ? 1 : 0);		// do residual stopping test
	ipar[14] = nrestart;	                    // number of non-restarted iterations
	ipar[ 9] = 0;								// do not request for the user defined stopping test
	ipar[11] = (m_doZeroNormTest ? 1 : 0);		// do the check of the norm of the next generated vector automatically
	if (m_tol > 0) dpar[0] = m_tol;				// set the relative tolerance

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request!=0) { MKL_Free_Buffers(); return false; }

	// solve the problem
	bool bdone = false;
	bool bconverged = false;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);

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
		case 4:
			break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// get the solution. 
	MKL_INT itercount;
	dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp, &itercount);
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
	return BackSolve(x, b);
}
