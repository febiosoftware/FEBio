#include "stdafx.h"
#include "RCICGSolver.h"

//-----------------------------------------------------------------------------
// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef MKL_ISS
#ifdef PARDISO
#undef PARDISO
#endif
#include "mkl_solver.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#endif // MKL_ISS

//-----------------------------------------------------------------------------
RCICGSolver::RCICGSolver() : m_pA(0)
{
}

//-----------------------------------------------------------------------------
SparseMatrix* RCICGSolver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	if (ntype != SPARSE_SYMMETRIC) return 0;
	m_pA = new CompactSymmMatrix(1);
	return m_pA;
#else
	return 0;
#endif
}

//-----------------------------------------------------------------------------
bool RCICGSolver::PreProcess()
{
	return true;
}

//-----------------------------------------------------------------------------
bool RCICGSolver::Factor()
{
	return true;
}

//-----------------------------------------------------------------------------
bool RCICGSolver::BackSolve(vector<double>& x, vector<double>& b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// get number of equations
	MKL_INT n = m_pA->Size();

	// zero solution vector
	zero(x);

	// get pointers to solution and RHS vector
	double* px = &x[0];
	double* pb = &b[0];

	// output parameters
	MKL_INT rci_request;
	MKL_INT ipar[128];
	double dpar[128];
	vector<double> tmp(n*4);
	double* ptmp = &tmp[0];

	// initialize parameters
	dcg_init(&n, px, pb, &rci_request, ipar, dpar, ptmp);
	if (rci_request != 0) return false;

	// set the desired parameters:
	// - do residual stopping test
	// - do not request for the user defined stopping test
	// - set the relative tolerance to 1.0E-5;
	ipar[8] = 1;
	ipar[9] = 0;
	dpar[0] = 1e-5;

	// check the consistency of the newly set parameters
	dcg_check(&n, px, pb, &rci_request, ipar, dpar, ptmp);
	if (rci_request != 0) return false;

	// loop until converged
	bool bsuccess = false;
	bool bdone = false;
	do
	{
		// compute the solution by RCI
		dcg(&n, px, pb, &rci_request, ipar, dpar, ptmp);

		switch (rci_request)
		{
		case 0: // solution converged! 
			bsuccess = true;
			bdone = true;
			break;
		case 1: // compute vector A*tmp[0] and store in tmp[n]
			{
				// NOTE: It seems that this blas operation has a memory leak for large problems (+1,500,000). 
				//       The solution is to set the environment variable MKL_DISABLE_FAST_MM to 1
				char tr = 'u';
				double* a = m_pA->Values();
				int* ia = m_pA->Pointers();
				int* ja = m_pA->Indices();
				mkl_dcsrsymv(&tr, &n, a, ia, ja, ptmp, ptmp+n);
			}
			break;
		default:
			bsuccess = false;
			bdone = true;
			break;
		}
	}
	while (!bdone);

	// get convergence information
	int niter;
	dcg_get(&n, px, pb, &rci_request, ipar, dpar, ptmp, &niter);

	// release internal MKL buffers
	MKL_FreeBuffers();

	return bsuccess;
#else
	return false;
#endif // MKL_ISS
}

//-----------------------------------------------------------------------------
void RCICGSolver::Destroy()
{
}
