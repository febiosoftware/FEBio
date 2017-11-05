#include "stdafx.h"
#include "RCICGSolver.h"

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
RCICGSolver::RCICGSolver() : m_pA(0)
{
	m_maxiter = 0;
	m_tol = 1e-5;
	m_precond = 1;
	m_print_level = 0;
}

//-----------------------------------------------------------------------------
SparseMatrix* RCICGSolver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	if (ntype != REAL_SYMMETRIC) return 0;
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
	if (m_pA == 0) return false;

	// get number of equations
	MKL_INT n = m_pA->Size();

	if (m_precond == 1)
	{
		// calculate pre-conditioner
		m_W.assign(n, 0.0);
		for (int i=0; i<n; ++i)
		{
			double di = m_pA->diag(i);
			if (di == 0.0) di = 1.0;
			m_W[i] = 1.0 / di;
		}
	}
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
	if (m_maxiter > 0) ipar[4] = m_maxiter;	// max nr of iterations
	ipar[8] = 1;			// do residual stopping test
	ipar[9] = 0;			// do not request for the user defined stopping test
	ipar[10] = m_precond;	// do preconditioning
	dpar[0] = m_tol;		// set the relative tolerance

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
		case 3:
			{
				assert(m_precond != 0);
				for (int i = 0; i<n; ++i) ptmp[3 * n + i] = m_W[i] * ptmp[2 * n + i];
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

	if (m_print_level != 0)
	{
		printf("CG iterations = %d\n", niter);
	}

	// release internal MKL buffers
	MKL_Free_Buffers();

	return bsuccess;
#else
	return false;
#endif // MKL_ISS
}

//-----------------------------------------------------------------------------
void RCICGSolver::Destroy()
{
}
