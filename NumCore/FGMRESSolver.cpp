#include "stdafx.h"
#include "FGMRESSolver.h"

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
FGMRESSolver::FGMRESSolver() : m_pA(0)
{
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRESSolver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	if (ntype != REAL_UNSYMMETRIC) return 0;
	m_pA = new CompactUnSymmMatrix(1);
	return m_pA;
#else
	return 0;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::PreProcess()
{
	return true;
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::Factor()
{
	return true;
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::BackSolve(vector<double>& x, vector<double>& b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// number of equations
	MKL_INT N = m_pA->Size();

	// data allocation
	MKL_INT ipar[128];
	double dpar[128];
	MKL_INT RCI_request;
	vector<double> tmp(N*(2*N+1)+(N*(N+9))/2+1);
	double* ptmp = &tmp[0];
	MKL_INT ivar = N;

	// initialize the solver
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// Set the desired parameters:
	ipar[8]=1;		// do residual stopping test
	ipar[9]=0;		// do not request for the user defined stopping test
	ipar[11]=1;		// do the check of the norm of the next generated vector automatically
	dpar[0]=1.0E-3;	// set the relative tolerance to 1.0D-3 instead of default value 1.0D-6

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
				char cvar = 'N'; // multiply with unmodified A
				double* pa = m_pA->Values();
				int* ia = m_pA->Pointers();
				int* ja = m_pA->Indices();
				mkl_dcsrgemv(&cvar, &ivar, pa, ia, ja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
			}
			break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// if converged get the solution. 
	if (bconverged)
	{
		MKL_INT itercount;
		dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp, &itercount);
	}

	MKL_Free_Buffers();
	return bconverged;

#else
	return false;
#endif // MKL_ISS
}

//-----------------------------------------------------------------------------
void FGMRESSolver::Destroy()
{
}
