#include "stdafx.h"
#include "LUPreconditioner.h"

// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef MKL_ISS
#ifdef PARDISO
#undef PARDISO
#endif
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#endif // MKL_ISS

LUPreconditioner::LUPreconditioner(FEModel* fem) : Preconditioner(fem), m_solver(fem)
{

}

bool LUPreconditioner::Create(SparseMatrix* A)
{
	// make sure we have work to do
	if ((A == nullptr) || (A->Rows() == 0)) return false;
	CompactMatrix* K = dynamic_cast<CompactMatrix*>(A);
	if (K == nullptr) return false;

	m_solver.SetSparseMatrix(K);
	if (m_solver.PreProcess() == false) return false;
	if (m_solver.Factor() == false) return false;

	return true;
}

bool LUPreconditioner::mult_vector(double* x, double* y)
{
	return m_solver.BackSolve(x, y);
}
