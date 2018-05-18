#include "stdafx.h"
#include "SkylineSolver.h"

//-----------------------------------------------------------------------------
FECORE_API void colsol_factor(int N, double* values, int* pointers);
FECORE_API void colsol_solve(int N, double* values, int* pointers, double* R);

//-----------------------------------------------------------------------------
SkylineSolver::SkylineSolver() : m_pA(0)
{
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* SkylineSolver::CreateSparseMatrix(Matrix_Type ntype)
{ 
	return (m_pA = (ntype == REAL_SYMMETRIC? new SkylineMatrix() : 0)); 
}

//-----------------------------------------------------------------------------
bool SkylineSolver::PreProcess()
{
	// We don't need to do any preprocessing for this solver
	return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool SkylineSolver::Factor()
{
	colsol_factor(m_pA->Size(), m_pA->values(), m_pA->pointers());
	return true;
}

//-----------------------------------------------------------------------------
bool SkylineSolver::BackSolve(vector<double>& x, vector<double>& R)
{
	// we need to make a copy of R since colsol overwrites the right hand side vector
	// with the solution
	int neq = m_pA->Size();
	for (int i=0; i<neq; ++i) x[i] = R[i];
	colsol_solve(m_pA->Size(), m_pA->values(), m_pA->pointers(), &x[0]);

	return true;
}

//-----------------------------------------------------------------------------
void SkylineSolver::Destroy()
{
	// Nothing to destroy
	LinearSolver::Destroy();
}
