#include "stdafx.h"
#include "SchurPreconditioner.h"

SchurPreconditioner::SchurPreconditioner(FEModel* fem) : Preconditioner(fem), m_solver(fem)
{
	m_nsize = 0;
	m_solver.UseMassMatrix(true);
	m_solver.SetRelativeResidualTolerance(1e-5);
	m_solver.SetMaxIterations(250);
	m_solver.FailOnMaxIterations(false);
}

void SchurPreconditioner::SetMaxIterations(int n)
{
	m_solver.SetMaxIterations(n);
}

void SchurPreconditioner::ZeroDBlock(bool b)
{
	m_solver.ZeroDBlock(b);
}

bool SchurPreconditioner::Create(SparseMatrix* A)
{
	m_nsize = A->Rows();
	if (m_solver.SetSparseMatrix(A) == false) return false;
	if (m_solver.PreProcess() == false) return false;
	if (m_solver.Factor() == false) return false;
	return true;
}

// apply to vector P x = y
bool SchurPreconditioner::mult_vector(double* x, double* y)
{
	return m_solver.BackSolve(y, x);
}
