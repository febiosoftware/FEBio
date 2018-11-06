#include "stdafx.h"
#include "SchurPreconditioner.h"

SchurPreconditioner::SchurPreconditioner(FEModel* fem) : Preconditioner(fem), m_solver(fem)
{
	m_nsize = 0;
	m_solver.UseMassMatrix(true);
	m_solver.SetConvergenceTolerance(1e-5);
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
	return m_solver.SetSparseMatrix(A);
}

// apply to vector P x = y
bool SchurPreconditioner::mult_vector(double* x, double* y)
{
	int neq = m_nsize;
	// copy x vector
	vector<double> b(neq), tmp(neq, 0);
	for (int i = 0; i<neq; ++i) b[i] = x[i];
	if (m_solver.Solve(tmp, b) == false) return false;
	// copy solution
	for (int i = 0; i<neq; ++i) y[i] = tmp[i];
	return true;
}
