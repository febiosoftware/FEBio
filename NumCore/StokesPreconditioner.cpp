#include "stdafx.h"
#include "StokesPreconditioner.h"

StokesPreconditioner::StokesPreconditioner(FEModel* fem) : m_solver(fem)
{
	m_nsize = 0;
}

bool StokesPreconditioner::Create(SparseMatrix* A)
{
	m_nsize = A->Rows();
	return m_solver.SetSparseMatrix(A);
}

// apply to vector P x = y
void StokesPreconditioner::mult_vector(double* x, double* y)
{
	int neq = m_nsize;
	// copy x vector
	vector<double> b(neq), tmp(neq, 0);
	for (int i=0; i<neq; ++i) b[i] = x[i];
	m_solver.Solve(tmp, b);
	// copy solution
	for (int i=0; i<neq; ++i) y[i] = tmp[i];
}
