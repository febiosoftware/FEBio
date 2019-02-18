#include "stdafx.h"
#include "SchurPreconditioner.h"

SchurPreconditioner::SchurPreconditioner(FEModel* fem) : Preconditioner(fem), m_solver(fem)
{
	m_nsize = 0;
	m_solver.SetRelativeResidualTolerance(1e-5);
	m_solver.SetMaxIterations(150);
	m_solver.FailOnMaxIterations(true);

	m_solver.SetLinearSolver(SchurSolver::Diagonal_Solver_LU);
	m_solver.SetSchurSolver(SchurSolver::Schur_Solver_PC);
	m_solver.SetSchurPreconditioner(SchurSolver::Schur_PC_NONE);
	m_solver.SetPrintLevel(0);
	m_solver.DoJacobiPreconditioning(false);
}

void SchurPreconditioner::SetMaxIterations(int n)
{
	m_solver.SetMaxIterations(n);
}

void SchurPreconditioner::SetTolerance(double tol)
{
	m_solver.SetRelativeResidualTolerance(tol);
}

void SchurPreconditioner::ZeroDBlock(bool b)
{
	m_solver.ZeroDBlock(b);
}

void SchurPreconditioner::SetLinearSolver(int n)
{
	m_solver.SetLinearSolver(n);
}

void SchurPreconditioner::SetSchurSolver(int n)
{
	m_solver.SetSchurSolver(n);
}

void SchurPreconditioner::SetSchurPreconditioner(int n)
{
	m_solver.SetSchurPreconditioner(n);
}

bool SchurPreconditioner::Create()
{
	SparseMatrix* A = GetSparseMatrix();
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
