#include "stdafx.h"
#include "FGMRES_Jacobi_Solver.h"

FGMRES_Jacobi_Solver::FGMRES_Jacobi_Solver(FEModel* fem) : FGMRESSolver(fem)
{
	m_W = new DiagonalPreconditioner(fem);
	m_W->CalculateSquareRoot(true);

	// Set the left and right preconditioner
	SetPreconditioner(m_W);
	SetRightPreconditioner(m_W);
}

FGMRES_Jacobi_Solver::~FGMRES_Jacobi_Solver()
{
	delete m_W;
}

//! Return a sparse matrix compatible with this solver
SparseMatrix* FGMRES_Jacobi_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	SparseMatrix* A = FGMRESSolver::CreateSparseMatrix(ntype);
	if (A) m_W->SetSparseMatrix(A);
	return A;
}

// this is used to build the preconditioner
bool FGMRES_Jacobi_Solver::Factor()
{
	if (m_W->Create() == false) return false;
	return FGMRESSolver::Factor();
}
