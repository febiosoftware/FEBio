#include "stdafx.h"
#include "FGMRES_ILU0_Solver.h"
#include <FECore/CompactUnSymmMatrix.h>
#include "ILU0_Preconditioner.h"

//=============================================================================
FGMRES_ILU0_Solver::FGMRES_ILU0_Solver(FEModel* fem) : FGMRESSolver(fem)
{
	SetPreconditioner(m_PC = new ILU0_Preconditioner(fem));
}

//-----------------------------------------------------------------------------
// do the zero diagonal check during preconditioner
void FGMRES_ILU0_Solver::DoZeroDiagonalCheck(bool b)
{
	m_PC->m_checkZeroDiagonal = b;
}

//-----------------------------------------------------------------------------
// Set the zero diagonal tolerance value
void FGMRES_ILU0_Solver::SetZeroDiagonalTolerance(double tol)
{
	m_PC->m_zeroThreshold = tol;
}

//-----------------------------------------------------------------------------
// set the zero diagonal replacement value
void FGMRES_ILU0_Solver::SetZeroDiagonalReplacement(double val)
{
	m_PC->m_zeroReplace = val;
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_ILU0_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	// We can only support non-symmetric matrices on account of the preconditioner
	if (ntype != REAL_UNSYMMETRIC) return 0;
	SparseMatrix* A = FGMRESSolver::CreateSparseMatrix(ntype);
	m_PC->SetSparseMatrix(A);
	return A;
}

//-----------------------------------------------------------------------------
// this is used to build the preconditioner
bool FGMRES_ILU0_Solver::Factor()
{
	if (FGMRESSolver::Factor() == false) return false;
	return m_PC->Create();
}

//==================================================================================
ILU0_Solver::ILU0_Solver(FEModel* fem) : LinearSolver(fem)
{
	ILU0_Preconditioner* PC = new ILU0_Preconditioner(fem);
	m_PC = PC;
}

bool ILU0_Solver::PreProcess() { return true; }
bool ILU0_Solver::Factor() { return m_PC->Create(); }
bool ILU0_Solver::BackSolve(double* x, double* y)
{
	return m_PC->mult_vector(y, x);
}

SparseMatrix* ILU0_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (ntype == REAL_SYMMETRIC) return nullptr;
	CRSSparseMatrix* A = new CRSSparseMatrix(1);
	m_PC->SetSparseMatrix(A);
	return A;
}

bool ILU0_Solver::SetSparseMatrix(SparseMatrix* pA)
{
	CRSSparseMatrix* A = dynamic_cast<CRSSparseMatrix*>(pA);
	if ((A == nullptr) || (A->Offset() != 1)) return false;

	m_PC->SetSparseMatrix(A);
	return true;
}
