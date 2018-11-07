#include "stdafx.h"
#include "FGMRES_ILU0_Solver.h"
#include <FECore/CompactUnSymmMatrix.h>
#include "ILUT_Preconditioner.h"

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
	return FGMRESSolver::CreateSparseMatrix(ntype);
}

//-----------------------------------------------------------------------------
// this is used to build the preconditioner
bool FGMRES_ILU0_Solver::Factor()
{
	if (FGMRESSolver::Factor() == false) return false;
	return m_PC->Create(GetSparseMatrix());
}

//==================================================================================
ILU0_Solver::ILU0_Solver(FEModel* fem) : LinearSolver(fem), m_A(nullptr) 
{
	ILU0_Preconditioner* PC = new ILU0_Preconditioner(fem);
	m_PC = PC;
}

bool ILU0_Solver::PreProcess() { return true; }
bool ILU0_Solver::Factor() { return m_PC->Create(m_A); }
bool ILU0_Solver::BackSolve(double* x, double* y)
{
	return m_A->mult_vector(x, y);
}

SparseMatrix* ILU0_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (ntype != REAL_UNSYMMETRIC) return nullptr;
	m_A = new CRSSparseMatrix(1);
	return m_A;
}
bool ILU0_Solver::SetSparseMatrix(SparseMatrix* pA)
{
	m_A = dynamic_cast<CRSSparseMatrix*>(pA);
	return (m_A != nullptr);
}
