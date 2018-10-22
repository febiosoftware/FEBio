#include "stdafx.h"
#include "FGMRES_ILUT_Solver.h"

//=============================================================================
FGMRES_ILUT_Solver::FGMRES_ILUT_Solver(FEModel* fem) : FGMRESSolver(fem)
{
	// set the preconditioner
	SetPreconditioner(m_PC = new ILUT_Preconditioner);
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_ILUT_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	// We can only support non-symmetric matrices because of the ILUT preconditioner
	if (ntype != REAL_UNSYMMETRIC) return 0;
	return FGMRESSolver::CreateSparseMatrix(ntype);
}

//-----------------------------------------------------------------------------
// set the max fill value of the preconditioner
void FGMRES_ILUT_Solver::SetMaxFill(int n)
{
	m_PC->m_maxfill = n;
}

//-----------------------------------------------------------------------------
// Set the fill tolerance of the preconditioner
void FGMRES_ILUT_Solver::SetFillTolerance(double fillTol)
{
	m_PC->m_fillTol = fillTol;
}

//-----------------------------------------------------------------------------
// do the zero diagonal check during preconditioner
void FGMRES_ILUT_Solver::DoZeroDiagonalCheck(bool b)
{
	m_PC->m_checkZeroDiagonal = b;
}

//-----------------------------------------------------------------------------
// Set the zero diagonal tolerance value of the preconditioner
void FGMRES_ILUT_Solver::SetZeroDiagonalTolerance(double tol)
{
	m_PC->m_zeroThreshold = tol;
}

//-----------------------------------------------------------------------------
// set the zero diagonal replacement value of the preconditioner
void FGMRES_ILUT_Solver::SetZeroDiagonalReplacement(double val)
{
	m_PC->m_zeroReplace = val;
}
