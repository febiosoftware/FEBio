#include "stdafx.h"
#include "FGMRES_Schur_Solver.h"

//=============================================================================
FGMRES_Schur_Solver::FGMRES_Schur_Solver(FEModel* fem) : FGMRESSolver(fem)
{
	SetPreconditioner(m_PC = new SchurPreconditioner(fem));
}

//-----------------------------------------------------------------------------
//! Set the partition
void FGMRES_Schur_Solver::SetPartitions(const vector<int>& part)
{
	m_npart = part;
}

//-----------------------------------------------------------------------------
void FGMRES_Schur_Solver::ZeroDBlock(bool b)
{
	m_PC->ZeroDBlock(b);
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_Schur_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (m_npart.size() != 2) return 0;
	BlockMatrix* pA = new BlockMatrix();
	pA->Partition(m_npart, ntype);
	FGMRESSolver::SetSparseMatrix(pA);
	return pA;
}

//-----------------------------------------------------------------------------
// this is used to build the preconditioner
bool FGMRES_Schur_Solver::Factor()
{
	if (FGMRESSolver::Factor() == false) return false;
//	m_PC->SetMaxIterations(GetMaxIterations());
	return m_PC->Create(GetSparseMatrix());
}
