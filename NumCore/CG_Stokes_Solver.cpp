#include "stdafx.h"
#include "CG_Stokes_Solver.h"
#include "StokesPreconditioner.h"

CG_Stokes_Solver::CG_Stokes_Solver()
{
}

//! Set the partition
void CG_Stokes_Solver::SetPartitions(const vector<int>& part)
{
	m_part = part;
}


//! create a sparse matrix that can be used with this solver (must be overridden)
SparseMatrix* CG_Stokes_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (m_part.size() != 2) return 0;

	// create block matrix
	BlockMatrix* A = new BlockMatrix;
	A->Partition(m_part);
	m_pA = A;

	// creat the stokes preconditioner
	StokesPreconditioner* P = new StokesPreconditioner;
	SetPreconditioner(P);

	return m_pA;
}
