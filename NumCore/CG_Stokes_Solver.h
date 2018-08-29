#pragma once
#include "RCICGSolver.h"
#include "BlockMatrix.h"

// Implements a CG solver with Stokes preconditioner
class CG_Stokes_Solver : public RCICGSolver
{
public:
	CG_Stokes_Solver();

	//! create a sparse matrix that can be used with this solver (must be overridden)
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! Set the partition
	void SetPartitions(const vector<int>& part) override;

private:
	vector<int>		m_part;
};
