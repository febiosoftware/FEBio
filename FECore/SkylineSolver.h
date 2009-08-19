#pragma once

#include "LinearSolver.h"

//-----------------------------------------------------------------------------
//! Implements a linear solver that uses a skyline format

class SkylineSolver : public LinearSolver
{
public:
	bool PreProcess(SparseMatrix& K);
	bool Factor(SparseMatrix& K);
	bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	void Destroy();

	SparseMatrix* GetMatrix(int ntype) { return (ntype == SPARSE_SYMMETRIC? new SkylineMatrix() : 0); }
};
