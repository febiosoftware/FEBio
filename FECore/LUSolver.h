#pragma once

#include "LinearSolver.h"

//-----------------------------------------------------------------------------
//! LU decomposition solver

//! This solver performs an LU decomposition and uses a backsolving algorithm
//! to solve the equations.
//! This solver uses the FullMatrix class and therefore is not the preferred
//! solver. It should only be used for small problems and only when the other
//! solvers are not adequate.

class LUSolver : public LinearSolver
{
public:
	bool PreProcess(SparseMatrix& K);
	bool Factor(SparseMatrix& K);
	bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	void Destroy(SparseMatrix& K);

	SparseMatrix* GetMatrix(int ntype) { return new FullMatrix(); }

protected:
	vector<int>	indx;
};
