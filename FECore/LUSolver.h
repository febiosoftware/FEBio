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
	bool PreProcess();
	bool Factor();
	bool Solve(vector<double>& x, vector<double>& b);
	void Destroy();

	SparseMatrix* CreateSparseMatrix(int ntype) { return (m_pA = new FullMatrix()); }

protected:
	vector<int>	indx;
};
