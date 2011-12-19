#pragma once

#include "LinearSolver.h"
#include "SkylineMatrix.h"

namespace NumCore {

//-----------------------------------------------------------------------------
//! Implements a linear solver that uses a skyline format

class SkylineSolver : public LinearSolver
{
public:
	bool PreProcess();
	bool Factor();
	bool BackSolve(vector<double>& x, vector<double>& b);
	void Destroy();

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) { return (m_pA = (ntype == SPARSE_SYMMETRIC? new SkylineMatrix() : 0)); }
};

} // namespace NumCore
