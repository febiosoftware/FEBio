#pragma once

#include "LinearSolver.h"
#include "SkylineMatrix.h"

namespace NumCore {

//-----------------------------------------------------------------------------
//! Implements a linear solver that uses a skyline format

class SkylineSolver : public LinearSolver
{
public:
	//! Preprocess 
	bool PreProcess();

	//! Factor matrix
	bool Factor();

	//! Backsolve the linear system
	bool BackSolve(vector<double>& x, vector<double>& b);

	//! Clean up
	void Destroy();

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) { return (m_pA = (ntype == SPARSE_SYMMETRIC? new SkylineMatrix() : 0)); }
};

} // namespace NumCore
