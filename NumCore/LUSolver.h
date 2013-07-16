#pragma once
#include "LinearSolver.h"
#include "DenseMatrix.h"
namespace NumCore {

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
	//! Pre-process data
	bool PreProcess();

	//! Factor matrix
	bool Factor();

	//! solve using factored matrix
	bool BackSolve(vector<double>& x, vector<double>& b);

	//! Clean-up
	void Destroy();

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) { return (m_pA = new DenseMatrix()); }

protected:
	vector<int>	indx;	//!< indices
};

} // namespace NumCore
