#pragma once

#include "FECore/LinearSolver.h"
#include "SkylineMatrix.h"

//-----------------------------------------------------------------------------
//! Implements a linear solver that uses a skyline format

class SkylineSolver : public LinearSolver
{
public:
	//! constructor
	SkylineSolver(FEModel* fem);

	//! Preprocess 
	bool PreProcess() override;

	//! Factor matrix
	bool Factor() override;

	//! Backsolve the linear system
	bool BackSolve(double* x, double* b) override;

	//! Clean up
	void Destroy() override;

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

private:
	SkylineMatrix*	m_pA;
};
