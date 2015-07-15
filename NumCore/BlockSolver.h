#pragma once
#include "FECore/LinearSolver.h"
#include "BlockMatrix.h"
#include "PardisoSolver.h"

//-----------------------------------------------------------------------------
// This class implements solution strategies for solving linear systems by taking
// advantage of their block structure.
class BlockSolver : public LinearSolver
{
public:
	//! constructor
	BlockSolver();

	//! destructor
	~BlockSolver();

	//! Preprocess 
	bool PreProcess();

	//! Factor matrix
	bool Factor();

	//! Backsolve the linear system
	bool BackSolve(vector<double>& x, vector<double>& b);

	//! Clean up
	void Destroy();

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

private:
	BlockMatrix*			m_pA;		//!< block matrices
	vector<PardisoSolver*>	m_solver;	//!< solvers for solving diagonal blocks
};
