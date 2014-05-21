#pragma once

#include "SparseMatrix.h"
#include <vector>

//-----------------------------------------------------------------------------
//! Different matrix types. This is used when requesting a sparse matrix format
//! from a linear solver. 
//! \sa LinearSolver::CreateSparseMatrix.
enum Matrix_Type {
	SPARSE_SYMMETRIC,
	SPARSE_UNSYMMETRIC
};

//-----------------------------------------------------------------------------
//! Abstract base class for the linear solver classes. Linear solver classes
//! are derived from this class and must implement the abstract virtual methods.

//! This class assumes that a linear system is solved in two steps. First, the Factor()
//! method factorizes the matrix, and then BackSolve() solves the system for a given 
//! right hand side vector using the previously factored matrix. 

class LinearSolver
{
public:
	//! constructor
	LinearSolver();

	//! destructor
	virtual ~LinearSolver();

	//! create a sparse matrix that can be used with this solver (must be overridden)
	virtual SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) = 0;

	//! perform any preprocessing
	virtual bool PreProcess();

	//! factor the matrix (must be overridden)
	virtual bool Factor() = 0;

	//! do a backsolve, i.e. solve for a right-hand side vector b (must be overridden)
	virtual bool BackSolve(std::vector<double>& x, std::vector<double>& b) = 0;

	//! Do any cleanup
	virtual void Destroy();
};
