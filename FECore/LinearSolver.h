#pragma once

#include "SparseMatrix.h"
#include <vector>

//-----------------------------------------------------------------------------
//! Different matrix types. This is used when requesting a sparse matrix format
//! from a linear solver. 
//! \sa LinearSolver::CreateSparseMatrix.
enum Matrix_Type {
	REAL_SYMMETRIC,
	REAL_UNSYMMETRIC,
	COMPLEX_SYMMETRIC,
	COMPLEX_UNSYMMETRIC
};

//-----------------------------------------------------------------------------
//! Abstract base class for the linear solver classes. Linear solver classes
//! are derived from this class and must implement the abstract virtual methods.

//! This class assumes that a linear system is solved in two steps. First, the Factor()
//! method factorizes the matrix, and then BackSolve() solves the system for a given 
//! right hand side vector using the previously factored matrix. 

class FECORE_API LinearSolver
{
public:
	//! constructor
	LinearSolver();

	//! destructor
	virtual ~LinearSolver();

	//! create a sparse matrix that can be used with this solver (must be overridden)
	virtual SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) = 0;

	//! Set the sparse matrix
	virtual bool SetSparseMatrix(SparseMatrix* pA);

	//! Perform any preprocessing
	//! This is called after the structure of the stiffness matrix was determined. 
	//! At this point, we know the size of the matrix and its sparsity pattern.
	virtual bool PreProcess();

	//! Factor the matrix (must be overridden)
	//! Iterative solvers can use this function for creating a pre-conditioner.
	virtual bool Factor() = 0;

	//! do a backsolve, i.e. solve for a right-hand side vector b (must be overridden)
	virtual bool BackSolve(std::vector<double>& x, std::vector<double>& b) = 0;

	//! Do any cleanup
	virtual void Destroy();

	//! Used by block solvers do determine the block partition
	//! The partition is where the global matrix will be divided into blocks
	virtual void SetPartition(int nsplit);
	virtual void SetPartitions(const vector<int>& part);

	//! convenience function for solving linear systems
	bool Solve(vector<double>& x, vector<double>& y);
};

//-----------------------------------------------------------------------------
// base class for iterative solvers
class FECORE_API IterativeLinearSolver : public LinearSolver
{
public:
	// constructor
	IterativeLinearSolver(){}

	// return whether the iterative solver has a preconditioner or not
	virtual bool HasPreconditioner() const = 0;
};
