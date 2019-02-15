#pragma once

#include "SparseMatrix.h"
#include "FECoreBase.h"
#include "fecore_enum.h"
#include <vector>

class FEModel;
class Preconditioner;

//-----------------------------------------------------------------------------
struct LinearSolverStats
{
	int		backsolves;		// number of times backsolve was called
	int		iterations;		// total number of iterations
};

//-----------------------------------------------------------------------------
//! Abstract base class for the linear solver classes. Linear solver classes
//! are derived from this class and must implement the abstract virtual methods.

//! This class assumes that a linear system is solved in two steps. First, the Factor()
//! method factorizes the matrix, and then BackSolve() solves the system for a given 
//! right hand side vector using the previously factored matrix. 

class FECORE_API LinearSolver : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	//! constructor
	LinearSolver(FEModel* fem);

	//! destructor
	virtual ~LinearSolver();

	virtual void SetPrintLevel(int n) {}

public:

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

	//! do a backsolve, i.e. solve for a right-hand side vector y (must be overridden)
	virtual bool BackSolve(double* x, double* y) = 0;

	//! Do any cleanup
	virtual void Destroy();

	//! Used by block solvers do determine the block partition
	//! The partition is where the global matrix will be divided into blocks
	void SetPartitions(const vector<int>& part);
	void SetPartitions(int npart0, int npart1);

	//! version for std::vector
	bool BackSolve(std::vector<double>& x, std::vector<double>& b)
	{
		return BackSolve(&x[0], &b[0]);
	}

	//! convenience function for solving linear systems
	bool Solve(vector<double>& x, vector<double>& y);

public:
	LinearSolverStats GetStats() const;

	void ResetStats();

protected:
	// used by derived classes to update stats.
	// Should be called after each backsolve. Will increment backsolves by one and add iterations
	void UpdateStats(int iterations);

protected:
	std::vector<int>	m_part;		//!< partitions of linear system.

private:
	LinearSolverStats	m_stats;	//!< stats on how often linear solver was called.
};

//-----------------------------------------------------------------------------
// base class for iterative solvers
class FECORE_API IterativeLinearSolver : public LinearSolver
{
public:
	// constructor
	IterativeLinearSolver(FEModel* fem) : LinearSolver(fem) {}

	// return whether the iterative solver has a preconditioner or not
	virtual bool HasPreconditioner() const = 0;

	// set the preconditioner
	virtual void SetPreconditioner(Preconditioner* pc) {}

	// get the preconditioner
	virtual Preconditioner* GetPreconditioner() { return nullptr; }

public:
	// helper function for solving a linear system of equations
	bool Solve(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b, Preconditioner* pc = 0);
};
