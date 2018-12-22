#pragma once
#include <FECore/LinearSolver.h>
#include <FECore/CompactUnSymmMatrix.h>

//-----------------------------------------------------------------------------
// This class implements the HYPRE GMRES solver
class HypreGMRESsolver : public LinearSolver
{
	class Implementation;

public:
	HypreGMRESsolver(FEModel* fem);
	~HypreGMRESsolver();

	void SetPrintLevel(int n) override;

	void SetMaxIterations(int n);

	void SetConvergencTolerance(double tol);

public:
	// allocate storage
	bool PreProcess() override;

	//! Factor the matrix (for iterative solvers, this can be used for creating pre-conditioner)
	bool Factor() override;

	//! Calculate the solution of RHS b and store solution in x
	bool BackSolve(double* x, double* b) override;

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* A) override;

private:
	Implementation*	imp;
};
