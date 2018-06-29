#pragma once
#include <FECore/LinearSolver.h>
#include "BlockMatrix.h"
#include "PardisoSolver.h"

//-----------------------------------------------------------------------------
// This class implements a solution strategy for solving a linear system that is structured
// as a Stokes problem. That is, it is 2x2 block symmetric matrix, but the lower diagonal block
// is zero. 
class StokesSolver : public LinearSolver
{
public:
	//! constructor
	StokesSolver();

	//! destructor
	~StokesSolver();

public:
	//! Preprocess 
	bool PreProcess() override;

	//! Factor matrix
	bool Factor() override;

	//! Backsolve the linear system
	bool BackSolve(vector<double>& x, vector<double>& b) override;

	//! Clean up
	void Destroy() override;

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* A) override;

	//! Set the partition
	void SetPartitions(const vector<int>& part) override;

public:
	// Set the relative convergence tolerance
	void SetRelativeTolerance(double tol);

	// get the iteration count
	int GetIterations() const;

	// set the print level
	void SetPrintLevel(int n);

	// set max nr of iterations
	void SetMaxIterations(int n);

	// set convergence tolerance
	void SetConvergenceTolerance(double tol);

private:
	BlockMatrix*	m_pA;		//!< block matrix
	LinearSolver*	m_solver;	//!< solver for solving diagonal block

private:
	double	m_tol;			//!< convergence tolerance
	int		m_maxiter;		//!< max number of iterations
	int		m_iter;			//!< nr of iterations of last solve
	int		m_printLevel;	//!< set print level
	vector<int>		m_npart;	//!< where to partition the matrix
};
