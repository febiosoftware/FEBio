#pragma once
#include <FECore/LinearSolver.h>
#include "BlockMatrix.h"
#include "PardisoSolver.h"

//-----------------------------------------------------------------------------
// This class implements solution strategies for solving linear systems by taking
// advantage of their block structure.
class BlockSolver : public LinearSolver
{
public:
	//! constructor
	BlockSolver(FEModel* fem);

	//! destructor
	~BlockSolver();

	//! Preprocess 
	bool PreProcess() override;

	//! Factor matrix
	bool Factor() override;

	//! Backsolve the linear system
	bool BackSolve(double* x, double* y) override;

	//! Clean up
	void Destroy() override;

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

public:
	// Set the relative convergence tolerance
	void SetRelativeTolerance(double tol);

	// get the iteration count
	int GetIterations() const;

	// set the print level
	void SetPrintLevel(int n) override;

private:
	BlockMatrix*			m_pA;		//!< block matrices
	vector<PardisoSolver*>	m_solver;	//!< solvers for solving diagonal blocks

private:
	double	m_tol;			//!< convergence tolerance
	int		m_maxiter;		//!< max number of iterations
	int		m_iter;			//!< nr of iterations of last solve
	int		m_printLevel;	//!< set print level
};
