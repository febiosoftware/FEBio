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

public:
	// Set the relative convergence tolerance
	void SetRelativeTolerance(double tol);

	// get the iteration count
	int GetIterations() const;

	// set the print level
	void SetPrintLevel(int n);

private:
	BlockMatrix*	m_pA;		//!< block matrix
	PardisoSolver*	m_solver;	//!< solver for solving diagonal block

private:
	double	m_tol;			//!< convergence tolerance
	int		m_maxiter;		//!< max number of iterations
	int		m_iter;			//!< nr of iterations of last solve
	int		m_printLevel;	//!< set print level
};
