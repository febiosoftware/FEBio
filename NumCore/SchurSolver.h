#pragma once
#include <FECore/LinearSolver.h>
#include <FECore/BlockMatrix.h>

//-----------------------------------------------------------------------------
// This class implements a solution strategy for solving a linear system that is structured
// as a 2x2 block matrix. It makes no assumption on the symmetry of the global matrix or its blocks.
class SchurSolver : public LinearSolver
{
public:
	//! constructor
	SchurSolver(FEModel* fem);

	//! destructor
	~SchurSolver();

public:
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

	//! set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* A) override;

	//! Set the partition
	void SetPartitions(const vector<int>& part) override;

public:
	// Set the relative convergence tolerance
	void SetRelativeResidualTolerance(double tol);
	void SetAbsoluteResidualTolerance(double tol);

	// get the iteration count
	int GetIterations() const;

	// set the print level
	void SetPrintLevel(int n);

	// set max nr of iterations
	void SetMaxIterations(int n);

	// Set the schur block
	void SetSchurBlock(int n);

	bool BuildMassMatrix(CompactSymmMatrix* M);

	void UseMassMatrix(bool b);

	void SetLinearSolver(int n);

	void FailOnMaxIterations(bool b);

	void ZeroDBlock(bool b);

private:
	BlockMatrix*	m_pA;		//!< block matrix
	LinearSolver*	m_solver;	//!< solver for solving diagonal block
	IterativeLinearSolver*	m_schurSolver;	//!< solver of Schur complement
	Preconditioner*	m_PS;		//!< preconditioner for the Schur system

private:
	double	m_reltol;		//!< convergence tolerance
	double	m_abstol;		//!< absolute residual convergence tolerance
	int		m_maxiter;		//!< max number of iterations
	int		m_iter;			//!< nr of iterations of last solve
	int		m_printLevel;	//!< set print level
	int		m_nsolver;		//!< 0 = FGMRES+ILU0, 1 = HYPRE (FGMRES+AMG)
	int		m_schurBlock;	//!< the block that will be used to calculated schur complement
	bool	m_buildMassMatrix;
	bool	m_bfailMaxIters;
	bool	m_bzeroDBlock;

	vector<int>		m_npart;	//!< where to partition the matrix
};
