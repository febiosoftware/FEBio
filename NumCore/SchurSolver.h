#pragma once
#include <FECore/LinearSolver.h>
#include "BlockMatrix.h"

//-----------------------------------------------------------------------------
// This class implements a solution strategy for solving a linear system that is structured
// as a 2x2 block matrix. It makes no assumption on the symmetry of the global matrix or its blocks.
class SchurSolver : public LinearSolver
{
public:
	// options for A block solver
	enum A_Solver {
		A_Solver_LU,
		A_Solver_FGMRES,
		A_Solver_FGMRES_ILU0,
		A_Solver_ILU0,
		A_Solver_DIAGONAL,
		A_Solver_HYPRE
	};

	// options for Schur complement solver
	enum Schur_Solver {
		Schur_Solver_FGMRES,
		Schur_Solver_CG,
		Schur_Solver_PC
	};

	// options for Schur complement preconditioner
	enum Schur_PC {
		Schur_PC_NONE,
		Schur_PC_DIAGONAL_MASS,
		Schur_PC_ICHOL_MASS
	};

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

public:
	// Set the relative convergence tolerance
	void SetRelativeResidualTolerance(double tol);
	void SetAbsoluteResidualTolerance(double tol);

	// get the iteration count
	int GetIterations() const;

	// set the print level
	void SetPrintLevel(int n) override;

	// set max nr of iterations
	void SetMaxIterations(int n);

	void SetLinearSolver(int n);

	int GetLinearSolver();

	void SetSchurSolver(int n);

	void SetSchurPreconditioner(int n);

	void FailOnMaxIterations(bool b);

	void ZeroDBlock(bool b);

	void DoJacobiPreconditioning(bool b);

protected:
	// allocate solver for A block
	LinearSolver* BuildASolver(int nsolver);

	// allocate Schur complement solver
	IterativeLinearSolver* BuildSchurSolver(int nsolver);

	// allocate the preconditioner for the Schur complement solver
	Preconditioner* BuildSchurPreconditioner(int nopt);

private:
	double	m_reltol;		//!< convergence tolerance
	double	m_abstol;		//!< absolute residual convergence tolerance
	int		m_maxiter;		//!< max number of iterations
	int		m_iter;			//!< nr of iterations of last solve
	int		m_printLevel;	//!< set print level
	int		m_nAsolver;		//!< A block solver: 0 = PARDISO, 1 = FGMRES+ILU0, 2 = HYPRE (FGMRES+AMG)
	int		m_nSchurSolver;	//!< Schur solver: 0 = FGMRES, 1 = CG
	int		m_nSchurPreC;	//!< Schur preconditioner : 0 = none
	bool	m_bfailMaxIters;
	bool	m_bzeroDBlock;

private:
	BlockMatrix*	m_pK;					//!< block matrix
	LinearSolver*	m_Asolver;				//!< solver for solving diagonal A block
	IterativeLinearSolver*	m_schurSolver;	//!< solver of Schur complement
	Preconditioner*	m_PS;					//!< preconditioner for the Schur system

	bool			m_doJacobi;		//!< apply Jacobi preconditioner to global system
	vector<double>	m_Wu, m_Wp;		//!< inverse of diagonals of global system (used by Jacobi preconditioner)
};


//-----------------------------------------------------------------------------
// This is just a dummy solver class that multiplies the rhs with the preconditioner
class PCSolver : public IterativeLinearSolver
{
public:	
	PCSolver(FEModel* fem);

	void SetPreconditioner(Preconditioner* pc) override;

	bool PreProcess() override;

	bool Factor() override;

	bool BackSolve(double* x, double* b) override;

	bool HasPreconditioner() const override;

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	bool SetSparseMatrix(SparseMatrix* A) override;

private:
	Preconditioner*	m_PC;
	int		m_neq;
};
