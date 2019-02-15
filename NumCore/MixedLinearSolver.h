#pragma once
#include <FECore/LinearSolver.h>
#include "CompactUnSymmMatrix.h"

class FGMRES_ILU0_Solver;

// Experimental solver class that can switch between direct (Pardiso) and 
// iterative (FGMRES\ILU0) solver.
// Currently, only supports non-symmetric matrices
class MixedLinearSolver : public LinearSolver
{
public:
	enum Strategy {
		DIRECT_SOLVER,
		ITERATIVE_SOLVER
	};

public:
	MixedLinearSolver(FEModel* fem);
	~MixedLinearSolver();
	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* y) override;
	void Destroy() override;

	void SetSolverStrategy(Strategy n);

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
	bool SetSparseMatrix(SparseMatrix* pA) override;

public: // properties for iterative solver
	void SetMaxIterations(int nmax);
	void SetPrintLevel(int n) override;
	void SetRelativeConvergence(double tol);
	void SetAbsoluteConvergence(double tol);

protected:
	LinearSolver* currentSolver() { return m_solver[m_strategy]; }

protected:
	int		m_strategy;	// 0 = direct, 1 = iterative
	CRSSparseMatrix*	m_A;
	LinearSolver*	m_solver[2];
	FGMRES_ILU0_Solver*	m_fgmres;
};
