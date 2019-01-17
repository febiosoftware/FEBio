// This class implements an interface to the RCI CG iterative solver from the MKL math library.
#pragma once

#include <FECore/LinearSolver.h>
#include <FECore/Preconditioner.h>
#include <FECore/CompactSymmMatrix.h>

class RCICGSolver : public IterativeLinearSolver
{
public:
	RCICGSolver(FEModel* fem);
	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* b) override;
	void Destroy() override;

public:
	bool HasPreconditioner() const override;

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	bool SetSparseMatrix(SparseMatrix* A) override;

	void SetPreconditioner(Preconditioner* P) override;

	void SetMaxIterations(int n) { m_maxiter = n; }
	void SetTolerance(double tol) { m_tol = tol; }
	void SetPrintLevel(int n) override { m_print_level = n; }

protected:
	SparseMatrix*		m_pA;
	Preconditioner*		m_P;

	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level

	DECLARE_FECORE_CLASS();
};

class RCICG_ICHOL_Solver : public RCICGSolver
{
public:
	RCICG_ICHOL_Solver(FEModel* fem);
	bool Factor() override;
};


class RCICG_Preconditioner : public Preconditioner
{
public:
	RCICG_Preconditioner(FEModel* fem);

	bool Create() override;

	bool mult_vector(double* x, double* y) override;

protected:
	RCICG_ICHOL_Solver	m_cg;
	int	m_neq;
};