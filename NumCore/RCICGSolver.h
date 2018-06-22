// This class implements an interface to the RCI CG iterative solver from the MKL math library.
#pragma once

#include "FECore/LinearSolver.h"
#include "CompactSymmMatrix.h"
#include "Preconditioner.h"

class RCICGSolver : public LinearSolver
{
public:
	RCICGSolver();
	virtual bool PreProcess();
	virtual bool Factor();
	virtual bool BackSolve(vector<double>& x, vector<double>& b);
	virtual void Destroy();

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	void SetSparseMatrix(SparseMatrix* A) override;

	void SetPreconditioner(Preconditioner* P);

	void SetMaxIterations(int n) { m_maxiter = n; }
	void SetTolerance(double tol) { m_tol = tol; }
	void SetPrintLevel(int n) { m_print_level = n; }

	bool Solve(SparseMatrix* A, vector<double>& x, vector<double>& b, Preconditioner* P = 0);

private:
	SparseMatrix*		m_pA;
	Preconditioner*		m_P;

	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level
};
