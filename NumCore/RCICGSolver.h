// This class implements an interface to the RCI CG iterative solver from the MKL math library.
#pragma once

#include "FECore/LinearSolver.h"
#include "CompactSymmMatrix.h"

class RCICGSolver : public LinearSolver
{
public:
	RCICGSolver();
	virtual bool PreProcess();
	virtual bool Factor();
	virtual bool BackSolve(vector<double>& x, vector<double>& b);
	virtual void Destroy();

	virtual SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

	void SetMaxIterations(int n) { m_maxiter = n; }
	void SetTolerance(double tol) { m_tol = tol; }
	void SetPreconditioner(int n) { m_precond = n; }
	void SetPrintLevel(int n) { m_print_level = n; }

private:
	CompactSymmMatrix*	m_pA;
	vector<double>		m_W;	// pre-conditioner

	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_precond;		// pre-conditioner
	int		m_print_level;	// output level
};
