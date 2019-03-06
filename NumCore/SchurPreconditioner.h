#pragma once
#include <FECore/Preconditioner.h>
#include "SchurSolver.h"

// A preconditioner based on the Schur solver
class SchurPreconditioner : public Preconditioner
{
public:
	SchurPreconditioner(FEModel* fem);

	bool Create() override;

	// apply to vector P x = y
	bool mult_vector(double* x, double* y) override;

	void ZeroDBlock(bool b);

public: // solution strategies
	void SetLinearSolver(int n);
	void SetSchurSolver(int n);
	void SetSchurASolver(int n);
	void SetSchurBlock(int n);
	void SetSchurPreconditioner(int n);

	int GetLinearSolver();

	// Schur parameters
	void SetTolerance(double tol);
	void SetMaxIterations(int n);

private:
	SchurSolver		m_solver;
	int				m_nsize;
};
