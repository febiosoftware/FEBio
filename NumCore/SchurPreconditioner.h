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

	void SetMaxIterations(int n);

	void ZeroDBlock(bool b);

public: // solution strategies
	void SetLinearSolver(int n);
	void SetSchurSolver(int n);
	void SetSchurPreconditioner(int n);

private:
	SchurSolver		m_solver;
	int				m_nsize;
};
