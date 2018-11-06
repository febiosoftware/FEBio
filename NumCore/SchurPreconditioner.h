#pragma once
#include <FECore/Preconditioner.h>
#include "SchurSolver.h"

// A preconditioner based on the Schur solver
class SchurPreconditioner : public Preconditioner
{
public:
	SchurPreconditioner(FEModel* fem);

	bool Create(SparseMatrix* A) override;

	// apply to vector P x = y
	bool mult_vector(double* x, double* y) override;

	void SetMaxIterations(int n);

	void ZeroDBlock(bool b);

private:
	SchurSolver		m_solver;
	int				m_nsize;
};
