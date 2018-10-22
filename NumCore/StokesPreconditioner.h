#pragma once
#include "Preconditioner.h"
#include "StokesSolver.h"

// A preconditioner based on the Stokes solver
class StokesPreconditioner : public Preconditioner
{
public:
	StokesPreconditioner(FEModel* fem);

	bool Create(SparseMatrix* A) override;

	// apply to vector P x = y
	void mult_vector(double* x, double* y) override;

private:
	StokesSolver	m_solver;
	int				m_nsize;
};
