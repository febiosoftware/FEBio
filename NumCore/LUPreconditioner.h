#pragma once
#include <FECore/Preconditioner.h>
#include "PardisoSolver.h"

//-----------------------------------------------------------------------------
class LUPreconditioner : public Preconditioner
{
public:
	LUPreconditioner(FEModel* fem);

	bool Create(SparseMatrix* A) override;

	bool mult_vector(double* x, double* y) override;

private:
	PardisoSolver	m_solver;
};
