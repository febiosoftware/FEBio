#pragma once
#include "FGMRESSolver.h"
#include <FECore/Preconditioner.h>

class FGMRES_Jacobi_Solver : public FGMRESSolver
{
public:
	FGMRES_Jacobi_Solver(FEModel* fem);

	~FGMRES_Jacobi_Solver();

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	// this is used to build the preconditioner
	bool Factor() override;

protected:
	DiagonalPreconditioner*	m_W;
};
