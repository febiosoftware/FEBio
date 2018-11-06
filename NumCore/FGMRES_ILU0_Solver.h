#pragma once
#include "FGMRESSolver.h"
#include <FECore/Preconditioner.h>
#include "ILU0_Preconditioner.h"

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver with
//! ILU0 pre-conditioner for nonsymmetric indefinite matrices.
class FGMRES_ILU0_Solver : public FGMRESSolver
{
public:
	//! constructor
	FGMRES_ILU0_Solver(FEModel* fem);

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

	// this is used to build the preconditioner
	bool Factor() override;

public: // preconditioner settings

	// do the zero diagonal check during preconditioner
	void DoZeroDiagonalCheck(bool b);

	// Set the zero diagonal tolerance value
	void SetZeroDiagonalTolerance(double tol);

	// set the zero diagonal replacement value
	void SetZeroDiagonalReplacement(double val);

private:
	ILU0_Preconditioner*	m_PC;		//!< the preconditioner
};
