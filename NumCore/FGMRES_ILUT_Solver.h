#pragma once
#include "FGMRESSolver.h"
#include <FECore/Preconditioner.h>
#include "ILUT_Preconditioner.h"

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver with
//! ILUT pre-conditioner for nonsymmetric indefinite matrices.
class FGMRES_ILUT_Solver : public FGMRESSolver
{
public:
	//! constructor
	FGMRES_ILUT_Solver(FEModel* fem);

	//! Return a sparse matrix compatible with this solver
	//! This is overridden because this solver will only work with nonsymmetric matrices
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

public: // Preconditioner properties

	// set the max fill value
	void SetMaxFill(int n);

	// Set the fill tolerance
	void SetFillTolerance(double fillTol);

	// do the zero diagonal check during preconditioner
	void DoZeroDiagonalCheck(bool b);

	// Set the zero diagonal tolerance value
	void SetZeroDiagonalTolerance(double tol);

	// set the zero diagonal replacement value
	void SetZeroDiagonalReplacement(double val);

private:
	ILUT_Preconditioner*	m_PC;		//!< the preconditioner
};
