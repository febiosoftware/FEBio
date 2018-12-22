#pragma once
#include "FGMRESSolver.h"
#include <FECore/Preconditioner.h>
#include "ILU0_Preconditioner.h"

class CRSSparseMatrix;

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver with
//! ILU0 pre-conditioner for nonsymmetric indefinite matrices.
class FGMRES_ILU0_Solver : public FGMRESSolver
{
public:
	//! constructor
	FGMRES_ILU0_Solver(FEModel* fem);

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

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

//-----------------------------------------------------------------------------
class ILU0_Solver : public LinearSolver
{
public:
	ILU0_Solver(FEModel* fem);
	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* y) override;

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
	bool SetSparseMatrix(SparseMatrix* pA) override;

private:
	Preconditioner*		m_PC;
	CRSSparseMatrix*	m_A;
};
