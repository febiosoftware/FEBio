#pragma once
#include "FECore/LinearSolver.h"
#include "CompactSymmMatrix.h"
#include <vector>

//-----------------------------------------------------------------------------
//! this class implements an iterative conjugate gradient solver 
class ConjGradIterSolver : public LinearSolver
{
public:
	//! constructor
	ConjGradIterSolver(FEModel* fem);

	//! Pre-process data
	bool PreProcess();

	//! Factor matrix
	bool Factor();

	//! solve system given a rhs vector
	bool BackSolve(double* x, double* y);

	//! Clean up
	void Destroy();

	//! Create a sparse matrix for this linear solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

public:
	CompactSymmMatrix*	m_pA;

	double	m_tol;		//!< convergence tolerance
	int		m_kmax;		//!< max iterations
	int		m_nprint;	//!< printing level

	std::vector<double>	m_P;	//!< preconditioning vector
};
