#pragma once
#include <FECore/LinearSolver.h>
#include "CompactUnSymmMatrix.h"
#include "CompactSymmMatrix.h"

//! The Pardiso solver is included in the Intel Math Kernel Library (MKL).
//! It can also be installed as a shared object library from
//!		http://www.pardiso-project.org


class PardisoSolver : public LinearSolver
{
public:
	PardisoSolver(FEModel* fem);
	~PardisoSolver();
	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* y) override;
	void Destroy() override;

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
	bool SetSparseMatrix(SparseMatrix* pA) override;

	void PrintConditionNumber(bool b);

	double condition_number();

	void UseIterativeFactorization(bool b);

protected:

	CompactMatrix*	m_pA;
	int				m_mtype; // matrix type

	// Pardiso control parameters
	int m_iparm[64];
	int m_maxfct, m_mnum, m_msglvl;
	double m_dparm[64];

	bool m_iparm3;	// use direct-iterative method

	// Matrix data
	int m_n, m_nnz, m_nrhs;

	bool	m_print_cn;	// estimate and print the condition number

	void* m_pt[64]; // Internal solver memory pointer
};
