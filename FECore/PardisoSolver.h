//! This class implements the Pardiso solver.

//! The Pardiso solver is included in the Intel Math Kernel Library (MKL).
//! It can also be installed as a shared object library from
//!		http://www.pardiso-project.org

#pragma once

#include "LinearSolver.h"

	/* Pardiso Fortran prototypes */
#ifndef WIN32
extern "C"
{
	int pardisoinit_(void *, int *, int *);

	int pardiso_(void *, int *, int *, int *, int *, int *,
		double *, int *, int *, int *, int *, int *,
		int *, double *, double *, int *);
}
#endif

class PardisoSolver : public LinearSolver
{
public:
	bool PreProcess(SparseMatrix& K);
	bool Factor(SparseMatrix& K);
	bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	void Destroy();

	SparseMatrix* GetMatrix(int ntype)
	{
		m_bsymm = (ntype == SPARSE_SYMMETRIC);
		if (m_bsymm) return new CompactSymmMatrix(1);
		else return new CompactUnSymmMatrix(1, true);
	}

	PardisoSolver();

protected:

	bool m_bsymm; // use symmetric mode or not

	// Pardiso control parameters
	int m_iparm[64];
	int m_maxfct, m_mnum, m_error, m_msglvl;

	// Matrix data
	int m_mtype;
	int m_n, m_nnz, m_nrhs;

	void* m_pt[64]; // Internal solver memory pointer
};

