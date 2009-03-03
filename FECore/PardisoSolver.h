//! This class implements the Pardiso solver.

//! The Pardiso solver is included in the Intel Math Kernel Library (MKL).
//! It can also be installed as a shared object library from
//!		http://www.pardiso-project.org

#pragma once

#include "SparseMatrix.h"
#include "LinearSolver.h"
#include "vector.h"
#include "matrix.h"

	/* Pardiso Fortran prototypes */
extern "C"
{
	int pardisoinit_(void *, int *, int *);

	int pardiso_(void *, int *, int *, int *, int *, int *,
		double *, int *, int *, int *, int *, int *,
		int *, double *, double *, int *);
}

class PardisoSolver : public LinearSolver
{
public:
	virtual bool PreProcess(SparseMatrix& K);
	virtual bool Factor(SparseMatrix& K);
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	virtual void Destroy();

	virtual SparseMatrix* GetMatrix() { return new CompactMatrix(1); }

	PardisoSolver();

protected:
	/* Pardiso control parameters */
	int m_iparm[64];
	int m_maxfct, m_mnum, m_error, m_msglvl;

	/* Matrix data */
	int m_mtype;

	int m_n, m_nnz, m_nrhs;

	/* Internal solver memory pointer */
	void* m_pt[64];
};

