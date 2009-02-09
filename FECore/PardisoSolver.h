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

	virtual SparseMatrix* GetMatrix() { return new CompactMatrix(true); }

	PardisoSolver();

protected:
	/* Pardiso control parameters */
	int iparm[64];
	int maxfct, mnum, phase, error, msglvl;

	/* Matrix data */
	int mtype;

	int n, nnz, nrhs;

	/* Internal solver memory pointer */
	void* pt[64];

	/* Number of processors */
	int num_procs;

	/* Auxiliary variables */
	char* var;
	int i;
	double ddum; /* Double dummy */
	int idum; /* Integer dummy */
};

