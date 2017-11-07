//! This class implements the Pardiso solver.

//! The Pardiso solver is included in the Intel Math Kernel Library (MKL).
//! It can also be installed as a shared object library from
//!		http://www.pardiso-project.org

#pragma once

#include "FECore/LinearSolver.h"
#include "CompactMatrix.h"

#ifdef PARDISO
	/* Pardiso prototypes for MKL version */
extern "C"
{
	int pardisoinit_(void *, int *, int *);

	int pardiso_(void *, int *, int *, int *, int *, int *,
		double *, int *, int *, int *, int *, int *,
		int *, double *, double *, int *);
}

#else
/* Pardiso prototypes for shared object library version */

#ifdef WIN32

#define pardisoinit_ PARDISOINIT
#define pardiso_ PARDISO

#endif

extern "C"
{
int pardisoinit_(void *, int *, int *, int *, double*, int*);

int pardiso_(void *, int *, int *, int *, int *, int *,
	double *, int *, int *, int *, int *, int *,
	int *, double *, double *, int *, double*);
}
#endif

class PardisoSolver : public LinearSolver
{
public:
	PardisoSolver();
	bool PreProcess();
	bool Factor();
	bool BackSolve(vector<double>& x, vector<double>& b);
	void Destroy();

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);
	void SetSparseMatrix(CompactMatrix* pA);

protected:

	CompactMatrix*	m_pA;
	bool m_bsymm; // use symmetric mode or not

	// Pardiso control parameters
	int m_iparm[64];
	int m_maxfct, m_mnum, m_msglvl;
	double m_dparm[64];

	// Matrix data
	int m_mtype;
	int m_n, m_nnz, m_nrhs;

	void* m_pt[64]; // Internal solver memory pointer
};
