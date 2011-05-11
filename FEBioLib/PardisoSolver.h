//! This class implements the Pardiso solver.

//! The Pardiso solver is included in the Intel Math Kernel Library (MKL).
//! It can also be installed as a shared object library from
//!		http://www.pardiso-project.org

#pragma once

#include "FECore/LinearSolver.h"
#include "CompactMatrix.h"
using namespace FECore;

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
	bool PreProcess();
	bool Factor();
	bool BackSolve(vector<double>& x, vector<double>& b);
	void Destroy();

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype)
	{
		m_bsymm = (ntype == SPARSE_SYMMETRIC);
		if (m_bsymm) m_pA = new CompactSymmMatrix(1);
		else m_pA = new CompactUnSymmMatrix(1, true);

		return m_pA;
	}

	PardisoSolver();

protected:

	bool m_bsymm; // use symmetric mode or not

	// Pardiso control parameters
	int m_iparm[64];
	int m_maxfct, m_mnum, m_error, m_msglvl;
	double m_dparm[64];

	// Matrix data
	int m_mtype;
	int m_n, m_nnz, m_nrhs;

	void* m_pt[64]; // Internal solver memory pointer

	void print_err()
	{
		switch (-m_error)
		{
			case 1 : fprintf(stderr, "Inconsistent input\n"); break;
			case 2 : fprintf(stderr, "Not enough memory\n"); break;
			case 3 : fprintf(stderr, "Reordering problem\n"); break;
			case 4 : fprintf(stderr, "Zero pivot, numerical fact. or iterative refinement problem\n"); break;
			case 5 : fprintf(stderr, "Unclassified (internal) error\n"); break;
			case 6 : fprintf(stderr, "Preordering failed\n"); break;
			case 7 : fprintf(stderr, "Diagonal matrix problem\n"); break;
			case 8 : fprintf(stderr, "32-bit integer overflow problem\n"); break;
		}
	}
};
