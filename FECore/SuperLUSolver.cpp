#include "stdafx.h"
#include <stdio.h>
#include "SuperLUSolver.h"
#include "math.h"

//////////////////////////////////////////////////////////////////////
// SUPERLUSolver
//////////////////////////////////////////////////////////////////////

bool SuperLUSolver::PreProcess()
{
	// Make sure the solver is available
#ifndef SUPERLU
	fprintf(stderr, "FATAL ERROR: The SUPERLU solver is not available on this platform\n\n");
	return false;
#else

	// get a reference to the correct matrix type
	CompactUnSymmMatrix& rK = dynamic_cast<CompactUnSymmMatrix&> (*m_pA);

	int N = rK.Size();
	int nnz = rK.NonZeroes();

	// set custom options
	set_default_options(&options);
	options.ColPerm			= MMD_AT_PLUS_A;
	options.SymmetricMode	= (m_bsymm ? YES : NO);		// we use symmetric mode
	options.DiagPivotThresh	= 0.0;
	options.Equil			= NO;		// no equilibration

	// create the supermatrix A
    dCreate_CompCol_Matrix(&A, N, N, nnz, rK.values(), rK.indices(), rK.pointers(), SLU_NC, SLU_D, SLU_GE);

	// since we don't have any values yet, we don't supply any data to the matrix B and X
	dCreate_Dense_Matrix(&B, N, 0, NULL, N, SLU_DN, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&X, N, 0, NULL, N, SLU_DN, SLU_D, SLU_GE);


	// allocate storage for the permutation matrices
	perm_c.resize(N);
	perm_r.resize(N);

	// allocate storage for elimination tree
	etree.resize(N);

	m_balloc = true;

	return LinearSolver::PreProcess();
#endif
}

#ifdef SUPERLU
double SuperLUSolver::norm(SparseMatrix& K)
{
	double n = 0, nc;
	int l;

	// get a reference to the correct matrix type
	CompactUnSymmMatrix& A = dynamic_cast<CompactUnSymmMatrix&> (K);

	int* ptr = A.pointers();
	double* pval = A.values(), *pv;

	for (int i=0; i<A.Size(); ++i)
	{
		nc = 0;
		l = ptr[i+1] - ptr[i];
		pv = pval + (ptr[i]);
		for (int j=0; j<l; ++j, ++pv) nc += fabs(*pv);

		if (nc > n) n = nc;
	}

	return n;
}
#endif

bool SuperLUSolver::Factor()
{
	// Make sure the solver is available
#ifndef SUPERLU
	fprintf(stderr, "FATAL ERROR: The SUPERLU solver is not available on this platform\n\n");
	return false;
#else

	if (m_bfact)
	{
		// deallocate storage for  L and U
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
	}

	double normA = 0;
	if (m_bcond) normA = norm(*m_pA);

	// initialize stats
	StatInit(&stat);

	// set nr of columns of B to zero to make sure we don't solve anything here
	B.ncol = 0;
	X.ncol = 0;

	// factorize the matrix
	options.Fact = DOFACT;
    dgssvx(&options, &A, &perm_c[0], &perm_r[0], &etree[0], equed, NULL, NULL,
           &L, &U, NULL, 0, &B, &X, &rpg, &rcond, &ferr, &berr,
           &mem_usage, &stat, &info);

	m_bfact = true;

	// calculate condition number
	if (m_bcond)
	{
		double rcond;
		int info;
		char cnorm = '1';
		dgscon(&cnorm, &L, &U, normA, &rcond, &stat, &info);
		if (info == 0)
		{
			fprintf(stdout, " ESTIMATED CONDITION NUMBER : %lg\n", 1./rcond);
		}
		else
		{
			fprintf(stdout, " FAILED ESTIMAING CONDITION NUMBER\n");
		}
	}

	// free the stats
	StatFree(&stat);

	return true;

#endif
}

bool SuperLUSolver::Solve(vector<double>& x, vector<double>& b)
{
	// Make sure the solver is available
#ifndef SUPERLU
	fprintf(stderr, "FATAL ERROR: The SUPERLU solver is not available on this platform\n\n");
	return false;
#else

	int info;

	// set the data in the B matrix
	DNformat *Bstore = (DNformat*) B.Store;
	Bstore->nzval = &b[0];
	B.ncol = 1;

	// set the data in the X matrix
	DNformat *Xstore = (DNformat*) X.Store;
	Xstore->nzval = &x[0];
	X.ncol = 1;

	// initialize stats
	SuperLUStat_t stat;
	StatInit(&stat);

	// solve the system
	options.Fact = FACTORED;
    dgssvx(&options, &A, &perm_c[0], &perm_r[0], &etree[0], equed, NULL, NULL,
           &L, &U, NULL, 0, &B, &X, &rpg, &rcond, &ferr, &berr,
           &mem_usage, &stat, &info);

	// free the stats
	StatFree(&stat);

	return true;
#endif
}

void SuperLUSolver::Destroy()
{
	// Make sure the solver is available
#ifndef SUPERLU
	fprintf(stderr, "FATAL ERROR: The SUPERLU solver is not available on this platform\n\n");
	return;
#else

	// since superlu by default deallocates the memory for the matrix data
	// we can't use the destroy routines for A and B. In stead we deallocate the memory for A
	// ourselve
	if (m_balloc)
	{
		Destroy_SuperMatrix_Store(&A);
		Destroy_SuperMatrix_Store(&B);
		Destroy_SuperMatrix_Store(&X);

		m_balloc = false;
	}

	if (m_bfact)
	{
		// we can however for L and U
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);

		m_bfact = false;
	}

	LinearSolver::Destroy();

#endif
}
