#ifndef PARDISO

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include "PardisoSolver.h"

/*
#if defined(WIN32) && defined(PARDISO_DLL)
typedef int (*PARDISOINITFNC)(void*, int*, int*);
typedef int (*PARDISOFNC)(void *, int *, int *, int *, int *, int *,
		double *, int *, int *, int *, int *, int *,
		int *, double *, double *, int *);

PARDISOINITFNC pardisoinit_;
PARDISOFNC pardiso_;

#ifndef PARDISO
#define PARDISO
#endif

#include "windows.h"
HINSTANCE HPARDISO_DLL = 0;
#elif defined PARDISO
extern "C"
{
	int pardisoinit_(void *, int *, int *);

	int pardiso_(void *, int *, int *, int *, int *, int *,
		double *, int *, int *, int *, int *, int *,
		int *, double *, double *, int *);
}
#endif // PARDISO_DLL
*/

//////////////////////////////////////////////////////////////
// PardisoSolver
//////////////////////////////////////////////////////////////

PardisoSolver::PardisoSolver()
{
	/* Make sure the solver is available */
#ifndef PARDISODL
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	exit(1);
/*
#elif defined(WIN32) && defined(PARDISO_DLL)
	HPARDISO_DLL = LoadLibraryA("libpardiso.dll");
	if (HPARDISO_DLL)
		fprintf(stderr, "Pardiso library loaded successfully.\n");
	else
	{
		fprintf(stderr, "Failed loading pardiso library.\n");
		exit(1);
	}

	pardisoinit_ = (PARDISOINITFNC) GetProcAddress(HPARDISO_DLL, "pardisoinit_");
	if (pardisoinit_ == 0) exit(1);
	pardiso_ = (PARDISOFNC) GetProcAddress(HPARDISO_DLL, "pardiso_");
	if (pardiso_ == 0) exit(1);
*/
#endif
}

bool PardisoSolver::PreProcess()
{
	/* Make sure the solver is available */
#ifndef PARDISODL
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	return false;
#else
	m_mtype = (m_bsymm ? -2 : 11); /* Real symmetric matrix */
	m_iparm[0] = 0;
	int solver = 0; /* 0 - sparce direct solver, 1 - multi-recursive iterative solver */
	pardisoinit_(m_pt, &m_mtype, &solver, m_iparm, m_dparm, &m_error);

	if (m_error)
	{
		if (m_error == -10) fprintf(stderr, "No license file found\n\n");
		if (m_error == -11) fprintf(stderr, "License is expired\n\n");
		if (m_error == -12) fprintf(stderr, "Wrong username or hostname\n\n");
		return false;
	}
	else fprintf(stderr, "PARDISO license check was successful\n\n");

	m_n = m_pA->Size();
	m_nnz = m_pA->NonZeroes();
	m_nrhs = 1;

	// number of processors: use value of OMP_NUM_THREADS
	m_iparm[2] = m_numthreads;

	m_maxfct = 1;	/* Maximum number of numerical factorizations */
	m_mnum = 1;	/* Which factorization to use */

	m_msglvl = 0;	/* 0 Suppress printing, 1 Print statistical information */

	return LinearSolver::PreProcess();
#endif
}

bool PardisoSolver::Factor()
{
	/* Make sure the solver is available */
#ifndef PARDISODL
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	return false;
#else
	CompactMatrix* A = dynamic_cast<CompactMatrix*> (m_pA);

// ------------------------------------------------------------------------------
// Reordering and Symbolic Factorization.  This step also allocates all memory
// that is necessary for the factorization.
// ------------------------------------------------------------------------------

	int phase = 11;

	pardiso_(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, A->values(), A->pointers(), A->indices(),
		 NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &m_error, m_dparm);

	if (m_error)
	{
		fprintf(stderr, "\nERROR during symbolic factorization: ");
		print_err();
		exit(2);
	}

// ------------------------------------------------------------------------------
// This step does the factorization
// ------------------------------------------------------------------------------

	phase = 22;

#ifdef PRINTHB
	A->print_hb();
#endif

	pardiso_(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, A->values(), A->pointers(), A->indices(),
		 NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &m_error, m_dparm);

	if (m_error)
	{
		fprintf(stderr, "\nERROR during factorization: ");
		print_err();
		exit(2);
	}

	return true;
#endif
}

bool PardisoSolver::Solve(vector<double>& x, vector<double>& b)
{
	/* Make sure the solver is available */
#ifndef PARDISODL
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	return false;
#else
	CompactMatrix* A = dynamic_cast<CompactMatrix*> (m_pA);

	int phase = 33;

	m_iparm[7] = 1;	/* Maximum number of iterative refinement steps */

	pardiso_(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, A->values(), A->pointers(), A->indices(),
		 NULL, &m_nrhs, m_iparm, &m_msglvl, b, x, &m_error, m_dparm);

	if (m_error)
	{
		fprintf(stderr, "\nERROR during solution: ");
		print_err();
		exit(3);
	}

	return true;
#endif
}

void PardisoSolver::Destroy()
{
	/* Make sure the solver is available */
#ifndef PARDISODL
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	exit(1);
#else
	int phase = -1;

	pardiso_(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, NULL, NULL, NULL,
		 NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &m_error, m_dparm);

	LinearSolver::Destroy();

#endif
}
#endif
