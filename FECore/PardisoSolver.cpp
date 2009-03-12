#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include "PardisoSolver.h"

//////////////////////////////////////////////////////////////
// PardisoSolver
//////////////////////////////////////////////////////////////

PardisoSolver::PardisoSolver()
{
	/* Make sure the solver is available */
#ifndef PARDISO
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	exit(1);
#else
	m_mtype = -2; /* Real symmetric matrix */
	m_iparm[0] = 0;
	pardisoinit_(m_pt, &m_mtype, m_iparm);
#endif
}

bool PardisoSolver::PreProcess(SparseMatrix& K)
{
	/* Make sure the solver is available */
#ifndef PARDISO
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	return false;
#else
	/* Auxiliary variables */
	double ddum; /* Double dummy */
	int idum; /* Integer dummy */

	CompactMatrix* A = dynamic_cast<CompactMatrix*> (&K);
	m_n = A->Size();
	m_nnz = A->NonZeroes();
	m_nrhs = 1;

	/* Number of processors OMP_NUM_THREADS */
	char* var = getenv("OMP_NUM_THREADS");
	int num_procs;
	if(var) num_procs = atoi(var);
	else {
		fprintf(stderr, "Set environment OMP_NUM_THREADS to 1");
		exit(1);
	}
	m_iparm[2] = num_procs;

	m_maxfct = 1;	/* Maximum number of numerical factorizations */
	m_mnum = 1;	/* Which factorization to use */

	m_msglvl = 0;	/* 0 Suppress printing, 1 Print statistical information */
	m_error = 0;	/* Initialize m_error flag */

// ------------------------------------------------------------------------------
// Reordering and Symbolic Factorization.  This step also allocates all memory
// that is necessary for the factorization.
// ------------------------------------------------------------------------------

	int phase = 11;

	pardiso_(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, A->values(), A->pointers(), A->indices(),
		 &idum, &m_nrhs, m_iparm, &m_msglvl, &ddum, &ddum, &m_error);

	if (m_error)
	{
		fprintf(stderr, "\nERROR during symbolic factorization: %i", m_error);
		exit(2);
	}

//	fprintf(stderr, "\nReordering completed ...\n");
//	fprintf(stderr, "\nNumber of nonzeros in factors = %i", m_iparm[17]);
//	fprintf(stderr, "\nNumber of factorization MFLOPS = %i", m_iparm[18]);

	return LinearSolver::PreProcess(K);
#endif
}

bool PardisoSolver::Factor(SparseMatrix& K)
{
	/* Make sure the solver is available */
#ifndef PARDISO
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	return false;
#else
	/* Auxiliary variables */
	double ddum; /* Double dummy */
	int idum; /* Integer dummy */

	CompactMatrix* A = dynamic_cast<CompactMatrix*> (&K);

	int phase = 22;

#ifdef DEBUG

	int i, *pointers, *indices;
	double* values;

	pointers = A->pointers();
	indices = A->indices();
	values = A->values();
	fprintf(stdout, "\nPointers:");
	for (i=0; i<=m_n; i++) fprintf(stdout, "\n%d", pointers[i]);
	fprintf(stdout, "\nIndices, Values:");
	for (i=0; i<m_nnz; i++) fprintf(stdout, "\n%d, %g", indices[i], values[i]);
#endif

	pardiso_(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, A->values(), A->pointers(), A->indices(),
		 &idum, &m_nrhs, m_iparm, &m_msglvl, &ddum, &ddum, &m_error);

	if (m_error)
	{
		fprintf(stderr, "\nERROR during factorization: %i", m_error);
		exit(2);
	}

//	cout << "\nFactorization completed ...\n";

	return true;
#endif
}

bool PardisoSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& b)
{
	/* Make sure the solver is available */
#ifndef PARDISO
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	return false;
#else

	/* Auxiliary variables */
	int idum; /* Integer dummy */

	CompactMatrix* A = dynamic_cast<CompactMatrix*> (&K);

	int phase = 33;

	m_iparm[7] = 1;	/* Maximum number of iterative refinement steps */

	pardiso_(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, A->values(), A->pointers(), A->indices(),
		 &idum, &m_nrhs, m_iparm, &m_msglvl, b, x, &m_error);

	if (m_error)
	{
		fprintf(stderr, "\nERROR during solution: %i", m_error);
		exit(3);
	}

//	cout << "\nSolve completed ...\n";

	return true;
#endif
}

bool PardisoSolver::Solve(SparseMatrix& K, matrix& x, matrix& b)
{
	/* Make sure the solver is available */
#ifndef PARDISO
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	return false;
#else

	//TODO: implement this solver routine for this class

	return false;
#endif
}

void PardisoSolver::Destroy()
{
	/* Make sure the solver is available */
#ifndef PARDISO
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	exit(1);
#else
	/* Auxiliary variables */
	double ddum; /* Double dummy */
	int idum; /* Integer dummy */

	int phase = -1;

	pardiso_(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, &ddum, &idum, &idum,
		 &idum, &m_nrhs, m_iparm, &m_msglvl, &ddum, &ddum, &m_error);

	LinearSolver::Destroy();

#endif
}
