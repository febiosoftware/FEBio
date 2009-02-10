#include <stdafx.h>
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
	mtype = -2; /* Real symmetric matrix */
	iparm[0] = 0;
	pardisoinit_(pt, &mtype, iparm);
#endif
}

bool PardisoSolver::PreProcess(SparseMatrix& K)
{
	/* Make sure the solver is available */
#ifndef PARDISO
	fprintf(stderr, "FATAL ERROR: The Pardiso solver is not available on this platform\n\n");
	return false;
#else
	CompactMatrix* A = dynamic_cast<CompactMatrix*> (&K);
	n = A->Size();
	nnz = A->NonZeroes();
	nrhs = 1;

	/* Number of processors OMP_NUM_THREADS */
	var = getenv("OMP_NUM_THREADS");
	if(var) num_procs = atoi(var);
	else {
		fprintf(stderr, "Set environment OMP_NUM_THREADS to 1");
		exit(1);
	}
	iparm[2] = num_procs;

	maxfct = 1;	/* Maximum number of numerical factorizations */
	mnum = 1;	/* Which factorization to use */

	msglvl = 0;	/* 0 Suppress printing, 1 Print statistical information */
	error = 0;	/* Initialize error flag */

// ------------------------------------------------------------------------------
// Reordering and Symbolic Factorization.  This step also allocates all memory
// that is necessary for the factorization.
// ------------------------------------------------------------------------------

	phase = 11;

	pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, A->values(), A->pointers(), A->indices(),
		 &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

	if (error)
	{
		fprintf(stderr, "\nERROR during symbolic factorization: %i", error);
		exit(2);
	}

//	fprintf(stderr, "\nReordering completed ...\n");
//	fprintf(stderr, "\nNumber of nonzeros in factors = %i", iparm[17]);
//	fprintf(stderr, "\nNumber of factorization MFLOPS = %i", iparm[18]);

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

	CompactMatrix* A = dynamic_cast<CompactMatrix*> (&K);

	phase = 22;

#ifdef DEBUG
	printf("\nThis is a test");

	int i, *pointers, *indices;
	double* values;

	fflush(stdout);
	pointers = A->pointers();
	indices = A->indices();
	values = A->values();
	fprintf(stdout, "\nPointers:");
	for (i=0; i<n; i++) fprintf(stdout, "\n%d", pointers[i]);
	fprintf(stdout, "\nIndices, Values:");
	for (i=0; i<nnz; i++) fprintf(stdout, "\n%d, %g", indices[i], values[i]);
#endif

	pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, A->values(), A->pointers(), A->indices(),
		 &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

	if (error)
	{
		fprintf(stderr, "\nERROR during factorization: %i", error);
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

	CompactMatrix* A = dynamic_cast<CompactMatrix*> (&K);

	phase = 33;

	iparm[7] = 1;	/* Maximum number of iterative refinement steps */

	pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, A->values(), A->pointers(), A->indices(),
		 &idum, &nrhs, iparm, &msglvl, b, x, &error);

	if (error)
	{
		fprintf(stderr, "\nERROR during solution: %i", error);
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
	phase = -1;

	pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, &idum, &idum,
		 &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

#endif
}
