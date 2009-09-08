#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include "WSMPSolver.h"

//////////////////////////////////////////////////////////////
// WSMPSolver
//////////////////////////////////////////////////////////////

bool WSMPSolver::PreProcess(SparseMatrix& K)
{
	// Make sure the solver is available
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform\n\n");
	return false;
#else
	// Auxiliary variables
	int idum, nrhs=1, naux=0;
	double ddum;

	CompactSymmMatrix* A = dynamic_cast<CompactSymmMatrix*> (&K);
	m_n = A->Size();
	m_nnz = A->NonZeroes();

	// Initialize m_perm and m_invp
	m_perm.create(m_n); m_perm.zero();
	m_invp.create(m_n); m_invp.zero();
	m_b.create(m_n);    m_b.zero();


	// Number of processors OMP_NUM_THREADS
	char* var = getenv("OMP_NUM_THREADS");
	int num_procs;
	if(var) num_procs = -atoi(var); // edited 6/1/09 (added negative) per Anshul Gupta
	else {
		fprintf(stderr, "Set environment OMP_NUM_THREADS to 1");
		exit(1);
	}
	wsetmaxthrds_(&num_procs);

// ------------------------------------------------------------------------------
// This step initializes 'm_iparm'
// ------------------------------------------------------------------------------

//	wsmp_initialize_();
	m_iparm[0] = 0;
	m_iparm[1] = 0;
	m_iparm[2] = 0;

	wssmp_(&m_n, A->pointers(), A->indices(), A->values(), &ddum, m_perm, m_invp,
		 m_b, &m_n, &nrhs, &ddum, &naux, &idum, m_iparm, m_dparm);

	if (m_iparm[63])
	{
		fprintf(stderr, "\nERROR during initialization: %i", m_iparm[63]);
		exit(2);
	}

	return LinearSolver::PreProcess(K);
#endif
}

bool WSMPSolver::Factor(SparseMatrix& K)
{
	// Make sure the solver is available
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform\n\n");
	return false;
#else
	// Auxiliary variables
	int idum, nrhs=1, naux=0;
	double ddum;

	CompactSymmMatrix* A = dynamic_cast<CompactSymmMatrix*> (&K);


#ifdef PRINTHB
	A->print_hb(); // Write Harwell-Boeing matrix to file
#endif

// ------------------------------------------------------------------------------
// This step performs matrix ordering
// ------------------------------------------------------------------------------

	m_iparm[1] = 1;
	m_iparm[2] = 1;
	m_dparm[9] = 1.0e-18; // matrix singularity threshold

	wssmp_(&m_n, A->pointers(), A->indices(), A->values(), &ddum, m_perm, m_invp,
		 m_b, &m_n, &nrhs, &ddum, &naux, &idum, m_iparm, m_dparm);

	if (m_iparm[63])
	{
		fprintf(stderr, "\nERROR during ordering: %i", m_iparm[63]);
		exit(2);
	}

// ------------------------------------------------------------------------------
// This step performs symbolic factorization
// ------------------------------------------------------------------------------

	m_iparm[1] = 2;
	m_iparm[2] = 2;

	wssmp_(&m_n, A->pointers(), A->indices(), A->values(), &ddum, m_perm, m_invp,
		 m_b, &m_n, &nrhs, &ddum, &naux, &idum, m_iparm, m_dparm);

	if (m_iparm[63])
	{
		fprintf(stderr, "\nERROR during ordering: %i", m_iparm[63]);
		exit(2);
	}
// ------------------------------------------------------------------------------
// This step performs Cholesky or LDLT factorization
// ------------------------------------------------------------------------------

	m_iparm[1] = 3;
	m_iparm[2] = 3;
	m_iparm[30] = 1; // 0: Cholesky factorization

	wssmp_(&m_n, A->pointers(), A->indices(), A->values(), &ddum, m_perm, m_invp,
		 m_b, &m_n, &nrhs, &ddum, &naux, &idum, m_iparm, m_dparm);

	if (m_iparm[63])
	{
		fprintf(stderr, "\nERROR during Cholesky factorization: %i", m_iparm[63]);

		if (m_iparm[63] > 0) // Try LDL factorization
		{
			m_iparm[1] = 3;
			m_iparm[2] = 3;
			m_iparm[30] = 1;

			wssmp_(&m_n, A->pointers(), A->indices(), A->values(), &ddum, m_perm, m_invp,
				 &ddum, &m_n, &nrhs, &ddum, &naux, &idum, m_iparm, m_dparm);

			if (m_iparm[63])
			{
				fprintf(stderr, "\nERROR during LDL factorization: %i", m_iparm[63]);
				exit(2);
			}
		}
	}


	return true;
#endif
}

bool WSMPSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& b)
{
	/* Make sure the solver is available */
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform\n\n");
	return false;
#else

	/* Auxiliary variables */
	int i, idum, nrhs=1, naux=0;
	double ddum;

	CompactSymmMatrix* A = dynamic_cast<CompactSymmMatrix*> (&K);

// ------------------------------------------------------------------------------
// This step performs back substitution
// ------------------------------------------------------------------------------

	m_iparm[1] = 4;
	m_iparm[2] = 4;

	wssmp_(&m_n, A->pointers(), A->indices(), A->values(), &ddum, m_perm, m_invp,
		 b, &m_n, &nrhs, &ddum, &naux, &idum, m_iparm, m_dparm);

	if (m_iparm[63])
	{
		fprintf(stderr, "\nERROR during ordering: %i", m_iparm[63]);
		exit(2);
	}

	for (i=0; i<m_n; i++) x[i] = b[i];

	return true;
#endif
}

bool WSMPSolver::Solve(SparseMatrix& K, matrix& x, matrix& b)
{
	/* Make sure the solver is available */
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform\n\n");
	return false;
#else

	//TODO: implement this solver routine for this class

	return false;
#endif
}

void WSMPSolver::Destroy()
{
	/* Make sure the solver is available */
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform\n\n");
	exit(1);
#else

	wsmp_clear_();
	LinearSolver::Destroy();

#endif
}
