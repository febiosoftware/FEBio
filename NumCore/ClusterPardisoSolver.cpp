/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "ClusterPardisoSolver.h"
#include <FECore/CompactMatrix.h>
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/CompactSymmMatrix.h>

#ifdef PARDISO
#undef PARDISO
#endif

#include <mpi.h>
#include <mkl.h>
#include <mkl_cluster_sparse_solver.h>


struct ClusterPardisoSolver::Imp
{
	int comm; // MPI communicator
	MKL_INT mtype = 0; // matrix type for cpardiso
	CompactMatrix* A = nullptr;
	MKL_INT iparm[64] = { 0 };
	void* pt[64] = { 0 };
	bool isFactored = false;
};

ClusterPardisoSolver::ClusterPardisoSolver(FEModel* fem) : LinearSolver(fem), m(*new Imp)
{
	m.comm = MPI_Comm_c2f(MPI_COMM_WORLD);
}

ClusterPardisoSolver::~ClusterPardisoSolver()
{
	delete& m;
}

SparseMatrix* ClusterPardisoSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	// allocate the correct matrix format depending on matrix symmetry type
	switch (ntype)
	{
	case REAL_SYMMETRIC: m.mtype = -2; m.A = new CompactSymmMatrix(1); break;
	case REAL_UNSYMMETRIC: m.mtype = 11; m.A = new CRSSparseMatrix(1); break;
	case REAL_SYMM_STRUCTURE: m.mtype = 1; m.A = new CRSSparseMatrix(1); break;
	default:
		assert(false);
		m.A = nullptr;
	}

	return m.A;
}

bool ClusterPardisoSolver::SetSparseMatrix(SparseMatrix* pA)
{
	return false;
}

bool ClusterPardisoSolver::PreProcess()
{
	m.iparm[0] = 1; /* Solver default parameters overridden with provided by iparm */
	m.iparm[1] = 2; /* Use METIS for fill-in reordering */
	m.iparm[5] = 0; /* Write solution into x */
	m.iparm[7] = 2; /* Max number of iterative refinement steps */
	m.iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	m.iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	m.iparm[12] = 1; /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
	m.iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	m.iparm[18] = -1; /* Output: Mflops for LU factorization */
	m.iparm[26] = 1; /* Check input data for correctness */
	m.iparm[39] = 0; /* Input: matrix/rhs/solution stored on master */

	return true;
}

bool ClusterPardisoSolver::Factor()
{
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates   */
	/* all memory that is necessary for the factorization.                  */
	/* -------------------------------------------------------------------- */
	MKL_INT maxfct, mnum, msglvl, error;
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Print statistical information in file */
	error = 0; /* Initialize error flag */

	MKL_INT n = m.A->Rows(); // Number of equations
	MKL_INT nrhs = 1; // Number of right hand sides.
	double* a = m.A->Values();
	MKL_INT* ia = m.A->Pointers();
	MKL_INT* ja = m.A->Indices();
	MKL_INT idum; /* Dummy integer */
	double ddum; /* Dummy double */

	// tell all workers do the symbolic factorization and memory allocation
	int cmd = ClusterPardisoRuntime::SYMBOLIC;
	MPI_Bcast(&cmd, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MKL_INT phase = 11;
	cluster_sparse_solver(m.pt, &maxfct, &mnum, &m.mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, m.iparm, &msglvl, &ddum, &ddum, &m.comm, &error);
	if (error != 0)
	{
		fprintf(stderr, "\nERROR during symbolic factorization: %lli", (long long int)error);
		return false;
	}
	fprintf(stderr, "\nReordering completed ... ");

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization.                                          */
	/* -------------------------------------------------------------------- */

	// tell all workers do the numerical factorization
	cmd = ClusterPardisoRuntime::FACTOR;
	MPI_Bcast(&cmd, 1, MPI_INT, 0, MPI_COMM_WORLD);

	phase = 22;
	cluster_sparse_solver(m.pt, &maxfct, &mnum, &m.mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, m.iparm, &msglvl, &ddum, &ddum, &m.comm, &error);
	if (error != 0)
	{
		fprintf(stderr, "\nERROR during numerical factorization: %lli", (long long int)error);
		return false;
	}
	fprintf(stderr, "\nFactorization completed ... ");
	m.isFactored = true;

	return true;
}

bool ClusterPardisoSolver::BackSolve(double* x, double* y)
{
	MKL_INT n = m.A->Rows(); // Number of equations
	MKL_INT nrhs = 1; // Number of right hand sides.
	double* a = m.A->Values();
	MKL_INT* ia = m.A->Pointers();
	MKL_INT* ja = m.A->Indices();

	MKL_INT maxfct, mnum, msglvl, error;
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Print statistical information in file */
	error = 0; /* Initialize error flag */

	MKL_INT idum; /* Dummy integer */
	double ddum; /* Dummy double */

	// let all workers do the backsolve
	int cmd = ClusterPardisoRuntime::BACKSOLVE;
	MPI_Bcast(&cmd, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement.                       */
	/* -------------------------------------------------------------------- */
	MKL_INT phase = 33;
	fprintf(stderr, "\nSolving system...");
	cluster_sparse_solver(m.pt, &maxfct, &mnum, &m.mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, m.iparm, &msglvl, y, x, &m.comm, &error);
	if (error != 0)
	{
		fprintf(stderr, "\nERROR during solution: %lli", (long long int)error);
		return false;
	}

	return true;
}

void ClusterPardisoSolver::Destroy()
{
	if ((m.A == nullptr) || !m.isFactored) return;

	MKL_INT maxfct, mnum, msglvl, error;
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Print statistical information in file */
	error = 0; /* Initialize error flag */

	MKL_INT idum; /* Dummy integer */
	double ddum; /* Dummy double */

	MKL_INT n = m.A->Rows(); // Number of equations
	MKL_INT nrhs = 1; // Number of right hand sides.
	double* a = m.A->Values();
	MKL_INT* ia = m.A->Pointers();
	MKL_INT* ja = m.A->Indices();

	// let all workers do the destroy
	int cmd = ClusterPardisoRuntime::DESTROY;
	MPI_Bcast(&cmd, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MKL_INT phase = -1; /* Release internal memory. */
	cluster_sparse_solver(m.pt, &maxfct, &mnum, &m.mtype, &phase,
		&n, &ddum, ia, ja, &idum, &nrhs, m.iparm, &msglvl, &ddum, &ddum, &m.comm, &error);
	if (error != 0)
	{
		fprintf(stderr, "\nERROR during release memory: %lli", (long long int)error);
	}
	m.isFactored = false;
}

void ClusterPardisoRuntime::WorkerLoop()
{
	MKL_INT iparm[64] = { 0 };
	void* pt[64] = { 0 };
	int comm = MPI_Comm_c2f(MPI_COMM_WORLD);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	MKL_INT maxfct, mnum, msglvl, error;
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Print statistical information in file */
	error = 0; /* Initialize error flag */
	int n = 0;

	MKL_INT idum; /* Dummy integer */
	double ddum; /* Dummy double */

	while (true)
	{
		fprintf(stderr, "[%d] waiting for command\n", rank);
		fflush(stderr);

		int cmd = 0;
		MPI_Bcast(&cmd, 1, MPI_INT, 0, MPI_COMM_WORLD);

		fprintf(stderr, "[%d] received command %d\n", rank, cmd);
		fflush(stderr);

		MKL_INT phase;
		switch (cmd)
		{
		case SYMBOLIC : phase = 11; break;
		case FACTOR   : phase = 22; break;
		case BACKSOLVE: phase = 33; break;
		case DESTROY  : phase = -1; break;
		case FINISH:
			return;
		}

		fprintf(stderr, "worker rank %d calling cpardiso with phase %d.\n", rank, phase);
		fflush(stderr);

		cluster_sparse_solver(pt, &maxfct, &mnum, &idum, &phase,
			&n, nullptr, nullptr, nullptr, &idum, &idum, iparm, &msglvl, &ddum, &ddum, &comm, &error);

		if (error != 0)
		{
			fprintf(stderr, "worker rank %d: cpardiso error %d\n", rank, error);
			fflush(stderr);
			MPI_Abort(MPI_COMM_WORLD, error);
		}
	}
}

void ClusterPardisoRuntime::Finish()
{
	int cmd = ClusterPardisoRuntime::FINISH;
	MPI_Bcast(&cmd, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
