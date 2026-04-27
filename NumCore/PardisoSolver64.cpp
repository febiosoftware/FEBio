/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include "PardisoSolver64.h"
#include "MatrixTools.h"
#include <FECore/CompactSymmMatrix64.h>
#include <FECore/log.h>

#ifdef PARDISO
#undef PARDISO
#include <mkl.h>
#include <mkl_pardiso.h>
#define PARDISO
#endif

class PardisoSolver64::Imp
{
public:

	CompactSymmMatrix64*	pA = nullptr;
	long long		mtype = -2; // matrix type

	// Pardiso control parameters
	long long iparm[64] = { 0 };
	long long maxfct = 0, mnum = 0;
	int msglvl = 0;
	double dparm[64] = { 0 };

	bool iparm3 = false;	// use direct-iterative method

	// Matrix data
	long long n = 0, nnz = 0, nrhs = 0;

	bool	isFactored = false;

	void* pt[64] = { nullptr }; // Internal solver memory pointer

	void print_err(long long nerror);
};

//-----------------------------------------------------------------------------
// print pardiso error message
void PardisoSolver64::Imp::print_err(long long nerror)
{
	switch (-nerror)
	{
	case 1: fprintf(stderr, "Inconsistent input\n"); break;
	case 2: fprintf(stderr, "Not enough memory\n"); break;
	case 3: fprintf(stderr, "Reordering problem\n"); break;
	case 4: fprintf(stderr, "Zero pivot, numerical fact. or iterative refinement problem\n"); break;
	case 5: fprintf(stderr, "Unclassified (internal) error\n"); break;
	case 6: fprintf(stderr, "Preordering failed\n"); break;
	case 7: fprintf(stderr, "Diagonal matrix problem\n"); break;
	case 8: fprintf(stderr, "32-bit integer overflow problem\n"); break;
	default:
		fprintf(stderr, " Unknown\n");
	}
}

BEGIN_FECORE_CLASS(PardisoSolver64, LinearSolver)
	ADD_PARAMETER(m.iparm3  , "precondition");
	ADD_PARAMETER(m.msglvl  , "msglvl");
END_FECORE_CLASS();

PardisoSolver64::PardisoSolver64(FEModel* fem) : LinearSolver(fem), m(*(new Imp))
{

}

PardisoSolver64::~PardisoSolver64()
{
#ifdef PARDISO
	MKL_Free_Buffers();
#endif
}

void PardisoSolver64::UseIterativeFactorization(bool b)
{
	m.iparm3 = b;
}

SparseMatrix* PardisoSolver64::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef PARDISO
	// allocate the correct matrix format depending on matrix symmetry type
	switch (ntype)
	{
	case REAL_SYMMETRIC     : m.mtype = -2; m.pA = new CompactSymmMatrix64(1); break;
	default:
		assert(false);
		m.pA = nullptr;
	}

	return m.pA;
#else
	return nullptr;
#endif
}

bool PardisoSolver64::SetSparseMatrix(SparseMatrix* pA)
{
#ifdef PARDISO
	if (m.pA && m.isFactored) Destroy();
	m.pA = dynamic_cast<CompactSymmMatrix64*>(pA);
	m.mtype = -2;
	return (m.pA != nullptr);
#else
	return false;
#endif
}

bool PardisoSolver64::PreProcess()
{
#ifdef PARDISO
	assert(m.isFactored == false);
	m.iparm[0] = 1; // supply all values
	m.iparm[1] = 2; // nested dissection algorithm from METIS
	m.iparm[9] = 8; // The default value for symmetric indefinite matrices (mtype =-2, mtype=-4, mtype=6), eps = 1e-8
	m.iparm[17] = 0; // Report the number of non-zero elements in the factors.(disable)
	m.iparm[18] = 0; // Disable report.
	m.iparm[20] = 1; // Pivoting for symmetric indefinite matrices. ( for matrices of mtype=-2, mtype=-4, or mtype=6.)

	m.n = m.pA->Rows();
	m.nnz = m.pA->NonZeroes();
	m.nrhs = 1;

	m.maxfct = 1;	/* Maximum number of numerical factorizations */
	m.mnum = 1;	/* Which factorization to use */

	return LinearSolver::PreProcess();
#else
	return false;
#endif
}

bool PardisoSolver64::Factor()
{
#ifdef PARDISO
	// make sure we have work to do
	if (m.pA->Rows() == 0) return true;

// ------------------------------------------------------------------------------
// Reordering and Symbolic Factorization.  This step also allocates all memory
// that is necessary for the factorization.
// ------------------------------------------------------------------------------

	long long phase = 11;
	long long error = 0;
	long long msglvl = m.msglvl;
	pardiso_64(m.pt, &m.maxfct, &m.mnum, &m.mtype, &phase, &m.n, m.pA->Values(), m.pA->Pointers64(), m.pA->Indices64(),
		 NULL, &m.nrhs, m.iparm, &msglvl, NULL, NULL, &error);

	if (error)
	{
		fprintf(stderr, "\nERROR during symbolic factorization: ");
		m.print_err(error);
		exit(2);
	}

	if (m.msglvl == 1)
	{
		long long* ip = m.iparm;
		fprintf(stdout, "\nMemory info:\n");
		fprintf(stdout, "============\n");
		fprintf(stdout, "Peak memory on symbolic factorization ............. : %d KB\n", (int)ip[14]);
		fprintf(stdout, "Permanent memory on symbolic factorization ........ : %d KB\n", (int)ip[15]);
		fprintf(stdout, "Peak memory on numerical factorization and solution : %d KB\n", (int)ip[16]);
		fprintf(stdout, "Total peak memory ................................. : %d KB\n\n", (int)max(ip[14], ip[15]+ip[16]));
	}

// ------------------------------------------------------------------------------
// This step does the factorization
// ------------------------------------------------------------------------------

	m.iparm[3] = (m.iparm3 ? 61 : 0);
	phase = 22;
	error = 0;
	pardiso_64(m.pt, &m.maxfct, &m.mnum, &m.mtype, &phase, &m.n, m.pA->Values(), m.pA->Pointers64(), m.pA->Indices64(),
		 NULL, &m.nrhs, m.iparm, &msglvl, NULL, NULL, &error);

	if (error)
	{
		fprintf(stderr, "\nERROR during factorization: ");
		m.print_err(error);
		return false;
	}

	m.isFactored = true;

	return true;
#else
	return false;
#endif
}

bool PardisoSolver64::BackSolve(double* x, double* b)
{
#ifdef PARDISO
	// make sure we have work to do
	if (m.pA->Rows() == 0) return true;

	long long phase = 33;

	m.iparm[7] = 1;	/* Maximum number of iterative refinement steps */

	long long error = 0;
	long long msglvl = m.msglvl;
	pardiso_64(m.pt, &m.maxfct, &m.mnum, &m.mtype, &phase, &m.n, m.pA->Values(), m.pA->Pointers64(), m.pA->Indices64(),
		 NULL, &m.nrhs, m.iparm, &msglvl, b, x, &error);

	if (error)
	{
		fprintf(stderr, "\nERROR during solution: ");
		m.print_err(error);
		exit(3);
	}

	// update stats
	UpdateStats(1);

	return true;
#else
	return false;
#endif
}

void PardisoSolver64::Destroy()
{
#ifdef PARDISO
	if (m.pA && m.isFactored)
	{
		long long phase = -1;
		long long error = 0;
		long long msglvl = m.msglvl;
		pardiso_64(m.pt, &m.maxfct, &m.mnum, &m.mtype, &phase, &m.n, NULL, m.pA->Pointers64(), m.pA->Indices64(),
			NULL, &m.nrhs, m.iparm, &msglvl, NULL, NULL, &error);
	}
	m.isFactored = false;
#endif
}
