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
#include "SuperLU_MT.h"
#include "MatrixTools.h"
#include <FECore/log.h>
#ifdef HAVE_SUPERLU_MT
#define __OPENMP
#undef MAX
#undef MIN
#include <slu_mt_ddefs.h>
#endif // HAVE_SUPERLU_MT

#ifdef HAVE_SUPERLU_MT

//-----------------------------------------------------------------------------
class SuperLU_MT_Solver::Impl
{
public:
	SuperMatrix A, L, U;
	SuperMatrix B, X;
	SuperMatrix AC;
	vector<int> perm_c;
	vector<int> perm_r;
	bool isFactored = false;
	superlumt_options_t options;
	Gstat_t  Gstat;
	int	nprocs = 1;
	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering
	 *   permc_spec = 1: minimum degree ordering on structure of A'*A
	 *   permc_spec = 2: minimum degree ordering on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */
	int permc_spec = 3;
};

BEGIN_FECORE_CLASS(SuperLU_MT_Solver, LinearSolver)
	ADD_PARAMETER(m->nprocs, "nprocs");
	ADD_PARAMETER(m->permc_spec, "permc_spec");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
SuperLU_MT_Solver::SuperLU_MT_Solver(FEModel* fem) : LinearSolver(fem), m_pA(nullptr)
{
	m = new SuperLU_MT_Solver::Impl;
}

//-----------------------------------------------------------------------------
SuperLU_MT_Solver::~SuperLU_MT_Solver()
{
	Destroy();
}

//-----------------------------------------------------------------------------
SparseMatrix* SuperLU_MT_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	// allocate the correct matrix format depending on matrix symmetry type
	switch (ntype)
	{
	case REAL_SYMMETRIC: 
	case REAL_UNSYMMETRIC: 
	case REAL_SYMM_STRUCTURE: 
		m_pA = new CCSSparseMatrix;
		break;
	default:
		assert(false);
		m_pA = nullptr;
	}

	return m_pA;
}

//-----------------------------------------------------------------------------
bool SuperLU_MT_Solver::SetSparseMatrix(SparseMatrix* pA)
{
	return (m_pA != nullptr);
}

//-----------------------------------------------------------------------------
bool SuperLU_MT_Solver::PreProcess()
{
	if (m_pA == nullptr) return false;

	int n = m_pA->Rows();
	int nnz = m_pA->NonZeroes();
	dCreate_CompCol_Matrix(&m->A, n, n, nnz, m_pA->Values(), m_pA->Indices(), m_pA->Pointers(), SLU_NC, SLU_D, SLU_GE);

	m->perm_c.resize(n);
	m->perm_r.resize(n);

	return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool SuperLU_MT_Solver::Factor()
{
	// make sure we have work to do
	if (m_pA->Rows() == 0) return true;

	if (m->isFactored)
	{
		Destroy_SuperNode_SCP(&m->L);
		Destroy_CompCol_NCP(&m->U);
		m->isFactored = false;
	}

	int n = m_pA->Rows();
	int nnz = m_pA->NonZeroes();

	int panel_size = sp_ienv(1);
	int relax = sp_ienv(2);
	StatAlloc(n, m->nprocs, panel_size, relax, &m->Gstat);
	StatInit(n, m->nprocs, &m->Gstat);

	get_perm_c(m->permc_spec, &m->A, m->perm_c.data());

	/* ------------------------------------------------------------
	   Initialize the option structure pdgstrf_options using the
	   user-input parameters;
	   Apply perm_c to the columns of original A to form AC.
	   ------------------------------------------------------------*/
	fact_t fact = fact_t::DOFACT;
	trans_t trans = trans_t::NOTRANS;
	yes_no_t refact = yes_no_t::NO;
	double diag_pivot_thresh = 1.0;
	yes_no_t usepr = NO;
	double drop_tol = 0.0;
	pdgstrf_init(m->nprocs, fact, trans, refact, panel_size, relax,
		diag_pivot_thresh, usepr, drop_tol, m->perm_c.data(), m->perm_r.data(),
		NULL, 0, &m->A, &m->AC, &m->options, &m->Gstat);

	/* ------------------------------------------------------------
	   Compute the LU factorization of A.
	   The following routine will create nprocs threads.
	   ------------------------------------------------------------*/
	int info = 0;
	pdgstrf(&m->options, &m->AC, m->perm_r.data(), &m->L, &m->U, &m->Gstat, &info);

	m->isFactored = true;

	return (info == 0);
}

//-----------------------------------------------------------------------------
bool SuperLU_MT_Solver::BackSolve(double* x, double* b)
{
	// make sure we have work to do
	if (m_pA->Rows() == 0) return true;

	int n = m_pA->Rows();
	memcpy(x, b, n * sizeof(double));
	dCreate_Dense_Matrix(&m->B, n, 1, x, n, SLU_DN, SLU_D, SLU_GE);

	/* ------------------------------------------------------------
		Solve the system A*X=B, overwriting B with X.
		------------------------------------------------------------*/
	int info = 0;
	dgstrs(trans_t::NOTRANS, &m->L, &m->U, m->perm_r.data(), m->perm_c.data(), &m->B, &m->Gstat, &info);

	// update stats
	UpdateStats(1);

	return true;
}

//-----------------------------------------------------------------------------
void SuperLU_MT_Solver::Destroy()
{
//	pdgstrf_finalize(&m->options, &m->AC);
//	StatFree(&m->Gstat);
	if (m->isFactored)
	{
		Destroy_SuperNode_SCP(&m->L);
		Destroy_CompCol_NCP(&m->U);
		m->isFactored = false;
	}
}
#else // HAVE_SUPERLU_MT
BEGIN_FECORE_CLASS(SuperLU_MT_Solver, LinearSolver)
END_FECORE_CLASS();

SuperLU_MT_Solver::SuperLU_MT_Solver(FEModel* fem) : LinearSolver(fem), m_pA(nullptr), m(nullptr) {}
SuperLU_MT_Solver::~SuperLU_MT_Solver() {}
bool SuperLU_MT_Solver::PreProcess() { return false; }
bool SuperLU_MT_Solver::Factor() { return false; }
bool SuperLU_MT_Solver::BackSolve(double* x, double* y) { return false; }
void SuperLU_MT_Solver::Destroy() {}
SparseMatrix* SuperLU_MT_Solver::CreateSparseMatrix(Matrix_Type ntype) { return nullptr; }
bool SuperLU_MT_Solver::SetSparseMatrix(SparseMatrix* pA) { return false; }
#endif // HAVE_SUPERLU_MT
