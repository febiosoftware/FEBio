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
#include "MKLDSSolver.h"
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/CompactSymmMatrix.h>

#ifdef PARDISO
#undef PARDISO
#include "mkl_dss.h"
#include "mkl_types.h"
#define PARDISO
#endif

BEGIN_FECORE_CLASS(MKLDSSolver, LinearSolver)
END_FECORE_CLASS();

#ifdef PARDISO
class MKLDSSolver::Imp 
{
public:
	CompactMatrix* A = nullptr;
	_MKL_DSS_HANDLE_t handle = nullptr;
	MKL_INT opt = MKL_DSS_DEFAULTS;
	MKL_INT sym = MKL_DSS_SYMMETRIC;
};

//-----------------------------------------------------------------------------
MKLDSSolver::MKLDSSolver(FEModel* fem) : LinearSolver(fem), m(new MKLDSSolver::Imp)
{
}

//-----------------------------------------------------------------------------
MKLDSSolver::~MKLDSSolver()
{
	Destroy();
	delete m;
}

//-----------------------------------------------------------------------------
SparseMatrix* MKLDSSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	// allocate the correct matrix format depending on matrix symmetry type
	switch (ntype)
	{
	case REAL_SYMMETRIC     : m->A = new CompactSymmMatrix(1); m->sym = MKL_DSS_SYMMETRIC; break;
	case REAL_UNSYMMETRIC   : m->A = new CRSSparseMatrix(1); m->sym = MKL_DSS_NON_SYMMETRIC; break;
	case REAL_SYMM_STRUCTURE: m->A = new CRSSparseMatrix(1); m->sym = MKL_DSS_SYMMETRIC_STRUCTURE; break;
	default:
		assert(false);
		m->A = nullptr;
	}

	return m->A;
}

//-----------------------------------------------------------------------------
bool MKLDSSolver::SetSparseMatrix(SparseMatrix* pA)
{
	return (m->A != nullptr);
}

//-----------------------------------------------------------------------------
bool MKLDSSolver::PreProcess()
{
	SparseMatrix& A = *m->A;
	int nRows = A.Rows();
	int nCols = A.Columns();
	int nnz = A.NonZeroes();

	// create handle
	if (m->handle == nullptr)
	{
		MKL_INT error = dss_create(m->handle, m->opt);
		if (error != MKL_DSS_SUCCESS) return false;
	}

	// define the structure
	MKL_INT error = dss_define_structure(m->handle, m->sym, A.Pointers(), nRows, nCols, A.Indices(), nnz);
	if (error != MKL_DSS_SUCCESS) return false;

	// reorder
//	m->opt = MKL_DSS_AUTO_ORDER;
//	m->opt = MKL_DSS_METIS_ORDER;
	m->opt = MKL_DSS_METIS_OPENMP_ORDER;
	error = dss_reorder(m->handle, m->opt, 0);
	if (error != MKL_DSS_SUCCESS) return false;

	return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool MKLDSSolver::Factor()
{
	// make sure we have work to do
	if ((m->A == nullptr) || (m->A->Rows() == 0)) return true;

	SparseMatrix& A = *m->A;
//	MKL_INT type = MKL_DSS_POSITIVE_DEFINITE;
	MKL_INT type = MKL_DSS_INDEFINITE;
	MKL_INT error = dss_factor_real(m->handle, type, A.Values());
	if (error != MKL_DSS_SUCCESS) return false;

	return true;
}

//-----------------------------------------------------------------------------
bool MKLDSSolver::BackSolve(double* x, double* b)
{
	// make sure we have work to do
	if ((m->A == nullptr) || (m->A->Rows() == 0)) return true;

	int nRhs = 1;
	m->opt = MKL_DSS_DEFAULTS;
	MKL_INT error = dss_solve_real(m->handle, m->opt, b, nRhs, x);
	if (error != MKL_DSS_SUCCESS) return false;

	// update stats
	UpdateStats(1);

	return true;
}

//-----------------------------------------------------------------------------
void MKLDSSolver::Destroy()
{
	if (m->handle)
	{
		m->opt = MKL_DSS_DEFAULTS;
		MKL_INT error = dss_delete(m->handle, m->opt);
		assert(error == MKL_DSS_SUCCESS);
		m->handle = nullptr;
	}
}
#else
//-----------------------------------------------------------------------------
class MKLDSSolver::Imp {};
MKLDSSolver::MKLDSSolver(FEModel* fem) : LinearSolver(fem), m(nullptr) {}
MKLDSSolver::~MKLDSSolver() {}
SparseMatrix* MKLDSSolver::CreateSparseMatrix(Matrix_Type ntype) {	return nullptr; }
bool MKLDSSolver::SetSparseMatrix(SparseMatrix* pA) { return false; }
bool MKLDSSolver::PreProcess() { return false; }
bool MKLDSSolver::Factor() { return false; }
bool MKLDSSolver::BackSolve(double* x, double* b) { return false; }
void MKLDSSolver::Destroy() {}
#endif
