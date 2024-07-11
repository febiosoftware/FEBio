/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.

 See Copyright-FEBio.txt for details.

 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include "UMFPACKSolver.h"
#include <FECore/CompactUnSymmMatrix.h>
#ifdef HAS_SUITESPARSE
#include <umfpack.h>
#endif

class UMFPACKSolver::Imp
{
public:
	SparseMatrix* A = nullptr;
	void* symbolic = nullptr;
	void* numeric = nullptr;
	std::vector<int64_t>	ind;
	std::vector<int64_t>	ptr;
};

BEGIN_FECORE_CLASS(UMFPACKSolver, LinearSolver)
END_FECORE_CLASS();

UMFPACKSolver::UMFPACKSolver(FEModel* fem) : LinearSolver(fem), m(new Imp)
{

}

UMFPACKSolver::~UMFPACKSolver()
{
	delete m;
}

SparseMatrix* UMFPACKSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	m->A = new CCSSparseMatrix(0);
	return m->A;
}

bool UMFPACKSolver::SetSparseMatrix(SparseMatrix* pA)
{
	return false;
}

bool UMFPACKSolver::PreProcess()
{
#ifdef HAS_SUITESPARSE
	SparseMatrix* A = m->A;
	if (A == nullptr) return false;

	int nr = A->Rows();
	int nc = A->Columns();

	int* pA = A->Pointers();
	int* iA = A->Indices();

	int nnz = A->NonZeroes();
	m->ptr.resize(nr + 1);
	for (size_t n = 0; n <= nr; ++n) m->ptr[n] = pA[n];
	m->ind.resize(nnz);
	for (size_t n = 0; n < nnz;  ++n) m->ind[n] = iA[n];

	int status = umfpack_dl_symbolic(nr, nc, m->ptr.data(), m->ind.data(), A->Values(), &m->symbolic, nullptr, nullptr);
	return (status == UMFPACK_OK);
#else
	return false;
#endif
}

bool UMFPACKSolver::Factor()
{
#ifdef HAS_SUITESPARSE
	SparseMatrix* A = m->A;
	if (A == nullptr) return false;

	int status = umfpack_dl_numeric(m->ptr.data(), m->ind.data(), A->Values(), m->symbolic, &m->numeric, nullptr, nullptr);
	return (status == UMFPACK_OK);
#else
	return false;
#endif
}

bool UMFPACKSolver::BackSolve(double* x, double* y)
{
#ifdef HAS_SUITESPARSE
	SparseMatrix* A = m->A;
	if (A == nullptr) return false;

	int status = umfpack_dl_solve(UMFPACK_A, m->ptr.data(), m->ind.data(), A->Values(), x, y, m->numeric, nullptr, nullptr);
	return (status == UMFPACK_OK);
#else
	return false;
#endif
}

void UMFPACKSolver::Destroy()
{
#ifdef HAS_SUITESPARSE
	umfpack_di_free_symbolic(&m->symbolic);
	umfpack_di_free_numeric(&m->numeric);
#endif
}
