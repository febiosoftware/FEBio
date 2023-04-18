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
#include "TestSolver.h"
#include <FECore/log.h>


BEGIN_FECORE_CLASS(TestSolver, LinearSolver)
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
TestSolver::TestSolver(FEModel* fem) : LinearSolver(fem), m_pA(0)
{
	m_mtype = -2;
}

//-----------------------------------------------------------------------------
TestSolver::~TestSolver()
{
	Destroy();
}

//-----------------------------------------------------------------------------
SparseMatrix* TestSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	// allocate the correct matrix format depending on matrix symmetry type
	switch (ntype)
	{
	case REAL_SYMMETRIC: m_mtype = -2; m_pA = new CompactSymmMatrix(1); break;
	case REAL_UNSYMMETRIC: m_mtype = 11; m_pA = new CRSSparseMatrix(1); break;
	case REAL_SYMM_STRUCTURE: m_mtype = 1; m_pA = new CRSSparseMatrix(1); break;
	default:
		assert(false);
		m_pA = nullptr;
	}

	return m_pA;
}

//-----------------------------------------------------------------------------
bool TestSolver::SetSparseMatrix(SparseMatrix* pA)
{
	if (m_pA) Destroy();
	m_pA = dynamic_cast<CompactMatrix*>(pA);
	m_mtype = -2;
	if (dynamic_cast<CRSSparseMatrix*>(pA)) m_mtype = 11;
	return (m_pA != nullptr);
}

//-----------------------------------------------------------------------------
bool TestSolver::PreProcess()
{
	m_n = m_pA->Rows();
	m_nnz = m_pA->NonZeroes();
	m_nrhs = 1;
	return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool TestSolver::Factor()
{
	return true;
}

//-----------------------------------------------------------------------------
bool TestSolver::BackSolve(double* x, double* b)
{
	// make sure we have work to do
	if (m_pA->Rows() == 0) return true;

	for (int i = 0; i < m_n; ++i) x[i] = 0.0;

	// update stats
	UpdateStats(1);

	return true;
}

//-----------------------------------------------------------------------------
void TestSolver::Destroy()
{

}
