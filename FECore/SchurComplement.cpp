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
#include "SchurComplement.h"

SchurComplementA::SchurComplementA(LinearSolver* A, SparseMatrix* B, SparseMatrix* C, SparseMatrix* D)
{
	m_print_level = 0;

	m_A = A;
	m_B = B;
	m_C = C;
	m_D = D;

	int n0 = m_B->Rows();
	int n1 = m_B->Columns();
	assert(n0 == m_C->Columns());
	assert(n1 == m_C->Rows());

	m_tmp1.resize(n0, 0.0);
	m_tmp2.resize(n0, 0.0);
	m_tmp3.resize(n1, 0.0);

	m_nrow = n1;
	m_ncol = n1;

	m_bnegate = false;
}

// negate schur complement
void SchurComplementA::NegateSchur(bool b)
{
	m_bnegate = b;
}

// set the print level
void SchurComplementA::SetPrintLevel(int printLevel)
{
	m_print_level = printLevel;
}

//! multiply with vector
bool SchurComplementA::mult_vector(double* x, double* r)
{
	m_B->mult_vector(x, &m_tmp1[0]);

	if (m_print_level != 0) printf("backsolving in SchurComplement\n");
	if (m_A) { if (m_A->BackSolve(m_tmp2, m_tmp1) == false) return false; }
	else m_tmp2 = m_tmp1;
	m_C->mult_vector(&m_tmp2[0], r);

	if (m_D)
	{
		if (m_D->mult_vector(x, &m_tmp3[0]) == false) return false;

		size_t n = m_tmp3.size();
		for (size_t i = 0; i<n; ++i) r[i] -= m_tmp3[i];
	}

	if (m_bnegate)
	{
		size_t n = m_tmp3.size();
		for (int i = 0; i < n; ++i) r[i] = -r[i];
	}

	return true;
}


SchurComplementD::SchurComplementD(SparseMatrix* A, SparseMatrix* B, SparseMatrix* C, LinearSolver* D)
{
	m_print_level = 0;

	m_A = A;
	m_B = B;
	m_C = C;
	m_D = D;

	int n0 = m_B->Rows();
	int n1 = m_B->Columns();
	assert(n0 == m_C->Columns());
	assert(n1 == m_C->Rows());

	m_tmp1.resize(n1, 0.0);
	m_tmp2.resize(n1, 0.0);
	m_tmp3.resize(n0, 0.0);

	m_nrow = n0;
	m_ncol = n0;
}

// set the print level
void SchurComplementD::SetPrintLevel(int printLevel)
{
	m_print_level = printLevel;
}

//! multiply with vector
bool SchurComplementD::mult_vector(double* x, double* r)
{
	m_C->mult_vector(x, &m_tmp1[0]);

	if (m_print_level != 0) printf("backsolving in SchurComplement\n");
	if (m_D->BackSolve(m_tmp2, m_tmp1) == false) return false;
	m_B->mult_vector(&m_tmp2[0], r);

	if (m_A)
	{
		if (m_A->mult_vector(x, &m_tmp3[0]) == false) return false;

		size_t n = m_tmp3.size();
		for (size_t i = 0; i<n; ++i) r[i] -= m_tmp3[i];
	}

	return true;
}
