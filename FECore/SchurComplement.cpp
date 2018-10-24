#include "stdafx.h"
#include "SchurComplement.h"

SchurComplement::SchurComplement(LinearSolver* A, SparseMatrix* B, SparseMatrix* C, SparseMatrix* D)
{
	m_print_level = 0;

	m_A = A;
	m_B = B;
	m_C = C;
	m_D = D;

	int n0 = m_B->Columns();
	int n1 = m_B->Rows();
	assert(n0 == m_C->Rows());
	assert(n1 == m_C->Columns());

	m_tmp1.resize(n1, 0.0);
	m_tmp2.resize(n1, 0.0);

	m_nrow = n0;
	m_ncol = n0;
}

// set the print level
void SchurComplement::SetPrintLevel(int printLevel)
{
	m_print_level = printLevel;
}

//! multiply with vector
bool SchurComplement::mult_vector(double* x, double* r)
{
	m_B->mult_vector(x, &m_tmp1[0]);

	if (m_print_level != 0) printf("backsolving in SchurComplement\n");
	if (m_A->BackSolve(m_tmp2, m_tmp1) == false) return false;
	m_C->mult_vector(&m_tmp2[0], r);

	if (m_D)
	{
		if (m_D->mult_vector(x, &m_tmp1[0]) == false) return false;

		size_t n = m_tmp1.size();
		for (size_t i = 0; i<n; ++i) r[i] -= m_tmp1[i];
	}

	return true;
}
