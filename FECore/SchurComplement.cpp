#include "stdafx.h"
#include "SchurComplement.h"

SchurComplement::SchurComplement(LinearSolver* A, SparseMatrix* B, SparseMatrix* C, SparseMatrix* D)
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
		if (m_D->mult_vector(x, &m_tmp3[0]) == false) return false;

		size_t n = m_tmp3.size();
		for (size_t i = 0; i<n; ++i) r[i] -= m_tmp3[i];
	}

	return true;
}


SchurComplement2::SchurComplement2(SparseMatrix* A, SparseMatrix* B, SparseMatrix* C, LinearSolver* D)
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
void SchurComplement2::SetPrintLevel(int printLevel)
{
	m_print_level = printLevel;
}

//! multiply with vector
bool SchurComplement2::mult_vector(double* x, double* r)
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
