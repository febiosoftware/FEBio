#include "stdafx.h"
#include "Preconditioner.h"

REGISTER_SUPER_CLASS(Preconditioner, FEPRECONDITIONER_ID);

//=================================================================================================
Preconditioner::Preconditioner(FEModel* fem) : FECoreBase(fem, FEPRECONDITIONER_ID)
{
	m_K = nullptr;
}

Preconditioner::~Preconditioner()
{
}

// return the sparse matrix
SparseMatrix* Preconditioner::GetSparseMatrix()
{
	return m_K;
}

// set the sparse matrix
void Preconditioner::SetSparseMatrix(SparseMatrix* A)
{
	m_K = A;
}

//=================================================================================================
DiagonalPreconditioner::DiagonalPreconditioner(FEModel* fem) : Preconditioner(fem)
{
}

// create a preconditioner for a sparse matrix
bool DiagonalPreconditioner::Create()
{
	SparseMatrix* A = GetSparseMatrix();
	if (A == nullptr) return false;

	int N = A->Rows();
	if (A->Columns() != N) return false;

	m_D.resize(N);
	for (int i=0; i<N; ++i)
	{
		double dii = A->diag(i);
		if (dii == 0.0) return false;
		m_D[i] = 1.0 / dii;
	}

	return true;
}

// apply to vector P x = y
bool DiagonalPreconditioner::mult_vector(double* x, double* y)
{
	int N = (int)m_D.size();

#pragma omp parallel for
	for (int i=0; i<N; ++i)
	{
		y[i] = x[i]*m_D[i];
	}

	return true;
}
