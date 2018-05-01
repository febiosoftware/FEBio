#include "stdafx.h"
#include "JFNKMatrix.h"
#include "FENewtonSolver.h"

JFNKMatrix::JFNKMatrix(FENewtonSolver* pns, SparseMatrix* K) : m_pns(pns), m_K(K)
{
	m_ndim = pns->m_neq;
	m_nsize = 0;

	// TODO: For contact problems we'll need some mechanism to change the array size
	m_v.resize(m_ndim);
	m_R.resize(m_ndim);
}

//! Create a sparse matrix from a sparse-matrix profile
void JFNKMatrix::Create(SparseMatrixProfile& MP) 
{ 
	m_K->Create(MP); 
	m_ndim = m_K->Size();
	m_nsize = m_K->NonZeroes(); 
}

void JFNKMatrix::mult_vector(double* x, double* r)
{
	double eps = 0.001;
	int neq = (int)m_pns->m_ui.size();

	for (int i = 0; i<neq; ++i) m_v[i] = eps*x[i];

	m_pns->Update(m_v);
	m_pns->Residual(m_R);

	for (int i = 0; i<neq; ++i)
	{
		r[i] = (m_pns->m_R0[i] - m_R[i]) / eps;
	}
}
