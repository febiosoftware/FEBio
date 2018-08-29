#include "stdafx.h"
#include "DenseMatrix.h"

//-----------------------------------------------------------------------------
DenseMatrix::DenseMatrix()
{
	m_pr = 0;
	m_pd = 0;
}

//-----------------------------------------------------------------------------
DenseMatrix::~DenseMatrix()
{
	delete [] m_pd; m_pd = 0;
	delete [] m_pr; m_pr = 0;
}

//-----------------------------------------------------------------------------
void DenseMatrix::Zero()
{
	memset(m_pd, 0, m_nsize*sizeof(double)); 
}

//-----------------------------------------------------------------------------
void DenseMatrix::Clear()
{
	if (m_pd) delete [] m_pd; m_pd = 0;
	if (m_pr) delete [] m_pr; m_pr = 0;

	SparseMatrix::Clear();
}

//-----------------------------------------------------------------------------
void DenseMatrix::Create(SparseMatrixProfile& mp)
{
	Create(mp.Rows(), mp.Columns()); 
}

//-----------------------------------------------------------------------------
// Creat a dense matrix of size N x N
void DenseMatrix::Create(int rows, int cols)
{
	if ((rows != m_nrow) || (cols != m_ncol))
	{
		if (m_pd) delete [] m_pd;
		if (m_pr) delete [] m_pr;

		m_pd = new double[rows*cols];
		m_pr = new double*[rows];

		for (int i=0; i<rows; ++i) m_pr[i] = m_pd + i*cols;

		m_nrow = rows;
		m_ncol = cols;
		m_nsize = rows*cols;
	}
}

//-----------------------------------------------------------------------------
//! This function assembles the local stiffness matrix
//! into the global stiffness matrix which is in dense format
//!
void DenseMatrix::Assemble(matrix& ke, vector<int>& lm)
{
	int I, J;
	const int N = ke.rows();
	const int M = ke.columns();

	for (int i=0; i<N; ++i)
	{
		if ((I = lm[i])>=0)
		{
			for (int j=0; j<M; ++j)
			{
				if ((J = lm[j]) >= 0) m_pr[I][J] += ke[i][j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
void DenseMatrix::Assemble(matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	int I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	for (int i=0; i<N; ++i)
	{
		if ((I = LMi[i])>=0)
		{
			for (int j=0; j<M; ++j)
			{
				if ((J = LMj[j]) >= 0) m_pr[I][J] += ke[i][j];
			}
		}
	}
}
