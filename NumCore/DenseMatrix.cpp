#include "stdafx.h"
#include "DenseMatrix.h"

//-----------------------------------------------------------------------------
NumCore::DenseMatrix::DenseMatrix()
{
	m_pr = 0;
}

//-----------------------------------------------------------------------------
NumCore::DenseMatrix::~DenseMatrix()
{
	delete [] m_pd; m_pd = 0;
	delete [] m_pr; m_pr = 0;
}

//-----------------------------------------------------------------------------
// Creat a dense matrix of size N x N
void NumCore::DenseMatrix::Create(int N)
{
	if (N != m_ndim)
	{
		if (m_pd) delete [] m_pd;
		if (m_pr) delete [] m_pr;

		m_pd = new double[N*N];
		m_pr = new double*[N];

		for (int i=0; i<N; ++i) m_pr[i] = m_pd + i*N;

		m_ndim = N;
		m_nsize = N*N;
	}
}


//-----------------------------------------------------------------------------
//! This function assembles the local stiffness matrix
//! into the global stiffness matrix which is in dense format
//!
void NumCore::DenseMatrix::Assemble(matrix& ke, vector<int>& lm)
{
	int i, j, I, J;

	const int N = ke.rows();

	for (i=0; i<N; ++i)
	{
		if ((I = lm[i])>=0)
		{
			for (j=0; j<N; ++j)
			{
				if ((J = lm[j]) >= 0) m_pr[I][J] += ke[i][j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
void NumCore::DenseMatrix::Assemble(matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	int i, j, I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	for (i=0; i<N; ++i)
	{
		if ((I = LMi[i])>=0)
		{
			for (j=0; j<M; ++j)
			{
				if ((J = LMj[j]) >= 0) m_pr[I][J] += ke[i][j];
			}
		}
	}
}
