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
#include "DenseMatrix.h"
using namespace FECore;
using namespace std;

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
void DenseMatrix::Assemble(const matrix& ke, const vector<int>& lm)
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
void DenseMatrix::Assemble(const matrix& ke, const vector<int>& LMi, const vector<int>& LMj)
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
