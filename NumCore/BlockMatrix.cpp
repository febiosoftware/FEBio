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
#include "BlockMatrix.h"
#include <assert.h>

//-----------------------------------------------------------------------------
BlockMatrix::BlockMatrix()
{
}

//-----------------------------------------------------------------------------
BlockMatrix::~BlockMatrix()
{
	// clear memory allocations
	Clear();

	// delete blocks
	const int n = (int) m_Block.size();
	for (int i=0; i<n; ++i) delete m_Block[i].pA;
}

//-----------------------------------------------------------------------------
//! This function sets the partitions for the blocks. 
//! The partition list contains the number of rows (and cols) for each partition.
//! So for instance, if part = {10,10}, a 2x2=4 partition is created, where
//! each block is a 10x10 matrix.
//
//! TODO: I want to put the partition information in the matrix profile structure
//!       so that the Create function can be used to create all the blocks.
void BlockMatrix::Partition(const vector<int>& part, Matrix_Type mtype, int offset)
{
	// copy the partitions, but store equation numbers instead of number of equations
	const int n = (int)part.size();
	m_part.resize(n+1);
	m_part[0] = 0;
	for (int i=0; i<n; ++i) m_part[i+1] = m_part[i] + part[i];

	// create the block structure for all the partitions
	m_Block.resize(n*n);
	int nrow = 0;
	for (int i=0; i<n; ++i) // loop over rows
	{
		int ncol = 0;
		for (int j=0; j<n; ++j) // loop over cols
		{
			BLOCK& Bij = m_Block[i*n+j];
			Bij.nstart_row = nrow;
			Bij.nend_row   = Bij.nstart_row + part[i] - 1;

			Bij.nstart_col = ncol;
			Bij.nend_col   = Bij.nstart_col + part[j] - 1;

			// Note the parameters in the constructors.
			// This is because we are using Pardiso for this
			if (i==j)
			{
				if (mtype == REAL_SYMMETRIC)
					Bij.pA = new CompactSymmMatrix(offset);
				else
					Bij.pA = new CRSSparseMatrix(offset);
			}
			else
				Bij.pA = new CRSSparseMatrix(offset);
			
			ncol += part[j];
		}
		nrow += part[i];
	}
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix from a sparse-matrix profile
void BlockMatrix::Create(SparseMatrixProfile& MP)
{
	m_nrow = MP.Rows();
	m_ncol = MP.Columns();
	m_nsize = 0;
	const int N = (int) m_Block.size();
	for (int i=0; i<N; ++i)
	{
		BLOCK& Bi = m_Block[i];
		SparseMatrixProfile MPi = MP.GetBlockProfile(Bi.nstart_row, Bi.nstart_col, Bi.nend_row, Bi.nend_col);
		Bi.pA->Create(MPi);
		m_nsize += Bi.pA->NonZeroes();
	}
}

//-----------------------------------------------------------------------------
//! assemble a matrix into the sparse matrix
void BlockMatrix::Assemble(const matrix& ke, const std::vector<int>& lm)
{
	int I, J;
	const int N = ke.rows();
	for (int i=0; i<N; ++i)
	{
		if ((I = lm[i])>=0)
		{
			for (int j=0; j<N; ++j)
			{
				if ((J = lm[j]) >= 0) add(I,J, ke[i][j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! assemble a matrix into the sparse matrix
void BlockMatrix::Assemble(const matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj)
{
	int I, J;
	const int N = ke.rows();
	const int M = ke.columns();
	for (int i=0; i<N; ++i)
	{
		if ((I = lmi[i])>=0)
		{
			for (int j=0; j<M; ++j)
			{
				if ((J = lmj[j]) >= 0) add(I,J, ke[i][j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
// helper function for finding partitions
int BlockMatrix::find_partition(int i)
{
	const int N = (int)m_part.size() - 1;
	int n = 0;
	for (; n<N; ++n)
		if (m_part[n+1] > i) break;
	assert(n<N);
	return n;
}

//-----------------------------------------------------------------------------
// helper function for finding a block
BlockMatrix::BLOCK& BlockMatrix::Block(int i, int j)
{
	const int N = (int)m_part.size() - 1;
	return m_Block[i*N+j];
}

//-----------------------------------------------------------------------------
//! set entry to value
bool BlockMatrix::check(int i, int j)
{
	int nr = find_partition(i);
	int nc = find_partition(j);
	return Block(nr, nc).pA->check(i - m_part[nr], j - m_part[nc]);
}

//-----------------------------------------------------------------------------
//! set entry to value
void BlockMatrix::set(int i, int j, double v)
{
	int nr = find_partition(i);
	int nc = find_partition(j);
	Block(nr, nc).pA->set(i - m_part[nr], j - m_part[nc], v);
}

//-----------------------------------------------------------------------------
//! add value to entry
void BlockMatrix::add(int i, int j, double v)
{
	int nr = find_partition(i);
	int nc = find_partition(j);
	Block(nr, nc).pA->add(i - m_part[nr], j - m_part[nc], v);
}

//-----------------------------------------------------------------------------
//! retrieve value
double BlockMatrix::get(int i, int j)
{
	int nr = find_partition(i);
	int nc = find_partition(j);
	return Block(nr, nc).pA->get(i - m_part[nr], j - m_part[nc]);
}

//-----------------------------------------------------------------------------
//! get the diagonal value
double BlockMatrix::diag(int i)
{
	int n = find_partition(i);
	return Block(n, n).pA->diag(i - m_part[n]);
}

//-----------------------------------------------------------------------------
//! release memory for storing data
void BlockMatrix::Clear()
{
	// clear the blocks
	const int n = (int) m_Block.size();
	for (int i=0; i<n; ++i) m_Block[i].pA->Clear();
}

//-----------------------------------------------------------------------------
//! Zero all matrix elements
void BlockMatrix::Zero()
{
	// zero the blocks
	const int n = (int) m_Block.size();
	for (int i=0; i<n; ++i) m_Block[i].pA->Zero();
}

//-----------------------------------------------------------------------------
//! multiply with vector
bool BlockMatrix::mult_vector(double* x, double* r)
{
	int nr = Rows();
	vector<double> tmp(nr, 0);
	for (int i=0; i<nr; ++i) r[i] = 0.0;
	int NP = Partitions();
	for (int i=0; i<NP; ++i)
	{
		int n0 = m_part[i];
		for (int j=0; j<NP; ++j)
		{
			int m0 = m_part[j];

			BLOCK& bij = Block(i, j);

			bij.pA->mult_vector(x + m0, &tmp[0] + n0);

			int nj = bij.Rows();
			for (int k=0; k<nj; ++k) r[n0 + k] += tmp[n0 + k];
		}
	}

	return true;
}

//! row and column scale
void BlockMatrix::scale(const vector<double>& L, const vector<double>& R)
{
	vector<double> Li, Rj;
	for (int n = 0; n < m_Block.size(); ++n)
	{
		BLOCK& bn = m_Block[n];
		if (bn.pA)
		{
			int NR = bn.Rows();
			int NC = bn.Cols();

			Li.resize(NR);
			Rj.resize(NC);

			for (int i = 0; i < NR; ++i) Li[i] = L[i + bn.nstart_row];
			for (int i = 0; i < NC; ++i) Rj[i] = R[i + bn.nstart_col];

			bn.pA->scale(Li, Rj);
		}
	}
}
