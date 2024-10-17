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
#include "CompactSymmMatrix64.h"
using namespace std;

CompactSymmMatrix64::CompactSymmMatrix64(int offset) : m_offset(offset)
{
	m_pvalues   = nullptr;
	m_pindices  = nullptr;
	m_ppointers = nullptr;
}

CompactSymmMatrix64::~CompactSymmMatrix64()
{
	Clear();
}

void CompactSymmMatrix64::Clear()
{
	delete[] m_pvalues; m_pvalues = nullptr;
	delete[] m_pindices; m_pindices = nullptr;
	delete[] m_ppointers; m_ppointers = nullptr;
}

void CompactSymmMatrix64::Create(SparseMatrixProfile& mp)
{
	// TODO: we should probably enforce that the matrix is square
	m_nrow = mp.Rows();
	m_ncol = mp.Columns();
	size_t nr = m_nrow;
	size_t nc = m_ncol;
	assert(nr == nc);

	Clear();

	// allocate pointers to column offsets
	m_ppointers = new long long[nc + 1];
	for (int i = 0; i <= nc; ++i) m_ppointers[i] = 0;

	size_t nsize = 0;
	for (size_t i = 0; i<nc; ++i)
	{
		SparseMatrixProfile::ColumnProfile& a = mp.Column(i);
		size_t n = a.size();
		for (size_t j = 0; j<n; j++)
		{
			size_t a0 = a[j].start;
			size_t a1 = a[j].end;

			// only grab lower-triangular
			if (a1 >= i)
			{
				if (a0 < i) a0 = i;
				size_t asize = a1 - a0 + 1;
				nsize += asize;
				m_ppointers[i] += asize;
			}
		}
	}

	// allocate indices which store row index for each matrix element
	m_pindices = new long long[nsize];
	size_t m = 0;
	for (size_t i = 0; i <= nc; ++i)
	{
		long long n = m_ppointers[i];
		m_ppointers[i] = m;
		m += n;
	}

	for (size_t i = 0; i<nc; ++i)
	{
		SparseMatrixProfile::ColumnProfile& a = mp.Column(i);
		size_t n = a.size();
		size_t nval = 0;
		for (size_t j = 0; j<n; j++)
		{
			size_t a0 = a[j].start;
			size_t a1 = a[j].end;

			// only grab lower-triangular
			if (a1 >= i)
			{
				if (a0 < i) a0 = i;
				for (size_t k = a0; k <= a1; ++k)
				{
					m_pindices[m_ppointers[i] + nval] = (long long)k;
					++nval;
				}
			}
		}
	}

	// offset the indicies for fortran arrays
	if (Offset())
	{
		for (size_t i = 0; i <= nc; ++i) m_ppointers[i]++;
		for (size_t i = 0; i<nsize; ++i) m_pindices[i]++;
	}

	// create the values array
	m_nsize = nsize;
	m_pvalues = new double[nsize];
}

//-----------------------------------------------------------------------------
// this sort function is defined in qsort.cpp
void qsort(int n, const int* arr, int* indx);

//-----------------------------------------------------------------------------
//! This function assembles the local stiffness matrix
//! into the global stiffness matrix which is in compact column storage
//!
void CompactSymmMatrix64::Assemble(const matrix& ke, const vector<int>& LM)
{
	// get the number of degrees of freedom
	const int N = ke.rows();

	// find the permutation array that sorts LM in ascending order
	// we can use this to speed up the row search (i.e. loop over n below)
	std::vector<int> P(N);
	qsort(N, &LM[0], &P[0]);

	// get the data pointers 
	long long* indices = Indices64();
	long long* pointers = Pointers64();
	double* pd = Values();

	// find the starting index
	size_t N0 = 0;
	while ((N0<N) && (LM[P[N0]]<0)) ++N0;

	// assemble element stiffness
	for (size_t m = N0; m<N; ++m)
	{
		size_t j = P[m];
		size_t J = LM[j];
		size_t n = 0;
		double* pm = pd + (pointers[J] - m_offset);
		long long* pi = indices + (pointers[J] - m_offset);
		long long l = pointers[J + 1] - pointers[J];
		size_t M0 = m;
		while ((M0>N0) && (LM[P[M0 - 1]] == J)) M0--;
		for (size_t k = M0; k<N; ++k)
		{
			size_t i = P[k];
			size_t I = LM[i] + (size_t)m_offset;
			for (; n<l; ++n)
				if (pi[n] == I)
				{
					#pragma omp atomic
					pm[n] += ke[i][j];
					break;
				}
		}
	}
}

void CompactSymmMatrix64::Assemble(const matrix& ke, const vector<int>& LMi, const vector<int>& LMj)
{
	const int N = ke.rows();
	const int M = ke.columns();

	long long* indices = Indices64();
	long long* pointers = Pointers64();
	double* values = Values();

	for (size_t i = 0; i<N; ++i)
	{
		int I = LMi[i];

		for (size_t j = 0; j<M; ++j)
		{
			int J = LMj[j];

			// only add values to lower-diagonal part of stiffness matrix
			if ((I >= J) && (J >= 0))
			{
				double* pv = values + (pointers[J] - m_offset);
				long long* pi = indices + (pointers[J] - m_offset);
				long long l = pointers[J + 1] - pointers[J];
				for (size_t n = 0; n<l; ++n) 
					if (pi[n] - m_offset == I)
					{
						#pragma omp atomic
						pv[n] += ke[i][j];
						break;
					}
			}
		}
	}
}

void CompactSymmMatrix64::add(int i, int j, double v)
{
	// only add to lower triangular part
	// since FEBio only works with the upper triangular part
	// we have to swap the indices
	i ^= j; j ^= i; i ^= j;

	if (j <= i)
	{
		double* pd = m_pvalues + (m_ppointers[j] - m_offset);
		long long* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		i += m_offset;
		long long n1 = m_ppointers[j + 1] - m_ppointers[j] - 1;
		long long n0 = 0;
		long long n = n1 / 2;
		do
		{
			long long m = pi[n];
			if (m == i)
			{
				#pragma omp atomic
				pd[n] += v;
				return;
			}
			else if (m < i)
			{
				n0 = n;
				n = (n0 + n1 + 1) >> 1;
			}
			else
			{
				n1 = n;
				n = (n0 + n1) >> 1;
			}
		} while (n0 != n1);

		assert(false);
	}
}

//! check fo a matrix item
bool CompactSymmMatrix64::check(int i, int j)
{
	// only the lower-triangular part is stored, so swap indices if necessary
	if (i<j) { i ^= j; j ^= i; i ^= j; }

	// only add to lower triangular part
	if (j <= i)
	{
		long long* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		long long l = m_ppointers[j + 1] - m_ppointers[j];
		for (long long n = 0; n<l; ++n)
			if (pi[n] == i + m_offset)
			{
				return true;
			}
	}

	return false;
}

//! set matrix item
void CompactSymmMatrix64::set(int i, int j, double v)
{
	if (j <= i)
	{
		long long* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		long long l = m_ppointers[j + 1] - m_ppointers[j];
		for (long long n = 0; n<l; ++n)
			if (pi[n] == i + m_offset)
			{
				long long k = m_ppointers[j] + n;
				k -= m_offset;
#pragma omp critical (CSM_set)
				m_pvalues[k] = v;
				return;
			}

		assert(false);
	}
}

//! get a matrix item
double CompactSymmMatrix64::get(int i, int j)
{
	// only the lower-triangular part is stored, so swap indices if necessary
	if (i<j) { i ^= j; j ^= i; i ^= j; }

	long long *pi = m_pindices + m_ppointers[j], k;
	pi -= m_offset;
	long long l = m_ppointers[j + 1] - m_ppointers[j];
	for (long long n = 0; n<l; ++n)
		if (pi[n] == i + m_offset)
		{
			k = m_ppointers[j] + n;
			k -= m_offset;
			return m_pvalues[k];
		}
	return 0;
}

void CompactSymmMatrix64::Zero()
{
	for (size_t i = 0; i < m_nsize; ++i) m_pvalues[i] = 0.0;
}
