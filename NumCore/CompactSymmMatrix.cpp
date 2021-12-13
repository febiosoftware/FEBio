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
#include "CompactSymmMatrix.h"
using namespace std;

//-----------------------------------------------------------------------------
//! constructor
CompactSymmMatrix::CompactSymmMatrix(int offset) : CompactMatrix(offset) 
{
	
}

//-----------------------------------------------------------------------------
bool CompactSymmMatrix::mult_vector(double* x, double* r)
{
	// get row count
	int N = Rows();
	int M = Columns();

	// zero result vector
	for (int j = 0; j<N; ++j) r[j] = 0.0;

	// loop over all columns
	for (int j = 0; j<M; ++j)
	{
		double* pv = m_pd + m_ppointers[j] - m_offset;
		int* pi = m_pindices + m_ppointers[j] - m_offset;
		int n = m_ppointers[j + 1] - m_ppointers[j];

		// add off-diagonal elements
		for (int i = 1; i<n - 7; i += 8)
		{
			// add lower triangular element
			r[pi[i    ] - m_offset] += pv[i    ] * x[j];
			r[pi[i + 1] - m_offset] += pv[i + 1] * x[j];
			r[pi[i + 2] - m_offset] += pv[i + 2] * x[j];
			r[pi[i + 3] - m_offset] += pv[i + 3] * x[j];
			r[pi[i + 4] - m_offset] += pv[i + 4] * x[j];
			r[pi[i + 5] - m_offset] += pv[i + 5] * x[j];
			r[pi[i + 6] - m_offset] += pv[i + 6] * x[j];
			r[pi[i + 7] - m_offset] += pv[i + 7] * x[j];
		}
		for (int i = 0; i<(n - 1) % 8; ++i)
			r[pi[n - 1 - i] - m_offset] += pv[n - 1 - i] * x[j];

		// add diagonal element
		double rj = pv[0] * x[j];

		// add upper-triangular elements
		for (int i = 1; i<n - 7; i += 8)
		{
			// add upper triangular element
			rj += pv[i    ] * x[pi[i    ] - m_offset];
			rj += pv[i + 1] * x[pi[i + 1] - m_offset];
			rj += pv[i + 2] * x[pi[i + 2] - m_offset];
			rj += pv[i + 3] * x[pi[i + 3] - m_offset];
			rj += pv[i + 4] * x[pi[i + 4] - m_offset];
			rj += pv[i + 5] * x[pi[i + 5] - m_offset];
			rj += pv[i + 6] * x[pi[i + 6] - m_offset];
			rj += pv[i + 7] * x[pi[i + 7] - m_offset];
		}
		for (int i = 0; i<(n - 1) % 8; ++i)
			rj += pv[n - 1 - i] * x[pi[n - 1 - i] - m_offset];

		r[j] += rj;
	}

	return true;
}

//-----------------------------------------------------------------------------
void CompactSymmMatrix::Create(SparseMatrixProfile& mp)
{
	// TODO: we should probably enforce that the matrix is square
	int nr = mp.Rows();
	int nc = mp.Columns();
	assert(nr==nc);

	// allocate pointers to column offsets
	int* pointers = new int[nc + 1];
	for (int i = 0; i <= nc; ++i) pointers[i] = 0;

	int nsize = 0;
	for (int i = 0; i<nc; ++i)
	{
		SparseMatrixProfile::ColumnProfile& a = mp.Column(i);
		int n = (int)a.size();
		for (int j = 0; j<n; j++)
		{
			int a0 = a[j].start;
			int a1 = a[j].end;

			// only grab lower-triangular
			if (a1 >= i)
			{
				if (a0 < i) a0 = i;
				int asize = a1 - a0 + 1;
				nsize += asize;
				pointers[i] += asize;
			}
		}
	}

	// allocate indices which store row index for each matrix element
	int* pindices = new int[nsize];
	int m = 0;
	for (int i = 0; i <= nc; ++i)
	{
		int n = pointers[i];
		pointers[i] = m;
		m += n;
	}

	for (int i = 0; i<nc; ++i)
	{
		SparseMatrixProfile::ColumnProfile& a = mp.Column(i);
		int n = (int)a.size();
		int nval = 0;
		for (int j = 0; j<n; j++)
		{
			int a0 = a[j].start;
			int a1 = a[j].end;

			// only grab lower-triangular
			if (a1 >= i)
			{
				if (a0 < i) a0 = i;
				for (int k = a0; k <= a1; ++k)
				{
					pindices[pointers[i] + nval] = k;
					++nval;
				}
			}
		}
	}

	// offset the indicies for fortran arrays
	if (Offset())
	{
		for (int i = 0; i <= nc; ++i) pointers[i]++;
		for (int i = 0; i<nsize; ++i) pindices[i]++;
	}

	// create the values array
	double* pvalues = new double[nsize];

	// create the stiffness matrix
	CompactMatrix::alloc(nr, nc, nsize, pvalues, pindices, pointers);
}

//-----------------------------------------------------------------------------
// this sort function is defined in qsort.cpp
void qsort(int n, const int* arr, int* indx);

//-----------------------------------------------------------------------------
//! This function assembles the local stiffness matrix
//! into the global stiffness matrix which is in compact column storage
//!
void CompactSymmMatrix::Assemble(const matrix& ke, const vector<int>& LM)
{
	// get the number of degrees of freedom
	const int N = ke.rows();

	// find the permutation array that sorts LM in ascending order
	// we can use this to speed up the row search (i.e. loop over n below)
	P.resize(N);
	qsort(N, &LM[0], &P[0]);

	// get the data pointers 
	int* indices = Indices();
	int* pointers = Pointers();
	double* pd = Values();
	int offset = Offset();

	// find the starting index
	int N0 = 0;
	while ((N0<N) && (LM[P[N0]]<0)) ++N0;

	// assemble element stiffness
	for (int m = N0; m<N; ++m)
	{
		int j = P[m];
		int J = LM[j];
		int n = 0;
		double* pm = pd + (pointers[J] - offset);
		int* pi = indices + (pointers[J] - offset);
		int l = pointers[J + 1] - pointers[J];
		int M0 = m;
		while ((M0>N0) && (LM[P[M0 - 1]] == J)) M0--;
		for (int k = M0; k<N; ++k)
		{
			int i = P[k];
			int I = LM[i] + offset;
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


//-----------------------------------------------------------------------------
void CompactSymmMatrix::Assemble(const matrix& ke, const vector<int>& LMi, const vector<int>& LMj)
{
	const int N = ke.rows();
	const int M = ke.columns();

	int* indices = Indices();
	int* pointers = Pointers();
	double* values = Values();

	for (int i = 0; i<N; ++i)
	{
		int I = LMi[i];

		for (int j = 0; j<M; ++j)
		{
			int J = LMj[j];

			// only add values to lower-diagonal part of stiffness matrix
			if ((I >= J) && (J >= 0))
			{
				double* pv = values + (pointers[J] - m_offset);
				int* pi = indices + (pointers[J] - m_offset);
				int l = pointers[J + 1] - pointers[J];
				for (int n = 0; n<l; ++n) 
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

//-----------------------------------------------------------------------------
//! add a matrix item
void CompactSymmMatrix::add(int i, int j, double v)
{
	// only add to lower triangular part
	// since FEBio only works with the upper triangular part
	// we have to swap the indices
	i ^= j; j ^= i; i ^= j;

	if (j <= i)
	{
		double* pd = m_pd + (m_ppointers[j] - m_offset);
		int* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		i += m_offset;
		int n1 = m_ppointers[j + 1] - m_ppointers[j] - 1;
		int n0 = 0;
		int n = n1 / 2;
		do
		{
			int m = pi[n];
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

//-----------------------------------------------------------------------------
//! check fo a matrix item
bool CompactSymmMatrix::check(int i, int j)
{
	// only the lower-triangular part is stored, so swap indices if necessary
	if (i<j) { i ^= j; j ^= i; i ^= j; }

	// only add to lower triangular part
	if (j <= i)
	{
		int* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		int l = m_ppointers[j + 1] - m_ppointers[j];
		for (int n = 0; n<l; ++n)
			if (pi[n] == i + m_offset)
			{
				return true;
			}
	}

	return false;
}


//-----------------------------------------------------------------------------
//! set matrix item
void CompactSymmMatrix::set(int i, int j, double v)
{
	if (j <= i)
	{
		int* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		int l = m_ppointers[j + 1] - m_ppointers[j];
		for (int n = 0; n<l; ++n)
			if (pi[n] == i + m_offset)
			{
				int k = m_ppointers[j] + n;
				k -= m_offset;
#pragma omp critical
				m_pd[k] = v;
				return;
			}

		assert(false);
	}
}

//-----------------------------------------------------------------------------
//! get a matrix item
double CompactSymmMatrix::get(int i, int j)
{
	// only the lower-triangular part is stored, so swap indices if necessary
	if (i<j) { i ^= j; j ^= i; i ^= j; }

	int *pi = m_pindices + m_ppointers[j], k;
	pi -= m_offset;
	int l = m_ppointers[j + 1] - m_ppointers[j];
	for (int n = 0; n<l; ++n)
		if (pi[n] == i + m_offset)
		{
			k = m_ppointers[j] + n;
			k -= m_offset;
			return m_pd[k];
		}
	return 0;
}

//-----------------------------------------------------------------------------
double CompactSymmMatrix::infNorm() const
{
	// get the matrix size
	const int N = Rows();

	// keep track of row sums
	vector<double> rowSums(N, 0.0);

	// loop over all columns
	for (int j = 0; j<N; ++j)
	{
		double* pv = m_pd + m_ppointers[j] - m_offset;
		int* pr = m_pindices + m_ppointers[j] - m_offset;
		int n = m_ppointers[j + 1] - m_ppointers[j];

		double ri = 0.0;
		for (int i = 0; i < n; ++i)
		{
			int irow = pr[i] - m_offset;
			double vij = fabs(pv[i]);
			ri += vij;
			if (irow != j) rowSums[irow] += vij;
		}

		rowSums[j] += ri;
	}

	// find the largest row sum
	double rmax = rowSums[0];
	for (int i = 1; i < N; ++i)
	{
		if (rowSums[i] > rmax) rmax = rowSums[i];
	}

	return rmax;
}

//-----------------------------------------------------------------------------
double CompactSymmMatrix::oneNorm() const
{
	// get the matrix size
	const int NR = Rows();
	const int NC = Columns();

	// keep track of col sums
	vector<double> colSums(NC, 0.0);

	// loop over all columns
	for (int j = 0; j<NC; ++j)
	{
		double* pv = m_pd + m_ppointers[j] - m_offset;
		int* pr = m_pindices + m_ppointers[j] - m_offset;
		int n = m_ppointers[j + 1] - m_ppointers[j];

		double cj = 0.0;
		for (int i = 0; i < n; ++i)
		{
			int irow = pr[i] - m_offset;
			double vij = fabs(pv[i]);
			cj += vij;
			if (irow != j) colSums[irow] += vij;
		}

		colSums[j] += cj;
	}

	// find the largest row sum
	double rmax = colSums[0];
	for (int i = 1; i < NC; ++i)
	{
		if (colSums[i] > rmax) rmax = colSums[i];
	}

	return rmax;
}

//-----------------------------------------------------------------------------
void CompactSymmMatrix::scale(const vector<double>& L, const vector<double>& R)
{
	// get the matrix size
	const int N = Columns();

	// loop over all columns
	for (int j = 0; j < N; ++j)
	{
		double* pv = m_pd + m_ppointers[j] - m_offset;
		int* pr = m_pindices + m_ppointers[j] - m_offset;
		int n = m_ppointers[j + 1] - m_ppointers[j];

		for (int i = 0; i < n; ++i)
		{
			pv[i] *= L[pr[i] - m_offset] * R[j];
		}
	}
}
