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
#include "CompactUnSymmMatrix.h"
using namespace std;

//-----------------------------------------------------------------------------
// this sort function is defined in qsort.cpp
void qsort(int n, const int* arr, int* indx);

//=================================================================================================
CRSSparseMatrix::Iterator::Iterator(CRSSparseMatrix* A) : m_A(A)
{
	reset();
}

bool CRSSparseMatrix::Iterator::valid()
{
	return (n != -1);
}

void CRSSparseMatrix::Iterator::next()
{
	if (valid())
	{
		int* pr = m_A->Pointers();
		int l = pr[r+1] - pr[r];
		if (n < l-1) n++;
		else
		{
			r++;
			if (r >= m_A->Rows()) n = -1;
			else n = 0;
		}
	}
	else assert(false);
}

void CRSSparseMatrix::Iterator::reset()
{
	r = 0;
	n = 0;
	if (m_A == nullptr) n = -1;
}

MatrixItem CRSSparseMatrix::Iterator::get()
{
	assert(valid());
	int* pr = m_A->Pointers();
	int* pi = m_A->Indices() + (pr[r] - m_A->Offset());
	double* pv = m_A->Values() + (pr[r] - m_A->Offset());

	MatrixItem m;
	m.row = r;
	m.col = pi[n] - m_A->Offset();
	m.val = pv[n];

	return m;
}

void CRSSparseMatrix::Iterator::set(double v)
{
	assert(valid());
	int* pr = m_A->Pointers();
	double* pv = m_A->Values() + (pr[r] - m_A->Offset());
	pv[n] = v;
}


//=================================================================================================
// CRSSparseMatrix
//=================================================================================================

//-----------------------------------------------------------------------------
//! Constructor for CRSSparseMatrix class 
CRSSparseMatrix::CRSSparseMatrix(int offset) : CompactMatrix(offset)
{

}

//-----------------------------------------------------------------------------
CRSSparseMatrix::CRSSparseMatrix(const CRSSparseMatrix& A) : CompactMatrix(A.m_offset)
{
	m_nrow = A.m_nrow;
	m_ncol = A.m_ncol;
	m_nsize = A.m_nsize;

	m_ppointers = new int[m_nrow + 1];
	m_pindices = new int[m_nsize];
	m_pd = new double[m_nsize];

	for (int i = 0; i <= m_nrow; ++i) m_ppointers[i] = A.m_ppointers[i];
	for (int i = 0; i<m_nsize; ++i) m_pindices[i] = A.m_pindices[i];
	for (int i = 0; i<m_nsize; ++i) m_pd[i] = A.m_pd[i];
}

//-----------------------------------------------------------------------------
void CRSSparseMatrix::Create(SparseMatrixProfile& mp)
{
	int nr = mp.Rows();
	int nc = mp.Columns();

	int* pointers = new int[nr + 1];
	for (int i = 0; i <= nr; ++i) pointers[i] = 0;

	int nsize = 0;
	for (int i = 0; i<nc; ++i)
	{
		SparseMatrixProfile::ColumnProfile& a = mp.Column(i);
		int n = (int)a.size();
		for (int j = 0; j<n; j++)
		{
			int asize = a[j].end - a[j].start + 1;
			nsize += asize;
			for (int k = a[j].start; k <= a[j].end; ++k) pointers[k]++;
		}
	}

	int* pindices = new int[nsize];
	int m = 0;
	for (int i = 0; i <= nr; ++i)
	{
		int n = pointers[i];
		pointers[i] = m;
		m += n;
	}
	assert(pointers[nr] == nsize);

	vector<int> pval(nr, 0);
	for (int i = 0; i<nc; ++i)
	{
		SparseMatrixProfile::ColumnProfile& a = mp.Column(i);
		int n = (int)a.size();
		for (int j = 0; j<n; j++)
		{
			for (int k = a[j].start; k <= a[j].end; ++k)
			{
				pindices[pointers[k] + pval[k]] = i;
				++pval[k];
			}
		}
	}

	// offset the indicies for fortran arrays
	if (Offset())
	{
		for (int i = 0; i <= nr; ++i) pointers[i]++;
		for (int i = 0; i<nsize; ++i) pindices[i]++;
	}

	// create the values array
	double* pvalues = new double[nsize];

	// create the stiffness matrix
	CompactMatrix::alloc(nr, nc, nsize, pvalues, pindices, pointers);

	// calculate and print matrix bandwidth
//	feLog("\tMatrix bandwidth .......................... : %d\n", bandWidth());
}

void CRSSparseMatrix::Assemble(const matrix& ke, const vector<int>& LM)
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
	for (int k = N0; k<N; ++k)
	{
		int i = P[k];
		int I = LM[i];
		int n = 0;

		double* pm = pd + (pointers[I] - offset);
		int* pi = indices + (pointers[I] - offset);
		int l = pointers[I + 1] - pointers[I];

		for (int m = N0; m<N; ++m)
		{
			int j = P[m];
			int J = LM[j] + offset;
			double kij = ke[i][j];
			for (; n<l; ++n)
				if (pi[n] == J)
				{
#pragma omp atomic
					pm[n] += kij;
					break;
				}
		}
	}
}

//-----------------------------------------------------------------------------
void CRSSparseMatrix::Assemble(const matrix& ke, const vector<int>& LMi, const vector<int>& LMj)
{
	int I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	for (int i = 0; i<N; ++i)
	{
		if ((I = LMi[i]) >= 0)
		{
			for (int j = 0; j<M; ++j)
			{
				if ((J = LMj[j]) >= 0) add(I, J, ke[i][j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
// This algorithm uses a binary search for locating the correct row index
// This assumes that the indices are ordered!
void CRSSparseMatrix::add(int i, int j, double v)
{
	assert((i >= 0) && (i<m_nrow));
	assert((j >= 0) && (j<m_ncol));

	int* pi = m_pindices + (m_ppointers[i] - m_offset);
	double* pd = m_pd + (m_ppointers[i] - m_offset);
	int n1 = m_ppointers[i + 1] - m_ppointers[i] - 1;
	int n0 = 0;
	int n = n1 / 2;
	j += m_offset;
	do
	{
		int m = pi[n];
		if (m == j)
		{
#pragma omp atomic
			pd[n] += v;
			return;
		}
		else if (m < j)
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


//-----------------------------------------------------------------------------
void CRSSparseMatrix::set(int i, int j, double v)
{
	int* pi = m_pindices + (m_ppointers[i] - m_offset);
	int l = m_ppointers[i + 1] - m_ppointers[i];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == j + m_offset)
		{
#pragma omp critical
			m_pd[m_ppointers[i] + n - m_offset] = v;
			return;
		}
	}
	assert(false);
}

//-----------------------------------------------------------------------------
double CRSSparseMatrix::get(int i, int j)
{
	int* pi = m_pindices + (m_ppointers[i] - m_offset);
	int l = m_ppointers[i + 1] - m_ppointers[i];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == j + m_offset) return m_pd[m_ppointers[i] + n - m_offset];
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool CRSSparseMatrix::check(int i, int j)
{
	int* pi = m_pindices + (m_ppointers[i] - m_offset);
	int l = m_ppointers[i + 1] - m_ppointers[i];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == j + m_offset) return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
double CRSSparseMatrix::diag(int i)
{
	int* pi = m_pindices + (m_ppointers[i] - m_offset);
	int l = m_ppointers[i + 1] - m_ppointers[i];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == i + m_offset)
		{
			return m_pd[m_ppointers[i] + n - m_offset];
		}
	}
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
bool CRSSparseMatrix::mult_vector(double* x, double* r)
{
	// get the matrix size
	const int N = Rows();

	// loop over all rows
	#pragma omp parallel for schedule(guided)
	for (int i = 0; i < N; ++i)
	{
		const double* pv = m_pd + (m_ppointers[i] - m_offset);
		const int* pi = m_pindices + (m_ppointers[i] - m_offset);
		const int n = m_ppointers[i + 1] - m_ppointers[i];
		r[i] = 0.0;
		for (int j = 0; j < n; j ++)
		{
			r[i] += (*pv++) * x[*pi++ - m_offset];
		}
	}

	return true;
}

//! calculate the abs row sum 
double CRSSparseMatrix::infNorm() const
{
	// get the matrix size
	const int N = Rows();

	double norm = 0.0;
	// loop over all rows
	for (int i = 0; i<N; ++i)
	{
		double ri = 0.0;
		double* pv = m_pd + m_ppointers[i] - m_offset;
		int* pi = m_pindices + m_ppointers[i] - m_offset;
		int n = m_ppointers[i + 1] - m_ppointers[i];
		for (int j = 0; j<n; ++j) ri += fabs(pv[j]);

		if (ri > norm) norm = ri;
	}

	return norm;
}

//! calculate the one norm
double CRSSparseMatrix::oneNorm() const
{
	// get the matrix size
	const int NR = Rows();
	const int NC = Columns();

	vector<double> colNorms(NC, 0.0);

	// loop over all rows
	for (int i = 0; i<NR; ++i)
	{
		double ri = 0.0;
		double* pv = m_pd + m_ppointers[i] - m_offset;
		int* pi = m_pindices + m_ppointers[i] - m_offset;
		int n = m_ppointers[i + 1] - m_ppointers[i];
		for (int j = 0; j<n; ++j) colNorms[pi[j]-m_offset] += fabs(pv[j]);
	}

	// find max value
	double rmax = 0;
	for (int i = 0; i < NC; ++i)
	{
		if (colNorms[i] > rmax) rmax = colNorms[i];
	}

	return rmax;
}

//! make the matrix a unit matrix (retains sparsity pattern)
void CRSSparseMatrix::makeUnit()
{
	// loop over all rows
	const int N = Rows();
	for (int i = 0; i<N; ++i)
	{
		double* pv = m_pd + m_ppointers[i] - m_offset;
		int* pi = m_pindices + m_ppointers[i] - m_offset;
		int n = m_ppointers[i + 1] - m_ppointers[i];
		for (int j = 0; j < n; ++j)
		{
			if (pi[j] - m_offset == i) pv[j] = 1.0;
			else pv[j] = 0.0;
		}
	}
}

void CRSSparseMatrix::scale(double s)
{
	int N = NonZeroes();
	for (int i = 0; i < N; ++i) m_pd[i] *= s;
}

void CRSSparseMatrix::scale(const vector<double>& L, const vector<double>& R)
{
	const int N = Rows();
	assert(L.size() == Rows());
	assert(R.size() == Columns());
	for (int i = 0; i<N; ++i)
	{
		double* pv = m_pd + m_ppointers[i] - m_offset;
		int* pi = m_pindices + m_ppointers[i] - m_offset;
		int n = m_ppointers[i + 1] - m_ppointers[i];
		for (int j = 0; j<n; ++j)
		{
			pv[j] *= L[i] * R[pi[j] - m_offset];
		}
	}
}

//! extract a block of this matrix
void CRSSparseMatrix::get(int i0, int j0, int nr, int nc, CSRMatrix& M)
{
	// create the matrix
	M.create(nr, nc, m_offset);

	vector<double>& val = M.values();
	vector<int>& ind = M.indices();
	vector<int>& pnt = M.pointers(); assert(pnt.size() == nr + 1);

	// count how many values we'll need to copy
	int nnz = 0;
	for (int i = i0; i<i0 + nr; ++i)
	{
		int* pi = m_pindices + m_ppointers[i] - m_offset;
		int n = m_ppointers[i + 1] - m_ppointers[i];
		for (int j = 0; j<n; ++j)
		{
			int colj = pi[j] - m_offset;
			if ((colj >= j0) && (colj < j0 + nc)) nnz++;
		}
		pnt[i - i0 + 1] = nnz + m_offset;
	}

	// allocate arrays
	val.resize(nnz);
	ind.resize(nnz);

	// copy data
	nnz = 0;
	for (int i = i0; i<i0 + nr; ++i)
	{
		double* pv = m_pd + m_ppointers[i] - m_offset;
		int* pi = m_pindices + m_ppointers[i] - m_offset;
		int n = m_ppointers[i + 1] - m_ppointers[i];
		for (int j = 0; j<n; ++j)
		{
			int colj = pi[j] - m_offset;
			if ((colj >= j0) && (colj < j0 + nc))
			{
				val[nnz] = pv[j];
				ind[nnz] = colj - j0 + m_offset;
				nnz++;
			}
		}
	}
}

//=================================================================================================
// CCSSparseMatrix
//=================================================================================================

//-----------------------------------------------------------------------------
//! Constructor for CCSSparseMatrix class 
CCSSparseMatrix::CCSSparseMatrix(int offset) : CompactMatrix(offset)
{

}

//-----------------------------------------------------------------------------
CCSSparseMatrix::CCSSparseMatrix(const CCSSparseMatrix& A) : CompactMatrix(A.m_offset)
{
	m_nrow = A.m_nrow;
	m_ncol = A.m_ncol;
	m_nsize = A.m_nsize;

	m_ppointers = new int[m_ncol + 1];
	m_pindices = new int[m_nsize];
	m_pd = new double[m_nsize];

	for (int i = 0; i <= m_ncol; ++i) m_ppointers[i] = A.m_ppointers[i];
	for (int i = 0; i<m_nsize; ++i) m_pindices[i] = A.m_pindices[i];
	for (int i = 0; i<m_nsize; ++i) m_pd[i] = A.m_pd[i];
}

//-----------------------------------------------------------------------------
void CCSSparseMatrix::Create(SparseMatrixProfile& mp)
{
	int nr = mp.Rows();
	int nc = mp.Columns();

	int* pointers = new int[nc + 1];
	for (int i = 0; i <= nc; ++i) pointers[i] = 0;

	int nsize = 0;
	for (int i = 0; i<nc; ++i)
	{
		SparseMatrixProfile::ColumnProfile& a = mp.Column(i);
		int n = (int)a.size();
		for (int j = 0; j<n; j++)
		{
			int asize = a[j].end - a[j].start + 1;
			nsize += asize;
			pointers[i] += asize;
		}
	}

	int* pindices = new int[nsize];
	int m = 0;
	for (int i = 0; i <= nc; ++i)
	{
		int n = pointers[i];
		pointers[i] = m;
		m += n;
	}
	assert(pointers[nc] == nsize);

	for (int i = 0; i<nc; ++i)
	{
		SparseMatrixProfile::ColumnProfile& a = mp.Column(i);
		int n = (int)a.size();
		int nval = 0;
		for (int j = 0; j<n; j++)
		{
			for (int k = a[j].start; k <= a[j].end; ++k)
			{
				pindices[pointers[i] + nval] = k;
				nval++;
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
void CCSSparseMatrix::Assemble(const matrix& ke, const vector<int>& LM)
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

	for (int m = N0; m<N; ++m)
	{
		int j = P[m];
		int J = LM[j];
		int n = 0;
		double* pm = pd + (pointers[J] - offset);
		int* pi = indices + (pointers[J] - offset);
		int l = pointers[J + 1] - pointers[J];
		for (int k = N0; k<N; ++k)
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
void CCSSparseMatrix::Assemble(const matrix& ke, const vector<int>& LMi, const vector<int>& LMj)
{
	int I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	for (int i = 0; i<N; ++i)
	{
		if ((I = LMi[i]) >= 0)
		{
			for (int j = 0; j<M; ++j)
			{
				if ((J = LMj[j]) >= 0) add(I, J, ke[i][j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
// This algorithm uses a binary search for locating the correct row index
// This assumes that the indices are ordered!
void CCSSparseMatrix::add(int i, int j, double v)
{
	assert((i >= 0) && (i<m_nrow));
	assert((j >= 0) && (j<m_ncol));

	int* pi = m_pindices + (m_ppointers[j] - m_offset);
	i += m_offset;
	double* pd = m_pd + (m_ppointers[j] - m_offset);
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

//-----------------------------------------------------------------------------
void CCSSparseMatrix::set(int i, int j, double v)
{
	int* pi = m_pindices + (m_ppointers[j] - m_offset);
	int l = m_ppointers[j + 1] - m_ppointers[j];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == i + m_offset)
		{
#pragma omp critical
			m_pd[m_ppointers[j] + n - m_offset] = v;
			return;
		}
	}
	assert(false);
}

//-----------------------------------------------------------------------------
double CCSSparseMatrix::get(int i, int j)
{
	int* pi = m_pindices + (m_ppointers[j] - m_offset);
	int l = m_ppointers[j + 1] - m_ppointers[j];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == i + m_offset) return m_pd[m_ppointers[j] + n - m_offset];
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool CCSSparseMatrix::check(int i, int j)
{
	int* pi = m_pindices + m_ppointers[j] - m_offset;
	int l = m_ppointers[j + 1] - m_ppointers[j];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == i + m_offset) return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
double CCSSparseMatrix::diag(int i)
{
	int* pi = m_pindices + m_ppointers[i] - m_offset;
	int l = m_ppointers[i + 1] - m_ppointers[i];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == i + m_offset)
		{
			return m_pd[m_ppointers[i] + n - m_offset];
		}
	}
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
bool CCSSparseMatrix::mult_vector(double* x, double* r)
{
	// get the matrix size
	const int N = Rows();
	const int M = Columns();

	// zero r
	for (int i=0; i<N; ++i) r[i] = 0.0;

	// loop over all columns
	for (int i = 0; i<M; ++i)
	{
		double* pv = m_pd + (m_ppointers[i] - m_offset);
		int* pi = m_pindices + (m_ppointers[i] - m_offset);
		int n = m_ppointers[i + 1] - m_ppointers[i];
		for (int j = 0; j<n; j++)  r[pi[j] - m_offset] += pv[j] * x[i];
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculate the inf norm
double CCSSparseMatrix::infNorm() const
{
	// get the matrix size
	const int NR = Rows();
	const int NC = Columns();

	// keep track of row sums
	vector<double> rowSums(NR, 0.0);

	// loop over all columns
	for (int j = 0; j<NC; ++j)
	{
		double* pv = m_pd + m_ppointers[j] - m_offset;
		int* pr = m_pindices + m_ppointers[j] - m_offset;
		int n = m_ppointers[j + 1] - m_ppointers[j];

		for (int i = 0; i < n; ++i)
		{
			int irow = pr[i] - m_offset;
			double vij = fabs(pv[i]);
			rowSums[irow] += vij;
		}
	}

	// find the largest row sum
	double rmax = rowSums[0];
	for (int i = 1; i < NR; ++i)
	{
		if (rowSums[i] > rmax) rmax = rowSums[i];
	}

	return rmax;
}

//-----------------------------------------------------------------------------
//! calculate the one norm
double CCSSparseMatrix::oneNorm() const
{
	// get the matrix size
	const int NR = Rows();
	const int NC = Columns();

	// max col sum
	double cmax = 0.0;

	// loop over all columns
	for (int j = 0; j<NC; ++j)
	{
		double* pv = m_pd + m_ppointers[j] - m_offset;
		int* pr = m_pindices + m_ppointers[j] - m_offset;
		int n = m_ppointers[j + 1] - m_ppointers[j];

		double cj = 0.0;
		for (int i = 0; i < n; ++i)
		{
			double vij = fabs(pv[i]);
			cj += vij;
		}

		if (cj > cmax) cmax = cj;
	}

	return cmax;
}

//-----------------------------------------------------------------------------
void CCSSparseMatrix::scale(const vector<double>& L, const vector<double>& R)
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

//-----------------------------------------------------------------------------
//! Create a copy of the matrix (does not copy values)
CRSSparseMatrix* CRSSparseMatrix::Copy(int offset)
{
	CRSSparseMatrix* A = new CRSSparseMatrix(offset);

	int nnz = NonZeroes();
	int nrow = Rows();
	int ncol = Columns();

	int nn = (isRowBased() ? nrow : ncol);

	// allocate memory for copy
	double* vd = new double[nnz];
	int* id = new int[nnz];
	int* pd = new int[nn + 1];
	A->alloc(nrow, ncol, nnz, vd, id, pd);

	int offs = Offset();

	// copy indices
	int* is = Indices();
	for (int i = 0; i < nnz; ++i) id[i] = (is[i] - offs) + offset;

	// copy pointers
	int* ps = Pointers();
	for (int i = 0; i < nn + 1; ++i) pd[i] = (ps[i] - offs) + offset;

	return A;
}

//-----------------------------------------------------------------------------
//! Copy the values from another matrix
void CRSSparseMatrix::CopyValues(CompactMatrix* A)
{
	assert(NonZeroes() == A->NonZeroes());
	memcpy(Values(), A->Values(), sizeof(double)*NonZeroes());
}

//-----------------------------------------------------------------------------
//! convert to another format (currently only offset can be changed)
bool CRSSparseMatrix::Convert(int newOffset)
{
	if (newOffset == m_offset) return true;

	int d_off = newOffset - m_offset;

	int nn = Rows();
	int nnz = NonZeroes();
	m_offset = newOffset;

	// copy indices
	int* is = Indices();
	for (int i = 0; i < nnz; ++i) is[i] += d_off;

	// copy pointers
	int* ps = Pointers();
	for (int i = 0; i < nn + 1; ++i) ps[i] += d_off;

	return true;
}
