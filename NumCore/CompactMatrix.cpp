#include "stdafx.h"
#include "CompactMatrix.h"
#include <assert.h>

//=============================================================================
// CompactMatrix
//=============================================================================

//-----------------------------------------------------------------------------
CompactMatrix::~CompactMatrix()
{ 
	Clear(); 
}

//-----------------------------------------------------------------------------
void CompactMatrix::zero()
{
	memset(m_pd, 0, m_nsize*sizeof(double)); 
}

//-----------------------------------------------------------------------------
void CompactMatrix::Clear()
{
	if (m_pd) delete [] m_pd; m_pd = 0;
	if (m_pindices) delete [] m_pindices; m_pindices = 0;
	if (m_ppointers) delete [] m_ppointers; m_ppointers = 0;
}

//-----------------------------------------------------------------------------
CompactMatrix::CompactMatrix(int offset)
{
	m_pd = 0;
	m_pindices = 0;
	m_ppointers = 0;
	m_offset = offset;
}

//-----------------------------------------------------------------------------
void CompactMatrix::alloc(int N, int nz, double* pv, int* pi, int* pp)
{
	if (m_pd  ) delete [] m_pd; m_pd = pv;
	if (m_pindices ) delete [] m_pindices; m_pindices = pi;
	if (m_ppointers) delete [] m_ppointers; m_ppointers = pp;

	m_ndim  = N;
	m_nsize = nz;
}

//=============================================================================
// CompactSymmMatrix
//=============================================================================

//-----------------------------------------------------------------------------
//! constructor
CompactSymmMatrix::CompactSymmMatrix(int offset) : CompactMatrix(offset) {}

//-----------------------------------------------------------------------------
void CompactSymmMatrix::mult_vector(const vector<double>& x, vector<double>& r)
{
	assert(Size() == (int) x.size());
	assert(Size() == (int) r.size());
	mult_vector(&x[0], &r[0]);
}

//-----------------------------------------------------------------------------
void CompactSymmMatrix::mult_vector(const double* x, double* r)
{
	int j, i, n;
	int N = Size();

	double* pv, rj;
	int* pi;

	// loop over all columns
	for (j=0; j<N; ++j)
	{
		pv = m_pd  + m_ppointers[j] - m_offset;
		pi = m_pindices + m_ppointers[j]  - m_offset;
		n = m_ppointers[j+1] - m_ppointers[j];

		// add off-diagonal elements
		for (i=1; i<n-7; i+=8)
		{
			// add lower triangular element
			r[pi[i  ] - m_offset] += pv[i  ]*x[j];
			r[pi[i+1] - m_offset] += pv[i+1]*x[j];
			r[pi[i+2] - m_offset] += pv[i+2]*x[j];
			r[pi[i+3] - m_offset] += pv[i+3]*x[j];
			r[pi[i+4] - m_offset] += pv[i+4]*x[j];
			r[pi[i+5] - m_offset] += pv[i+5]*x[j];
			r[pi[i+6] - m_offset] += pv[i+6]*x[j];
			r[pi[i+7] - m_offset] += pv[i+7]*x[j];
		}
		for (i=0; i<(n-1)%8; ++i)
			r[pi[n-1-i] - m_offset] += pv[n-1-i]*x[j];

		// add diagonal element
		rj = pv[0]*x[j]; 

		// add upper-triangular elements
		for (i=1; i<n-7; i+=8)
		{
			// add upper triangular element
			rj += pv[i  ]*x[pi[i  ] - m_offset];
			rj += pv[i+1]*x[pi[i+1] - m_offset];
			rj += pv[i+2]*x[pi[i+2] - m_offset];
			rj += pv[i+3]*x[pi[i+3] - m_offset];
			rj += pv[i+4]*x[pi[i+4] - m_offset];
			rj += pv[i+5]*x[pi[i+5] - m_offset];
			rj += pv[i+6]*x[pi[i+6] - m_offset];
			rj += pv[i+7]*x[pi[i+7] - m_offset];
		}
		for (i=0; i<(n-1)%8; ++i)
			rj += pv[n-1-i]*x[pi[n-1-i] - m_offset];

		r[j] += rj;
	}
}


//-----------------------------------------------------------------------------
void CompactSymmMatrix::Create(SparseMatrixProfile& mp)
{
	int i, j, k, n;

	int neq = mp.size();

	int* pointers = new int[neq + 1];
	for (i=0; i<=neq; ++i) pointers[i] = 0;

	int nsize = 0;
	for (i=0; i<neq; ++i)
	{
		vector<int>& a = mp.column(i);
		n = a.size();
		for (j=0; j<n; j += 2)
		{
			nsize += a[j+1] - a[j] + 1;
			for (k=a[j]; k<=a[j+1]; ++k) pointers[k]++;
		}
	}

	int* pindices = new int[nsize];
	int m = 0;
	for (i=0; i<=neq; ++i)
	{
		n = pointers[i];
		pointers[i] = m;
		m += n;
	}

	int* pval = new int[neq];
	for (i=0; i<neq; ++i) pval[i] = 0;

	for (i=0; i<neq; ++i)
	{
		vector<int>& a = mp.column(i);
		n = a.size();
		for (j=0; j<n; j += 2)
		{
			for (k=a[j]; k<=a[j+1]; ++k)
			{
				pindices[ pointers[k] + pval[k]] = i;
				++pval[k];
			}
		}
	}

	// cleanup
	delete [] pval;

	// offset the indicies for fortran arrays
	if(Offset())
	{
		for (i=0; i<=neq; ++i) pointers[i]++;
		for (i=0; i<nsize; ++i) pindices[i]++;
	}

	// create the values array
	double* pvalues = new double[nsize];
	if (pvalues == 0) throw MemException(sizeof(double)*nsize);

	// create the stiffness matrix
	CompactMatrix::alloc(neq, nsize, pvalues, pindices, pointers);
}

//-----------------------------------------------------------------------------
// this sort function is defined in qsort.cpp
void qsort(int n, int* arr, int* indx);

//-----------------------------------------------------------------------------
//! This function assembles the local stiffness matrix
//! into the global stiffness matrix which is in compact column storage
//!
void CompactSymmMatrix::Assemble(matrix& ke, vector<int>& LM)
{
	// get the number of degrees of freedom
	const int N = ke.rows();

	// find the permutation array that sorts LM in ascending order
	// we can use this to speed up the row search (i.e. loop over n below)
	static vector<int> P; P.resize(N);
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
	for (int m=N0; m<N; ++m)
	{
		int j = P[m];
		int J = LM[j];
		int n = 0;
		double* pm = pd+pointers[J]-offset;
		int* pi = indices + pointers[J] - offset;
		int l = pointers[J+1] - pointers[J];
		int M0 = m;
		while ((M0>N0) && (LM[P[M0-1]] == J)) M0--;
		for (int k=M0; k<N; ++k)
		{
			int i = P[k];
			int I = LM[i] + offset;
			for (;n<l; ++n) 
				if (pi[n] == I)
				{
					pm[n] += ke[i][j];
					break;
				}
		}
	}
}


//-----------------------------------------------------------------------------
void CompactSymmMatrix::Assemble(matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	int i, j, I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	int* indices = Indices();
	int* pointers = Pointers();
	double* pd = Values();

	int *pi, l, n;

	for (i=0; i<N; ++i)
	{
		I = LMi[i];

		for (j=0; j<M; ++j)
		{
			J = LMj[j];

			// only add values to lower-diagonal part of stiffness matrix
			if ((I>=J) && (J>=0))
			{
				pi = indices + pointers[J];
				l = pointers[J+1] - pointers[J];
				for (n=0; n<l; ++n) if (pi[n] == I)
				{
					pd[pointers[J] + n] += ke[i][j];
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

	if (j<=i)
	{
		int* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		int l = m_ppointers[j+1] - m_ppointers[j];
		for (int n=0; n<l; ++n)
			if (pi[n] == i + m_offset)
			{
				m_pd[ m_ppointers[j] + n - m_offset ] += v;
				return;
			}

		assert(false);
	}
}

//-----------------------------------------------------------------------------
//! set matrix item
void CompactSymmMatrix::set(int i, int j, double v)
{
	int k;

	// only add to lower triangular part
	// since FEBio only works with the upper triangular part
	// we have to swap the indices
	i ^= j; j ^= i; i ^= j;

	if (j<=i)
	{
		int* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		int l = m_ppointers[j+1] - m_ppointers[j];
		for (int n=0; n<l; ++n)
			if (pi[n] == i + m_offset)
			{
				k = m_ppointers[j] + n;
				k -= m_offset;
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
	if (j>i) { i ^= j; j ^= i; i ^= j; }

	int *pi = m_pindices + m_ppointers[j], k;
	pi -= m_offset;
	int l = m_ppointers[j+1] - m_ppointers[j];
	for (int n=0; n<l; ++n)
		if (pi[n] == i + m_offset)
		{
			k = m_ppointers[j] + n;
			k -= m_offset;
			return m_pd[k];
		}
	return 0;
}

//=============================================================================
// CompactUnSymmMatrix
//=============================================================================

//-----------------------------------------------------------------------------
//! Constructor for CompactUnSymmMatrix class 
CompactUnSymmMatrix::CompactUnSymmMatrix(int offset, bool row_based) : CompactMatrix(offset)
{
	m_brow_based = row_based;
}

//-----------------------------------------------------------------------------
void CompactUnSymmMatrix::Create(SparseMatrixProfile& mp)
{
	int i, j, k, n;

	int neq = mp.size();

	int* pointers = new int[neq + 1];
	for (i=0; i<=neq; ++i) pointers[i] = 0;

	int nsize = 0;
	for (i=0; i<neq; ++i)
	{
		vector<int>& a = mp.column(i);
		n = a.size();
		for (j=0; j<n; j += 2)
		{
			nsize += 2*(a[j+1] - a[j] + 1);
			pointers[i] += a[j+1] - a[j] + 1;
			for (k=a[j]; k<=a[j+1]; ++k) pointers[k]++;
		}
		--pointers[i]; // we double counted the diagonal
		--nsize;
	}

	int* pindices = new int[nsize];
	int m = 0;
	for (i=0; i<=neq; ++i)
	{
		n = pointers[i];
		pointers[i] = m;
		m += n;
	}
	assert(pointers[neq] == nsize);

	int* pval = new int[neq];
	for (i=0; i<neq; ++i) pval[i] = 0;

	for (i=0; i<neq; ++i)
	{
		vector<int>& a = mp.column(i);
		n = a.size();
		for (j=0; j<n; j += 2)
		{
			for (k=a[j]; k<=a[j+1]; ++k)
			{
				pindices[ pointers[i] + pval[i]] = k;
				++pval[i];
			}
		}
		for (j=0; j<n; j += 2)
		{
			for (k=a[j]; k<=a[j+1]; ++k)
			{
				if (k != i)
				{
					pindices[ pointers[k] + pval[k]] = i;
					++pval[k];
				}
			}
		}
	}

	// cleanup
	delete [] pval;

	// offset the indicies for fortran arrays
	if(Offset())
	{
		for (i=0; i<=neq; ++i) pointers[i]++;
		for (i=0; i<nsize; ++i) pindices[i]++;
	}

	// create the values array
	double* pvalues = new double[nsize];
	if (pvalues == 0) throw MemException(sizeof(double)*nsize);

	// create the stiffness matrix
	CompactMatrix::alloc(neq, nsize, pvalues, pindices, pointers);
}

//-----------------------------------------------------------------------------
//! This function assembles the local stiffness matrix
//! into the global stiffness matrix which is in compact column storage and
//! the matrix is unsymmetric
//!

/*
void CompactUnSymmMatrix::Assemble(matrix& ke, vector<int>& LM)
{
	int i, j, I, J;

	const int N = ke.rows();

	for (i=0; i<N; ++i)
	{
		if ((I = LM[i])>=0)
		{
			for (j=0; j<N; ++j)
			{
				if ((J = LM[j]) >= 0) add(I,J, ke[i][j]);
			}
		}
	}
}
*/

void CompactUnSymmMatrix::Assemble(matrix& ke, vector<int>& LM)
{
	// get the number of degrees of freedom
	const int N = ke.rows();

	// find the permutation array that sorts LM in ascending order
	// we can use this to speed up the row search (i.e. loop over n below)
	static vector<int> P; P.resize(N);
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
	if (m_brow_based)
	{
		for (int m = N0; m<N; ++m)
		{
			int j = P[m];
			int J = LM[j];
			int n = 0;
			double* pm = pd + pointers[J] - offset;
			int* pi = indices + pointers[J] - offset;
			int l = pointers[J + 1] - pointers[J];
			for (int k = N0; k<N; ++k)
			{
				int i = P[k];
				int I = LM[i] + offset;
				for (; n<l; ++n)
				if (pi[n] == I)
				{
					pm[n] += ke[j][i];	// (i,j) is swapped since the matrix is row-based
					break;
				}
			}
		}
	}
	else
	{
		for (int m = N0; m<N; ++m)
		{
			int j = P[m];
			int J = LM[j];
			int n = 0;
			double* pm = pd + pointers[J] - offset;
			int* pi = indices + pointers[J] - offset;
			int l = pointers[J + 1] - pointers[J];
			for (int k = N0; k<N; ++k)
			{
				int i = P[k];
				int I = LM[i] + offset;
				for (; n<l; ++n)
				if (pi[n] == I)
				{
					pm[n] += ke[i][j];
					break;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void CompactUnSymmMatrix::Assemble(matrix& ke, vector<int>& LMi, vector<int>& LMj)
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
				if ((J = LMj[j]) >= 0) add(I,J, ke[i][j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
/*
void CompactUnSymmMatrix::add(int i, int j, double v)
{
	if (m_brow_based)
	{
		i ^= j; j ^= i; i ^= j;
	}
	assert((i>=0) && (i<m_ndim));
	assert((j>=0) && (j<m_ndim));

	int* pi = m_pindices + m_ppointers[j];
	pi -= m_offset;
	i += m_offset;
	double* pd = m_pd + (m_ppointers[j] - m_offset);
	int l = m_ppointers[j+1] - m_ppointers[j];
	for (int n=0; n<l; ++n, ++pd, ++pi)
	{
		if (*pi == i)
		{
			*pd += v;
			return;
		}
	}
	assert(false);
}
*/

//-----------------------------------------------------------------------------
// This algorithm uses a binary search for locating the correct row index
// This assumes that the indices are ordered!
void CompactUnSymmMatrix::add(int i, int j, double v)
{
	if (m_brow_based)
	{
		i ^= j; j ^= i; i ^= j;
	}
	assert((i >= 0) && (i<m_ndim));
	assert((j >= 0) && (j<m_ndim));

	int* pi = m_pindices + m_ppointers[j];
	pi -= m_offset;
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
void CompactUnSymmMatrix::set(int i, int j, double v)
{
	if (m_brow_based)
	{
		i ^= j; j ^= i; i ^= j;
	}
	int* pi = m_pindices + m_ppointers[j];
	pi -= m_offset;
	int l = m_ppointers[j+1] - m_ppointers[j];
	for (int n=0; n<l; ++n)
	{
		if (pi[n] == i + m_offset)
		{
			m_pd[ m_ppointers[j] + n - m_offset ] = v;
			return;
		}
	}
	assert(false);
}

//-----------------------------------------------------------------------------
double CompactUnSymmMatrix::get(int i, int j)
{
	if (m_brow_based)
	{
		i ^= j; j ^= i; i ^= j;
	}
	int* pi = m_pindices + m_ppointers[j];
	pi -= m_offset;
	int l = m_ppointers[j+1] - m_ppointers[j];
	for (int n=0; n<l; ++n)
	{
		if (pi[n] == i + m_offset) return m_pd[ m_ppointers[j] + n - m_offset ];
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool CompactUnSymmMatrix::check(int i, int j)
{
	if (m_brow_based)
	{
		i ^= j; j ^= i; i ^= j;
	}
	int* pi = m_pindices + m_ppointers[j];
	pi -= m_offset;
	int l = m_ppointers[j + 1] - m_ppointers[j];
	for (int n = 0; n<l; ++n)
	{
		if (pi[n] == i + m_offset) return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
double CompactUnSymmMatrix::diag(int i)
{
	int* pi = m_pindices + m_ppointers[i] - m_offset;
	int l = m_ppointers[i+1] - m_ppointers[i];
	for (int n=0; n<l; ++n)
	{
		if (pi[n] == i + m_offset)
		{
			return m_pd[ m_ppointers[i] + n - m_offset ];
		}
	}

	assert(false);

	return 0;
}

//-----------------------------------------------------------------------------
void CompactUnSymmMatrix::mult_vector(const vector<double>& x, vector<double>& r)
{
	// get the matrix size
	const int N = Size();

	// loop over all columns
	if (m_brow_based)
	{
		// loop over all rows
		for (int i=0; i<N; ++i)
		{
			double ri = 0.0;
			double* pv = m_pd + m_ppointers[i] - m_offset;
			int* pi = m_pindices + m_ppointers[i] - m_offset;
			int n = m_ppointers[i+1] - m_ppointers[i];
			for (int j=0; j<n; ++j) ri += pv[j]*x[pi[j]-m_offset];
			r[i] = ri;
		}
	}
	else
	{
		// loop over all columns
		for (int i=0; i<N; ++i)
		{
			double* pv = m_pd  + m_ppointers[i] - m_offset;
			int* pi = m_pindices + m_ppointers[i]  - m_offset;
			int n = m_ppointers[i+1] - m_ppointers[i];
			for (int j=1; j<n; j++)  r[pi[j] - m_offset] += pv[j]*x[i];
		}
	}
}
