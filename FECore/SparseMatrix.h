// SparseMatrix.h: interface for the SparseMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SPARSEMATRIX_H__B6DFA524_679D_4A35_86F8_D7F080D0ACD5__INCLUDED_)
#define AFX_SPARSEMATRIX_H__B6DFA524_679D_4A35_86F8_D7F080D0ACD5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <memory.h>
#include <assert.h>
#include <stdio.h>
#include "vector.h"

//-----------------------------------------------------------------------------
//! Base class for sparse matrices

//! This is the base class for the sparse matrix classes and defines the interface
//! to the different matrix classes

class SparseMatrix
{
public:
	SparseMatrix();
	virtual ~SparseMatrix();

	// set entry to value
	virtual void set(int i, int j, double v) = 0;

	// add value to entry
	virtual void add(int i, int j, double v) = 0;

	// retrieve value
	virtual double get(int i, int j) { return 0; }

	// get the diagonal
	virtual double diag(int i) = 0;

	virtual void Clear()
	{
		if (m_pd) delete [] m_pd; m_pd = 0;
	}

	void zero() { memset(m_pd, 0, m_nsize*sizeof(double)); };

	int NonZeroes() { return m_nsize; };
	int Size() { return m_ndim; }

	void print(FILE* fp, int i0 = 0, int j0 = 0, int i1 = -1, int j1 = -1);

protected:
	int	m_ndim;	// dimension of matrix

	double*	m_pd;		// matrix values
	int	m_nsize;	// size of m_pd array
};

//-----------------------------------------------------------------------------
//! This class implements a full matrix

//! that is a matrix that stores all its elements.

class FullMatrix : public SparseMatrix
{
public:
	FullMatrix();
	virtual ~FullMatrix();

	void Clear()
	{
		if (m_pd) delete [] m_pd; m_pd = 0;
		if (m_pr) delete [] m_pr; m_pr = 0;
		m_ndim = 0;
	}

	void Create(int N);

	double& operator () (int i, int j) { return m_pr[i][j]; }

	void add(int i, int j, double v) { m_pr[i][j] += v; }
	void set(int i, int j, double v) { m_pr[i][j] = v;  }

	double diag(int i) { return m_pr[i][i]; }

protected:
	double**	m_pr;	// pointers to rows
};

//-----------------------------------------------------------------------------
//! Implements a sparse matrix using the skyline storage

//! This class implements a symmetric sparse matrix where only the values
//! below the skyline are stored.

class SkylineMatrix : public SparseMatrix
{
public:
	SkylineMatrix();
	virtual ~SkylineMatrix();

	void Clear()
	{
		if (m_pd) delete [] m_pd; m_pd = 0;
		if (m_ppointers) delete [] m_ppointers; m_ppointers = 0;
	}

	void Create(double* pv, int* pp, int N);

	void add(int i, int j, double v)
	{
#ifdef _DEBUG
		if (j<i) { j ^= i; i ^= j; j ^= i; }
		int l = m_ppointers[j+1] - m_ppointers[j];
		assert(j-i<l);
#endif // _DEBUG

		// only add to the upper triangular part
		if (j>=i) m_pd[ m_ppointers[j] + j-i] += v;
	}

	void set(int i, int j, double v)
	{
		// only add to the upper triangular part
		if (j>=i) m_pd[ m_ppointers[j] + j-i] = v;
	}

	double get(int i, int j)
	{
		if (i>j) { i^= j; j ^= i; i ^= j; }

		int l = m_ppointers[j+1] - m_ppointers[j];
		if (j-i < l) return m_pd[ m_ppointers[j] + j-i];
		return 0;
	}

	double diag(int i)
	{
		return m_pd[ m_ppointers[i] ];
	}


	double* values() { return m_pd; }
	int* pointers() { return m_ppointers; }

protected:
	int*	m_ppointers;	// arrays of indices to diagonal elements
};

//-----------------------------------------------------------------------------
//! This class stores a sparse matrix in Harwell-Boeing format.

//! This is the base class for the symmetric and unsymmetric classes

class CompactMatrix : public SparseMatrix
{
public:
	CompactMatrix( int offset );
	virtual ~CompactMatrix() { Clear(); }

	void Clear()
	{
		if (m_pd) delete [] m_pd; m_pd = 0;
		if (m_pindices) delete [] m_pindices; m_pindices = 0;
		if (m_ppointers) delete [] m_ppointers; m_ppointers = 0;
	}

	void Create(int N, int nz, double* pv, int *pi, int* pp);

	virtual void add(int i, int j, double v) = 0;

	virtual void set(int i, int j, double v) = 0;

	virtual double get(int i, int j) { return 0; }

	virtual double diag(int i) = 0;

public:
	double* values  () { return m_pd;   }
	int*    indices () { return m_pindices;  }
	int*    pointers() { return m_ppointers; }
	int     offset  () { return m_offset; }

	bool print_hb(); // Output Harwell-Boeing compact matrix

protected:
	int*	m_pindices;
	int*	m_ppointers;
	int		m_offset; // adjust array indices for fortran arrays
};

//-----------------------------------------------------------------------------
//! This class stores a sparse matrix in Harwell-Boeing format.

//! This class also assumes the matrix is symmetric and therefor only stores
//! the lower triangular matrix

class CompactSymmMatrix : public CompactMatrix
{
public:
	CompactSymmMatrix( int offset = 0 );

	void add(int i, int j, double v)
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

	void set(int i, int j, double v)
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

	double get(int i, int j)
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

	double diag(int i)
	{
		return m_pd[m_ppointers[i] - m_offset];
	}


	void mult_vector(const vector<double>& x, vector<double>& r);

};

//-----------------------------------------------------------------------------
//! This class stores a sparse matrix in Harwell-Boeing format

//! Unlike CompactMatrix does not assume the matrix is symmetric.
//! This still assumes that the sparsity pattern is symmetric for
//! row based formats.

class CompactUnSymmMatrix : public CompactMatrix
{
public:
	CompactUnSymmMatrix( int offset = 0, bool row_based = false );

	void add(int i, int j, double v)
	{
		if (m_brow_based)
		{
			i ^= j; j ^= i; i ^= j;
		}
		assert((i>=0) && (i<m_ndim));
		assert((j>=0) && (j<m_ndim));

		int* pi = m_pindices + m_ppointers[j];
		pi -= m_offset;
		int l = m_ppointers[j+1] - m_ppointers[j];
		for (int n=0; n<l; ++n)
		{
			if (pi[n] == i + m_offset)
			{
				m_pd[ m_ppointers[j] + n - m_offset ] += v;
				return;
			}
		}
		assert(false);
	}

	void set(int i, int j, double v)
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

	double get(int i, int j)
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

	double diag(int i)
	{
		int* pi = m_pindices + m_ppointers[i];
		int l = m_ppointers[i+1] - m_ppointers[i];
		for (int n=0; n<l; ++n)
		{
			if (pi[n] == i)
			{
				return m_pd[ m_ppointers[i] + n - m_offset ];
			}
		}

		assert(false);

		return 0;
	}

protected:
	bool m_brow_based;

};

#endif // !defined(AFX_SPARSEMATRIX_H__B6DFA524_679D_4A35_86F8_D7F080D0ACD5__INCLUDED_)
