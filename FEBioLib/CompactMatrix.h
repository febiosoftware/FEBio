#pragma once
#include "FECore/SparseMatrix.h"
using namespace FECore;

//=============================================================================
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

	virtual void Create(int N, int nz, double* pv, int *pi, int* pp);

	virtual void add(int i, int j, double v) = 0;

	virtual void set(int i, int j, double v) = 0;

	virtual double get(int i, int j) { return 0; }

	virtual double diag(int i) = 0;

public:
	double* Values  () { return m_pd;   }
	int*    Indices () { return m_pindices;  }
	int*    Pointers() { return m_ppointers; }
	int     Offset  () { return m_offset; }

protected:
	int*	m_pindices;
	int*	m_ppointers;
	int		m_offset; // adjust array indices for fortran arrays
};

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format.

//! This class also assumes the matrix is symmetric and therefor only stores
//! the lower triangular matrix

class CompactSymmMatrix : public CompactMatrix
{
public:
	CompactSymmMatrix( int offset = 0 );

	void Create(SparseMatrixProfile& mp);

	void Create(int N, int nz, double* pv, int *pi, int* pp) { CompactMatrix::Create(N, nz, pv, pi, pp); }

	void Assemble(matrix& ke, vector<int>& lm);

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj);

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
	void mult_vector(const double* x, double* r);
};

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format

//! Unlike CompactMatrix does not assume the matrix is symmetric.
//! This still assumes that the sparsity pattern is symmetric for
//! row based formats.

class CompactUnSymmMatrix : public CompactMatrix
{
public:
	CompactUnSymmMatrix( int offset = 0, bool row_based = false );

	void Create(SparseMatrixProfile& mp);

	void Assemble(matrix& ke, vector<int>& lm);

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj);

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

protected:
	bool m_brow_based;

};

//-----------------------------------------------------------------------------
//! Output Harwell-Boeing compact matrix
void write_hb(CompactMatrix& m, FILE* fp); 

//-----------------------------------------------------------------------------
//! read Symmetric compact matrix data
void read_hb(CompactSymmMatrix& m, FILE* fp);
