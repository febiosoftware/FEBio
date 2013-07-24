#pragma once
#include <FECore\SparseMatrix.h>

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format.

//! This is the base class for the symmetric and unsymmetric classes

class CompactMatrix : public SparseMatrix
{
public:
	//! constructor
	CompactMatrix( int offset );

	//! destructor
	virtual ~CompactMatrix();

	//! Clear 
	void Clear();

public:
	//! Create the matrix
	virtual void Create(int N, int nz, double* pv, int *pi, int* pp);

	//! add an item to the matrix
	virtual void add(int i, int j, double v) = 0;

	//! set the matrix item
	virtual void set(int i, int j, double v) = 0;

	//! get the matrix item
	virtual double get(int i, int j) { return 0; }

	//! get the diagonal entry
	virtual double diag(int i) = 0;

public:
	//! Pointer to matrix values
	double* Values  () { return m_pd;   }

	//! Pointer to matrix indices
	int*    Indices () { return m_pindices;  }

	//! pointer to matrix row pointers
	int*    Pointers() { return m_ppointers; }

	//! return the index offset (is 0 or 1)
	int     Offset  () { return m_offset; }

protected:
	int*	m_pindices;		//!< indices
	int*	m_ppointers;	//!< pointers
	int		m_offset;		//!< adjust array indices for fortran arrays
};

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format.

//! This class also assumes the matrix is symmetric and therefor only stores
//! the lower triangular matrix

class CompactSymmMatrix : public CompactMatrix
{
public:
	//! class constructor
	CompactSymmMatrix( int offset = 0 );

	//! Create the matrix structure from the SparseMatrixProfile.
	void Create(SparseMatrixProfile& mp);

	//! Allocate storage for matrix data
	void Create(int N, int nz, double* pv, int *pi, int* pp) { CompactMatrix::Create(N, nz, pv, pi, pp); }

	//! Assemble an element matrix into the global matrix
	void Assemble(matrix& ke, vector<int>& lm);

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj);

	//! add a matrix item
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

	//! set matrix item
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

	//! get a matrix item
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

	//! return the diagonal component
	double diag(int i)
	{
		return m_pd[m_ppointers[i] - m_offset];
	}

	//! multiply with vector
	void mult_vector(const vector<double>& x, vector<double>& r);

	//! multiply with vector
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
	//! constructor
	CompactUnSymmMatrix( int offset = 0, bool row_based = false );

	//! Create the matrix structure from the SparseMatrixProfile
	void Create(SparseMatrixProfile& mp);

	//! Assemble the element matrix into the global matrix
	void Assemble(matrix& ke, vector<int>& lm);

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj);

	//! add a value to the matrix item
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

	//! set the matrix item
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

	//! get a matrix item
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

	//! return the diagonal value
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
	bool m_brow_based;	//!< flag indicating whether the matrix is stored row-based on column-based

};

//-----------------------------------------------------------------------------
//! Output Harwell-Boeing compact matrix
void write_hb(CompactMatrix& m, FILE* fp); 

//-----------------------------------------------------------------------------
//! read Symmetric compact matrix data
void read_hb(CompactSymmMatrix& m, FILE* fp);
