#pragma once
#include "FECore/SparseMatrix.h"
#include "CSRMatrix.h"

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

	//! zero matrix elements
	void zero();

	//! Clear 
	void Clear();

public:
	//! multiply with vector
	virtual void mult_vector(const vector<double>& x, vector<double>& r) = 0;

public:
	//! Pointer to matrix values
	double* Values  () { return m_pd;   }

	//! Pointer to matrix indices
	int*    Indices () { return m_pindices;  }

	//! pointer to matrix row pointers
	int*    Pointers() { return m_ppointers; }

	//! return the index offset (is 0 or 1)
	int     Offset  () { return m_offset; }

public:
	//! Create the matrix
	void alloc(int N, int nz, double* pv, int *pi, int* pp);

	//! see if a matrix element is defined
	virtual bool check(int i, int j) = 0;

protected:
	double*	m_pd;			//!< matrix values
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

	//! Assemble an element matrix into the global matrix
	void Assemble(matrix& ke, vector<int>& lm);

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj);

	//! add a matrix item
	void add(int i, int j, double v);

	//! set matrix item
	void set(int i, int j, double v);

	//! get a matrix item
	double get(int i, int j);

	//! return the diagonal component
	double diag(int i) { return m_pd[m_ppointers[i] - m_offset]; }

	//! multiply with vector
	void mult_vector(const vector<double>& x, vector<double>& r);

	//! multiply with vector
	void mult_vector(const double* x, double* r);

	//! see if a matrix element is defined
	bool check(int i, int j);
};

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format

//! Unlike CompactMatrix does not assume the matrix is symmetric.
//! This still assumes that the sparsity pattern is symmetric for
//! row based formats.

class CompactUnSymmMatrix : public CompactMatrix
{
public:
	class Block
	{
	public:
		Block(CompactUnSymmMatrix* A, int i, int j, int nr, int nc);

		// multiply this block with a vector x and store in r
		void mult(const vector<double>& x, vector<double>& r);
		void mult(const double* x, double* r);

	private:
		CompactUnSymmMatrix*	m_A;	// the parent matrix
		int	m_i, m_j;					// start row,column indices
		int	m_nr, m_nc;					// row and column size of block
	};

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
	void add(int i, int j, double v);

	//! set the matrix item
	void set(int i, int j, double v);

	//! get a matrix item
	double get(int i, int j);

	//! return the diagonal value
	double diag(int i);

	//! multiply with vector
	void mult_vector(const vector<double>& x, vector<double>& r);

	//! see if a matrix element is defined
	bool check(int i, int j);

	// multiply a block of the matrix with a vector
	// This multiplies the vector x with a block of this matrix.
	// The block starts at index i,j, and contains rows and cols
	void blockmult(int i0, int j0, int rows, int cols, const vector<double>& x, vector<double>& r);
	void blockmult(int i0, int j0, int rows, int cols, const double* x, double* r);

	// scale matrix 
	void scale(const vector<double>& L, const vector<double>& R);

	//! extract a block of this matrix
	void get(int i0, int j0, int nr, int nc, CSRMatrix& M);

protected:
	bool m_brow_based;	//!< flag indicating whether the matrix is stored row-based on column-based
};
