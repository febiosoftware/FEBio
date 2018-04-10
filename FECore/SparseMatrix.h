// SparseMatrix.h: interface for the SparseMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SPARSEMATRIX_H__B6DFA524_679D_4A35_86F8_D7F080D0ACD5__INCLUDED_)
#define AFX_SPARSEMATRIX_H__B6DFA524_679D_4A35_86F8_D7F080D0ACD5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "MatrixProfile.h"
#include "matrix.h"
#include <vector>

//-----------------------------------------------------------------------------
//! Base class for sparse matrices

//! This is the base class for the sparse matrix classes and defines the interface
//! to the different matrix classes

class SparseMatrix
{
public:
	SparseMatrix();
	virtual ~SparseMatrix() {}

public:
	//! return number of nonzeros
	int NonZeroes() { return m_nsize; };

	//! return size, i.e. number of rows (or columns)
	int Size() { return m_ndim; }

public: // functions to be overwritten in derived classes

	//! set all matrix elements to zero
	virtual void zero() = 0;

	//! Create a sparse matrix from a sparse-matrix profile
	virtual void Create(SparseMatrixProfile& MP) = 0;

	//! assemble a matrix into the sparse matrix
	virtual void Assemble(matrix& ke, std::vector<int>& lm) = 0;

	//! assemble a matrix into the sparse matrix
	virtual void Assemble(matrix& ke, std::vector<int>& lmi, std::vector<int>& lmj) = 0;

	//! set entry to value
	virtual void set(int i, int j, double v) = 0;

	//! add value to entry
	virtual void add(int i, int j, double v) = 0;

	//! retrieve value
	virtual double get(int i, int j) { return 0; }

	//! get the diagonal value
	virtual double diag(int i) = 0;

	//! release memory for storing data
	virtual void Clear() = 0;

	//! multiply with vector
	virtual void mult_vector(double* x, double* r) { assert(false); }

protected:
	int	m_ndim;		//!< dimension of matrix
	int	m_nsize;	//!< size of m_pd array
};

#endif // !defined(AFX_SPARSEMATRIX_H__B6DFA524_679D_4A35_86F8_D7F080D0ACD5__INCLUDED_)
