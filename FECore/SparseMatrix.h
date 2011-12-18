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
#include <vector>
#include "MatrixProfile.h"
#include "NumCore/matrix.h"
using namespace std;

namespace FECore {

//=============================================================================
//! Base class for sparse matrices

//! This is the base class for the sparse matrix classes and defines the interface
//! to the different matrix classes

class SparseMatrix
{
public:
	SparseMatrix();
	virtual ~SparseMatrix() {}

	//! Create a sparse matrix from a sparse-matrix profile
	virtual void Create(SparseMatrixProfile& MP) = 0;

	//! assemble a matrix into the sparse matrix
	virtual void Assemble(matrix& ke, vector<int>& lm) = 0;

	//! assemble a matrix into the sparse matrix
	virtual void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj) = 0;

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

protected:
	int	m_ndim;	// dimension of matrix

	double*	m_pd;		// matrix values
	int	m_nsize;	// size of m_pd array
};

//-----------------------------------------------------------------------------
void print(SparseMatrix& A, FILE* fp, int i0 = 0, int j0 = 0, int i1 = -1, int j1 = -1);

} // namespace FECore

#endif // !defined(AFX_SPARSEMATRIX_H__B6DFA524_679D_4A35_86F8_D7F080D0ACD5__INCLUDED_)
