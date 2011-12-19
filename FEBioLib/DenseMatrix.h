#pragma once
#include "NumCore/SparseMatrix.h"
using namespace NumCore;

//=============================================================================
//! This class implements a full matrix

//! that is a matrix that stores all its elements.

class DenseMatrix : public SparseMatrix
{
public:
	// con/de-structor
	DenseMatrix();
	virtual ~DenseMatrix();

	// create a matrix from a spares matrix profile
	void Create(SparseMatrixProfile& mp) { Create(mp.size()); }

	// assemble matrix into sparse matrix
	void Assemble(matrix& ke, vector<int>& lm);

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj);

	// clear all data
	void Clear()
	{
		if (m_pd) delete [] m_pd; m_pd = 0;
		if (m_pr) delete [] m_pr; m_pr = 0;
		m_ndim = 0;
	}

	// create a matrix of size N x N
	void Create(int N);

	// retrieve matrix data
	double& operator () (int i, int j) { return m_pr[i][j]; }

	void add(int i, int j, double v) { m_pr[i][j] += v; }
	void set(int i, int j, double v) { m_pr[i][j] = v;  }

	double diag(int i) { return m_pr[i][i]; }

protected:
	double**	m_pr;	//!< pointers to rows
};
