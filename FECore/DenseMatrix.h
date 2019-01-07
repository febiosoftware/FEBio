#pragma once
#include "SparseMatrix.h"

//=============================================================================
//! This class implements a full matrix
//! that is a matrix that stores all its elements.

class FECORE_API DenseMatrix : public SparseMatrix
{
public:
	// con/de-structor
	DenseMatrix();
	~DenseMatrix();

	// create a matrix of particular size
	void Create(int rows, int cols);

	// retrieve matrix data
	double& operator () (int i, int j) { return m_pr[i][j]; }

public:
	// zero matrix elements
	void Zero() override;

	// create a matrix from a spares matrix profile
	void Create(SparseMatrixProfile& mp) override;

	// assemble matrix into sparse matrix
	void Assemble(matrix& ke, vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj) override;

	// clear all data
	void Clear() override;

	bool check(int i, int j) override { return true; }
	void add(int i, int j, double v) override { m_pr[i][j] += v; }
	void set(int i, int j, double v) override { m_pr[i][j] = v; }
	double diag(int i) override { return m_pr[i][i]; }

protected:
	double*		m_pd;	//!< matrix values
	double**	m_pr;	//!< pointers to rows
};
