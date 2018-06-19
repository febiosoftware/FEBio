#pragma once
#include "CompactMatrix.h"

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format (i.e. column major, lower triangular compact).

//! This class also assumes the matrix is symmetric and therefor only stores
//! the lower triangular matrix

class CompactSymmMatrix : public CompactMatrix
{
public:
	//! class constructor
	CompactSymmMatrix(int offset = 0);

	//! Create the matrix structure from the SparseMatrixProfile.
	void Create(SparseMatrixProfile& mp) override;

	//! Assemble an element matrix into the global matrix
	void Assemble(matrix& ke, vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj) override;

	//! add a matrix item
	void add(int i, int j, double v) override;

	//! set matrix item
	void set(int i, int j, double v) override;

	//! get a matrix item
	double get(int i, int j) override;

	//! return the diagonal component
	double diag(int i) override { return m_pd[m_ppointers[i] - m_offset]; }

	//! multiply with vector
	void mult_vector(double* x, double* r) override;

	//! see if a matrix element is defined
	bool check(int i, int j) override;

	//! is the matrix symmetric or not
	bool isSymmetric() override { return true; }

	//! this is a column based format
	bool isRowBased() override { return false; }
};
