#pragma once
#include "CompactMatrix.h"

//=============================================================================
//! This class stores a general, sparse matrix in Compact Row Storage format
class CRSSparseMatrix : public CompactMatrix
{
public:
	//! constructor
	CRSSparseMatrix(int offset = 0);

	//! copy constructor
	CRSSparseMatrix(const CRSSparseMatrix& A);

	//! Create the matrix structure from the SparseMatrixProfile
	void Create(SparseMatrixProfile& mp) override;

	//! Assemble the element matrix into the global matrix
	void Assemble(matrix& ke, vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj) override;

	//! add a value to the matrix item
	void add(int i, int j, double v) override;

	//! set the matrix item
	void set(int i, int j, double v) override;

	//! get a matrix item
	double get(int i, int j) override;

	//! return the diagonal value
	double diag(int i) override;

	//! multiply with vector
	void mult_vector(double* x, double* r) override;

	//! see if a matrix element is defined
	bool check(int i, int j) override;

	// scale matrix 
	void scale(const vector<double>& L, const vector<double>& R);

	//! extract a block of this matrix
	void get(int i0, int j0, int nr, int nc, CSRMatrix& M);

	//! is the matrix symmetric or not
	bool isSymmetric() override { return false; }

	//! is this a row-based format or not
	bool isRowBased() override { return true; }

	//! calculate the inf norm
	double infNorm() const;
};

//=============================================================================
//! This class stores a sparse matrix in Compact Column Storage format

class CCSSparseMatrix : public CompactMatrix
{
public:
	//! constructor
	CCSSparseMatrix(int offset = 0);

	//! copy constructor
	CCSSparseMatrix(const CCSSparseMatrix& A);

	//! Create the matrix structure from the SparseMatrixProfile
	void Create(SparseMatrixProfile& mp) override;

	//! Assemble the element matrix into the global matrix
	void Assemble(matrix& ke, vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj) override;

	//! add a value to the matrix item
	void add(int i, int j, double v) override;

	//! set the matrix item
	void set(int i, int j, double v) override;

	//! get a matrix item
	double get(int i, int j) override;

	//! return the diagonal value
	double diag(int i) override;

	//! multiply with vector
	void mult_vector(double* x, double* r) override;

	//! see if a matrix element is defined
	bool check(int i, int j) override;

	//! is the matrix symmetric or not
	bool isSymmetric() override { return false; }

	//! is this a row-based format or not
	bool isRowBased() override { return false; }
};
