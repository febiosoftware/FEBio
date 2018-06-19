#pragma once
#include "FECore/SparseMatrix.h"

//=============================================================================
//! Implements a sparse matrix using the skyline storage

//! This class implements a symmetric sparse matrix where only the values
//! below the skyline are stored.

class SkylineMatrix : public SparseMatrix
{
public:
	SkylineMatrix();
	virtual ~SkylineMatrix();

public: // from SparseMatrix

	void Zero() override;

	void Clear() override;

	void Create(SparseMatrixProfile& mp) override;

	void Assemble(matrix& ke, vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj) override;

	void add(int i, int j, double v) override;

	void set(int i, int j, double v) override;

	// NOTE: This is not implemented yet!
	bool check(int i, int j) override;

	double get(int i, int j);

	double diag(int i);

	double* values() { return m_pd; }
	int* pointers() { return m_ppointers; }

protected:
	void Create(double* pv, int* pp, int N);

protected:
	double*	m_pd;			//!< matrix values
	int*	m_ppointers;	//!< arrays of indices to diagonal elements
};
