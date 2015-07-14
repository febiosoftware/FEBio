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

	void zero();

	void Clear();

	void Create(SparseMatrixProfile& mp);

	void Assemble(matrix& ke, vector<int>& lm);

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj);

	void Create(double* pv, int* pp, int N);

	void add(int i, int j, double v);

	void set(int i, int j, double v);

	double get(int i, int j);

	double diag(int i);

	double* values() { return m_pd; }
	int* pointers() { return m_ppointers; }

protected:
	double*	m_pd;			//!< matrix values
	int*	m_ppointers;	//!< arrays of indices to diagonal elements
};
