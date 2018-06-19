#pragma once
#include "SparseMatrix.h"

class FENewtonSolver;

//-----------------------------------------------------------------------------
// This is a class that just mimics a sparse matrix.
// It is only used by the JFNK strategy. 
// The only function it implements is the mult_vector.
class FECORE_API JFNKMatrix : public SparseMatrix
{
public:
	JFNKMatrix(FENewtonSolver* pns, SparseMatrix* K = 0);

	//! override multiply with vector (Does not use sparse matrix m_K)
	void mult_vector(double* x, double* r) override;

public: // these functions use the actual sparse matrix m_K

	//! set all matrix elements to zero
	void Zero() override { m_K->Zero(); }

	//! Create a sparse matrix from a sparse-matrix profile
	void Create(SparseMatrixProfile& MP) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, std::vector<int>& lm) override { m_K->Assemble(ke, lm); }

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, std::vector<int>& lmi, std::vector<int>& lmj) override { m_K->Assemble(ke, lmi, lmj); }

	//! check if an entry was allocated
	bool check(int i, int j) override { return m_K->check(i, j); }

	//! set entry to value
	void set(int i, int j, double v) override { m_K->set(i, j, v); }

	//! add value to entry
	void add(int i, int j, double v) override { m_K->add(i, j, v); }

	//! retrieve value
	double get(int i, int j) override { return m_K->get(i, j); }

	//! get the diagonal value
	double diag(int i) override { return m_K->diag(i); }

	//! release memory for storing data
	void Clear() override { m_K->Clear(); }

	// interface to compact matrices
	double* Values() { return m_K->Values(); }
	int*    Indices() { return m_K->Indices(); }
	int*    Pointers() { return m_K->Pointers(); }
	int     Offset() const { return m_K->Offset(); }

private:
	SparseMatrix*	m_K;		// the actual sparse matrix (This is only used as a preconditioner and can be null)
	FENewtonSolver*	m_pns;
	vector<double>	m_v, m_R;
};
