#pragma once
#include <FECore/Preconditioner.h>

class CompactSymmMatrix;

class IncompleteCholesky : public Preconditioner
{
public:
	IncompleteCholesky(FEModel* fem);

	// create a preconditioner for a sparse matrix
	bool Create(SparseMatrix* A) override;

	// apply to vector P x = y
	void mult_vector(double* x, double* y) override;

public:
	CompactSymmMatrix* getMatrix();

private:
	CompactSymmMatrix*	m_L;
	vector<double>		z;
};
