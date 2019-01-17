#pragma once
#include "SparseMatrix.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class CRSSparseMatrix;
class CompactMatrix;

//-----------------------------------------------------------------------------
// Base class for preconditioners for iterative linear solvers
class FECORE_API Preconditioner : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	Preconditioner(FEModel* fem);
	virtual ~Preconditioner();

	// return the sparse matrix
	SparseMatrix* GetSparseMatrix();

	// set the sparse matrix
	void SetSparseMatrix(SparseMatrix* A);

	// create a preconditioner for a sparse matrix
	virtual bool Create() = 0;
	
	// apply to vector P x = y
	virtual bool mult_vector(double* x, double* y) = 0;

private:
	SparseMatrix*	m_K;
};

//-----------------------------------------------------------------------------
class FECORE_API DiagonalPreconditioner : public Preconditioner
{
public:
	DiagonalPreconditioner(FEModel* fem);

	// create a preconditioner for a sparse matrix
	bool Create() override;

	// apply to vector P x = y
	bool mult_vector(double* x, double* y) override;

private:
	vector<double>	m_D;
};
