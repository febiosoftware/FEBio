#pragma once
#include "SparseMatrix.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class CRSSparseMatrix;
class CompactMatrix;

//-----------------------------------------------------------------------------
// Base class for preconditioners for iterative linear solvers
class Preconditioner : public FECoreBase
{
	DECLARE_SUPER_CLASS(FEPRECONDITIONER_ID);

public:
	Preconditioner(FEModel* fem);
	virtual ~Preconditioner();

	// create a preconditioner for a sparse matrix
	virtual bool Create(SparseMatrix* A) = 0;
	
	// apply to vector P x = y
	virtual void mult_vector(double* x, double* y) = 0;
};

//-----------------------------------------------------------------------------
class DiagonalPreconditioner : public Preconditioner
{
public:
	DiagonalPreconditioner(FEModel* fem);

	// create a preconditioner for a sparse matrix
	bool Create(SparseMatrix* A) override;

	// apply to vector P x = y
	void mult_vector(double* x, double* y) override;

private:
	vector<double>	m_D;
};
