#pragma once
#include <FECore/Preconditioner.h>

//-----------------------------------------------------------------------------
class ILU0_Preconditioner : public Preconditioner
{
public:
	ILU0_Preconditioner(FEModel* fem);

	// create a preconditioner for a sparse matrix
	bool Create(SparseMatrix* A) override;

	// apply to vector P x = y
	void mult_vector(double* x, double* y) override;

public:
	bool	m_checkZeroDiagonal;	// check for zero diagonals
	double	m_zeroThreshold;		// threshold for zero diagonal check
	double	m_zeroReplace;			// replacement value for zero diagonal

private:
	vector<double>		m_bilu0;
	vector<double>		m_tmp;
	CRSSparseMatrix*	m_K;
};

