#pragma once
#include <FECore/Preconditioner.h>

//-----------------------------------------------------------------------------
class ILUT_Preconditioner : public Preconditioner
{
public:
	ILUT_Preconditioner(FEModel* fem);

	// create a preconditioner for a sparse matrix
	bool Create(SparseMatrix* A) override;

	// apply to vector P x = y
	bool mult_vector(double* x, double* y) override;

public:
	int		m_maxfill;
	double	m_fillTol;
	bool	m_checkZeroDiagonal;	// check for zero diagonals
	double	m_zeroThreshold;		// threshold for zero diagonal check
	double	m_zeroReplace;			// replacement value for zero diagonal

private:
	CRSSparseMatrix*	m_K;
	vector<double>	m_bilut;
	vector<int>		m_jbilut;
	vector<int>		m_ibilut;
	vector<double>	m_tmp;
};
