#pragma once
#include <FECore/SparseMatrix.h>

class CRSSparseMatrix;

// Base class for preconditioners for iterative linear solvers
class Preconditioner
{
public:
	Preconditioner();
	virtual ~Preconditioner();

	// create a preconditioner for a sparse matrix
	virtual bool Create(SparseMatrix* A) = 0;
	
	// apply to vector P x = y
	virtual void mult_vector(double* x, double* y) = 0;
};

class ILU0_Preconditioner : public Preconditioner
{
public:
	ILU0_Preconditioner();

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

class ILUT_Preconditioner : public Preconditioner
{
public:
	ILUT_Preconditioner();

	// create a preconditioner for a sparse matrix
	bool Create(SparseMatrix* A) override;

	// apply to vector P x = y
	void mult_vector(double* x, double* y) override;

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
