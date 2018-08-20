#include <FECore/SparseMatrix.h>
#include <FECore/LinearSolver.h>

class SchurComplement : public SparseMatrix
{
public:
	SchurComplement(LinearSolver* A, SparseMatrix* B, SparseMatrix* C)
	{
		m_A = A;
		m_B = B;
		m_C = C;

		int n0 = m_B->Columns();
		int n1 = m_B->Rows();
		assert(n0 == m_C->Rows());
		assert(n1 == m_C->Columns());

		m_tmp1.resize(n1, 0.0);
		m_tmp2.resize(n1, 0.0);

		m_nrow = n0;
		m_ncol = n0;
	}

	//! multiply with vector
	void mult_vector(double* x, double* r) override
	{
		m_B->mult_vector(x, &m_tmp1[0]);
		m_A->BackSolve(m_tmp2, m_tmp1);
		m_C->mult_vector(&m_tmp2[0], r);
	}

private: // we need to override these functions although we don't want to use them
	void Zero() { assert(false); }
	void Create(SparseMatrixProfile& MP) { assert(false); }
	void Assemble(matrix& ke, std::vector<int>& lm) { assert(false); }
	void Assemble(matrix& ke, std::vector<int>& lmi, std::vector<int>& lmj) { assert(false); }
	bool check(int i, int j) { assert(false); return false; }
	void set(int i, int j, double v)  { assert(false); }
	void add(int i, int j, double v)  { assert(false); }
	double diag(int i)  { assert(false); return 0.0; }

private:
	LinearSolver*	m_A;
	SparseMatrix*	m_B;
	SparseMatrix*	m_C;

	vector<double>	m_tmp1, m_tmp2;
};
