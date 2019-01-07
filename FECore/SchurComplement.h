#include "SparseMatrix.h"
#include "LinearSolver.h"

// This class implements a sparse matrix operator that represents the Schur complement of a matrix M
//
//       | A | B |
//  M =  | --+-- |
//       | C | D |
//
// The Schur complement of A, is given by 
//       
//  S\A = C*A^-1*B - D
//
// The only operator this class implements is the matrix-vector multiplication operator (in mult_vector function). 
// This function evaluates a product of the form S*x, where v is a vector, as follows:
// 1. evaluate u = B*x
// 2. evaluate v = A^-1*u, by solving A*v = u
// 3. evaluate r = C*v
// 4. if D is given, then subtract r <-- r - D*x
//
// If D = 0, it does not need to be specified, in which case step 4 is not done.

class FECORE_API SchurComplement : public SparseMatrix
{
public:
	SchurComplement(LinearSolver* A, SparseMatrix* B, SparseMatrix* C, SparseMatrix* D = 0);

	// set the print level
	void SetPrintLevel(int printLevel);

	//! multiply with vector
	bool mult_vector(double* x, double* r) override;

private: // we need to override these functions although we don't want to use them
	void Zero() override { assert(false); }
	void Create(SparseMatrixProfile& MP) override { assert(false); }
	void Assemble(matrix& ke, std::vector<int>& lm) override { assert(false); }
	void Assemble(matrix& ke, std::vector<int>& lmi, std::vector<int>& lmj) override { assert(false); }
	bool check(int i, int j) override { assert(false); return false; }
	void set(int i, int j, double v) override  { assert(false); }
	void add(int i, int j, double v) override  { assert(false); }
	double diag(int i) override  { assert(false); return 0.0; }

private:
	int	m_print_level;
	
	LinearSolver*	m_A;
	SparseMatrix*	m_B;
	SparseMatrix*	m_C;
	SparseMatrix*	m_D;

	vector<double>	m_tmp1, m_tmp2, m_tmp3;
};

//       | A | B |
//  M =  | --+-- |
//       | C | D |
//
// The Schur complement of A, is given by 
//       
//  S\D = B*D^-1*C - A

class SchurComplement2 : public SparseMatrix
{
public:
	SchurComplement2(SparseMatrix* A, SparseMatrix* B, SparseMatrix* C, LinearSolver* D);

	// set the print level
	void SetPrintLevel(int printLevel);

	//! multiply with vector
	bool mult_vector(double* x, double* r) override;

private: // we need to override these functions although we don't want to use them
	void Zero() override { assert(false); }
	void Create(SparseMatrixProfile& MP) override { assert(false); }
	void Assemble(matrix& ke, std::vector<int>& lm) override { assert(false); }
	void Assemble(matrix& ke, std::vector<int>& lmi, std::vector<int>& lmj) override { assert(false); }
	bool check(int i, int j) override { assert(false); return false; }
	void set(int i, int j, double v) override { assert(false); }
	void add(int i, int j, double v) override { assert(false); }
	double diag(int i) override { assert(false); return 0.0; }

private:
	int	m_print_level;

	LinearSolver*	m_D;
	SparseMatrix*	m_B;
	SparseMatrix*	m_C;
	SparseMatrix*	m_A;

	vector<double>	m_tmp1, m_tmp2, m_tmp3;
};
