#include <FECore/SparseMatrix.h>
#include <FECore/LinearSolver.h>

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

class SchurComplement : public SparseMatrix
{
public:
	SchurComplement(LinearSolver* A, SparseMatrix* B, SparseMatrix* C, SparseMatrix* D = 0);

	//! multiply with vector
	void mult_vector(double* x, double* r) override;

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
	SparseMatrix*	m_D;

	vector<double>	m_tmp1, m_tmp2;
};
