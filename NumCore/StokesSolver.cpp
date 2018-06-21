#include "StokesSolver.h"
#include "FGMRESSolver.h"

class SchurCompliment : public SparseMatrix
{
public:
	SchurCompliment(LinearSolver* A, SparseMatrix* B, SparseMatrix* C)
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

//-----------------------------------------------------------------------------
//! constructor
StokesSolver::StokesSolver()
{
	m_pA = 0;
	m_tol = 1e-12;
	m_maxiter = 150;
	m_iter = 0;
	m_printLevel = 0;
}

//-----------------------------------------------------------------------------
//! constructor
StokesSolver::~StokesSolver()
{
}

//-----------------------------------------------------------------------------
void StokesSolver::SetRelativeTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
// get the iteration count
int StokesSolver::GetIterations() const
{
	return m_iter;
}

//-----------------------------------------------------------------------------
// set the print level
void StokesSolver::SetPrintLevel(int n)
{
	m_printLevel = n;
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* StokesSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	return (m_pA = new BlockMatrix());
}

//-----------------------------------------------------------------------------
//! Preprocess 
bool StokesSolver::PreProcess()
{
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// get the number of partitions
	// and make sure we have two
	int NP = m_pA->Partitions();
	if (NP != 2) return false;

	// allocate solvers for diagonal blocks
	m_solver = new PardisoSolver();
	BlockMatrix::BLOCK& Bi = m_pA->Block(0, 0);
	m_solver->SetSparseMatrix(Bi.pA);
	if (m_solver->PreProcess() == false) return false;

	m_iter = 0;

	return true;
}

//-----------------------------------------------------------------------------
//! Factor matrix
bool StokesSolver::Factor()
{
	// factor the diagonal matrix
	return m_solver->Factor();
}

//-----------------------------------------------------------------------------
//! Backsolve the linear system
bool StokesSolver::BackSolve(vector<double>& x, vector<double>& b)
{
	// get the partition sizes
	int n0 = m_pA->PartitionEquations(0);
	int n1 = m_pA->PartitionEquations(1);
	assert(x.size() == (n0+n1));

	// Get the blocks
	BlockMatrix::BLOCK& A = m_pA->Block(0, 0);
	BlockMatrix::BLOCK& B = m_pA->Block(0, 1);
	BlockMatrix::BLOCK& C = m_pA->Block(1, 0);

	// split right hand side in two
	vector<double> F(n0), G(n1);
	for (int i=0; i<n0; ++i) F[i] = b[i];
	for (int i=0; i<n1; ++i) G[i] = b[i + n0];

	// step 1: solve Ay = F
	vector<double> y(n0);
	m_solver->BackSolve(y, F);

	// step 2: calculate H = Cy - G
	vector<double> H(n1);
	C.vmult(y, H);
	H -= G;

	// step 3: Solve Sv = H
	SchurCompliment S(m_solver, B.pA, C.pA);
	vector<double> v(n1);
	FGMRESSolver fgmres;
	fgmres.SetPrintLevel(m_printLevel);
	bool bconv = fgmres.Solve(&S, v, H);

	// step 4: calculate L = F - Bv
	vector<double> tmp(n0);
	B.vmult(v, tmp);
	vector<double> L = F - tmp;

	// step 5: solve Au = L
	vector<double> u(n0, 0.0);
	m_solver->BackSolve(u, L);

	// put it back together
	for (int i=0; i<n0; ++i) x[i   ] = u[i];
	for (int i=0; i<n1; ++i) x[i+n0] = v[i];

	return bconv;
}

//-----------------------------------------------------------------------------
//! Clean up
void StokesSolver::Destroy()
{
	m_solver->Destroy();
}
