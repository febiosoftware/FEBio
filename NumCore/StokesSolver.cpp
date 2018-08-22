#include "StokesSolver.h"
#include "RCICGSolver.h"
#include "SchurComplement.h"

//-----------------------------------------------------------------------------
//! constructor
StokesSolver::StokesSolver()
{
	m_pA = 0;
	m_tol = 1e-12;
	m_maxiter = 0;
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
// set max nr of iterations
void StokesSolver::SetMaxIterations(int n)
{
	m_maxiter = n;
}

//-----------------------------------------------------------------------------
// set convergence tolerance
void StokesSolver::SetConvergenceTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
//! Set the partition
void StokesSolver::SetPartitions(const vector<int>& part)
{
	m_npart = part;
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* StokesSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (ntype != REAL_SYMMETRIC) return 0;

	if (m_npart.size() != 2) return 0;
	m_pA = new BlockMatrix();
	m_pA->Partition(m_npart, ntype);
	return m_pA;
}

//-----------------------------------------------------------------------------
//! set the sparse matrix
bool StokesSolver::SetSparseMatrix(SparseMatrix* A)
{
	m_pA = dynamic_cast<BlockMatrix*>(A);
	if (m_pA == 0) return false;
	return true;
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
	RCICGSolver* cg = new RCICGSolver();
	cg->SetPreconditioner(new DiagonalPreconditioner);
	cg->SetMaxIterations(m_maxiter);
	cg->SetPrintLevel(m_printLevel);
	m_solver = cg;
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
	if (m_printLevel == 2) fprintf(stdout, "----------------------\nstep 1:\n");
	m_solver->BackSolve(y, F);

	// step 2: calculate H = Cy - G
	vector<double> H(n1);
	C.vmult(y, H);
	H -= G;

	// step 3: Solve Sv = H
	SchurComplement S(m_solver, B.pA, C.pA);
	vector<double> v(n1);
	RCICGSolver cg;
	cg.SetPrintLevel(m_printLevel);
	if (m_maxiter > 0) cg.SetMaxIterations(m_maxiter);
	if (m_printLevel == 2) fprintf(stdout, "step 3:\n");
	bool bconv = cg.Solve(&S, v, H);

	// step 4: calculate L = F - Bv
	vector<double> tmp(n0);
	B.vmult(v, tmp);
	vector<double> L = F - tmp;

	// step 5: solve Au = L
	vector<double> u(n0, 0.0);
	if (m_printLevel == 2) fprintf(stdout, "step 5:\n");
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
