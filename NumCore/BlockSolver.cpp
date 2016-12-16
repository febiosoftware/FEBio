#include "BlockSolver.h"

//-----------------------------------------------------------------------------
//! constructor
BlockSolver::BlockSolver()
{
	m_pA = 0;
}

//-----------------------------------------------------------------------------
//! constructor
BlockSolver::~BlockSolver()
{
}

//-----------------------------------------------------------------------------
//! Preprocess 
bool BlockSolver::PreProcess()
{
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// get the number of partitions
	int NP = m_pA->Partitions();

	// allocate solvers for diagonal blocks
	m_solver.resize(NP);
	for (int i=0; i<NP; ++i)
	{
		m_solver[i] = new PardisoSolver();
		BlockMatrix::BLOCK& Bi = m_pA->Block(i,i);
		m_solver[i]->SetSparseMatrix(Bi.pA);
		if (m_solver[i]->PreProcess() == false) return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Factor matrix
bool BlockSolver::Factor()
{
	// factor the diagonal matrices
	int N = (int) m_solver.size();
	for (int i=0; i<N; ++i) m_solver[i]->Factor();

	return true;
}

//-----------------------------------------------------------------------------
//! Backsolve the linear system
bool BlockSolver::BackSolve(vector<double>& x, vector<double>& b)
{
	// get partitions
	int NP = m_pA->Partitions();

	// split right-hand-side and solution vector in partitions
	vector< vector<double> > R(NP);
	vector< vector<double> > X(NP);
	int neq0 = 0;
	for (int i=0; i<NP; ++i)
	{
		int neq = m_pA->PartitionEquations(i);
		vector<double>& Ri = R[i];
		Ri.resize(neq);
		for (int j=0; j<neq; ++j) Ri[j] = b[j + neq0];
		neq0 += neq;

		// also allocate and initialize solution vectors
		X[i].resize(neq);
		zero(X[i]);
	}

	// temp storage for RHS
	vector< vector<double> > T = R;

	// solve the linear system iteratively
	const int MAX_ITERS = 5;
	for (int n=0; n<MAX_ITERS; ++n)
	{
		// loop over rows
		for (int i=0; i<NP; ++i)
		{
			// loop over columns
			for (int j=0; j<NP; ++j)
			{
				if (i != j)
				{
					// get the off-diagonal matrix
					CompactMatrix& Cij = *(m_pA->Block(i,j).pA);

					// multiply with X[j] and add to T[i]
					Cij.mult_vector(X[j], T[i]);
				}
			}

			// subtract temp from RHS
			int neq = m_pA->PartitionEquations(i);
			for (int j=0; j<neq; ++j) T[j] = R[j] - T[j];
		}

		// backsolve the equations
		for (int i=0; i<NP; ++i)
		{
			if (m_solver[i]->BackSolve(X[i], T[i]) == false)
				return false;
		}

		// TODO: check convergence
	}

	// combine solution into single solution vector
	neq0 = 0;
	for (int i=0; i<NP; ++i)
	{
		int neq = m_pA->PartitionEquations(i);
		vector<double>& Xi = X[i];
		for (int j=0; j<neq; ++j) x[neq0 + j] = Xi[j];
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Clean up
void BlockSolver::Destroy()
{
	int N = (int) m_solver.size();
	for (int i=0; i<N; ++i) m_solver[i]->Destroy();
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* BlockSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	return (m_pA = new BlockMatrix());
}
