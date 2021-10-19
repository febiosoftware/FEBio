/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "BlockSolver.h"
#include <FECore/log.h>

BEGIN_FECORE_CLASS(BlockIterativeSolver, IterativeLinearSolver)
	ADD_PARAMETER(m_maxiter   , "max_iter");
	ADD_PARAMETER(m_printLevel, "print_level");
	ADD_PARAMETER(m_tol        , "tol");
	ADD_PARAMETER(m_failMaxIter, "fail_max_iter");
	ADD_PARAMETER(m_method     , "solution_method");
	ADD_PARAMETER(m_zeroInitGuess, "zero_initial_guess");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
BlockIterativeSolver::BlockIterativeSolver(FEModel* fem) : IterativeLinearSolver(fem)
{
	m_pA = 0;
	m_tol = 1e-8;
	m_maxiter = 150;
	m_iter = 0;
	m_printLevel = 0;
	m_failMaxIter = true;
	m_method = JACOBI;
	m_zeroInitGuess = true;
}

//-----------------------------------------------------------------------------
//! constructor
BlockIterativeSolver::~BlockIterativeSolver()
{
}

//-----------------------------------------------------------------------------
// return whether the iterative solver has a preconditioner or not
bool BlockIterativeSolver::HasPreconditioner() const
{
	return false;
}

//-----------------------------------------------------------------------------
void BlockIterativeSolver::SetRelativeTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
// set the max nr of iterations
void BlockIterativeSolver::SetMaxIterations(int maxiter)
{
	m_maxiter = maxiter;
}

//-----------------------------------------------------------------------------
// get the iteration count
int BlockIterativeSolver::GetIterations() const
{
	return m_iter;
}

//-----------------------------------------------------------------------------
// set the print level
void BlockIterativeSolver::SetPrintLevel(int n)
{
	m_printLevel = n;
}

//-----------------------------------------------------------------------------
// set the solution method
void BlockIterativeSolver::SetSolutionMethod(int method)
{
	m_method = method;
}

//-----------------------------------------------------------------------------
// set fail on max iterations flag
void BlockIterativeSolver::SetFailMaxIters(bool b)
{
	m_failMaxIter = b;
}

//-----------------------------------------------------------------------------
// set the zero-initial-guess flag
void BlockIterativeSolver::SetZeroInitialGuess(bool b)
{
	m_zeroInitGuess = b;
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* BlockIterativeSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	m_pA = new BlockMatrix();
	m_pA->Partition(m_part, ntype);
	return m_pA;
}

//-----------------------------------------------------------------------------
//! set the sparse matrix
bool BlockIterativeSolver::SetSparseMatrix(SparseMatrix* m)
{
	m_pA = dynamic_cast<BlockMatrix*>(m);
	if (m_pA)
	{
		vector<int> p(m_pA->Partitions());
		for (int i = 0; i < m_pA->Partitions(); ++i) p[i] = m_pA->PartitionEquations(i);
		SetPartitions(p);
	}
	return (m_pA != nullptr);
}

//-----------------------------------------------------------------------------
//! Preprocess 
bool BlockIterativeSolver::PreProcess()
{
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// get the number of partitions
	int NP = m_pA->Partitions();

	// allocate solvers for diagonal blocks
	m_solver.resize(NP);
	for (int i=0; i<NP; ++i)
	{
		m_solver[i] = new PardisoSolver(GetFEModel());
		BlockMatrix::BLOCK& Bi = m_pA->Block(i,i);
		m_solver[i]->SetSparseMatrix(Bi.pA);
		if (m_solver[i]->PreProcess() == false) return false;
	}

	m_iter = 0;

	return true;
}

//-----------------------------------------------------------------------------
//! Factor matrix
bool BlockIterativeSolver::Factor()
{
	// factor the diagonal matrices
	int N = (int) m_solver.size();
	for (int i=0; i<N; ++i) m_solver[i]->Factor();

	return true;
}

//-----------------------------------------------------------------------------
//! Backsolve the linear system
bool BlockIterativeSolver::BackSolve(double* x, double* b)
{
	// get partitions
	int NP = m_pA->Partitions();
	assert(NP == Partitions());

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

		// also allocate and initialize solution vectors
		// we assume that the passed x is an initial guess
		X[i].resize(neq);
		if (m_zeroInitGuess)
		{
			zero(X[i]);
		}
		else
		{
			for (int j = 0; j < neq; ++j) X[i][j] = x[j + neq0];
		}

		neq0 += neq;
	}
	assert(neq0 == m_pA->Rows());

	// residual vector
	vector<double> res(neq0);

	// calculate initial norm
	m_pA->mult_vector(x, &res[0]);
	for (int i = 0; i<neq0; ++i) res[i] -= b[i];
	double norm0 = l2_norm(res);
	if (m_printLevel == 1) feLog("%d: %lg\n", 0, norm0);

	// temp storage for RHS
	vector< vector<double> > T(R.size());

	// solve the linear system iteratively
	bool bconv = false;
	m_iter = 0;
	double norm = 0.0;
	for (int n=0; n<m_maxiter; ++n)
	{
		// loop over rows
		for (int i=0; i<NP; ++i)
		{
			T[i].assign(R[i].size(), 0.0);

			// loop over columns
			for (int j=0; j<NP; ++j)
			{
				if (i != j)
				{
					// get the off-diagonal matrix
					CompactMatrix& Cij = *(m_pA->Block(i,j).pA);

					// multiply with X[j] and add to T[i]
					vector<double>& Xj = X[j];
					vector<double>& Ti = T[i];

					if (l2_sqrnorm(Xj) != 0.0)
						Cij.mult_vector(&Xj[0], &Ti[0]);
				}
			}

			// subtract temp from RHS
			int neq = m_pA->PartitionEquations(i);
			for (int j=0; j<neq; ++j) T[i][j] = R[i][j] - T[i][j];

			if (m_method == GAUSS_SEIDEL)
			{
				if (l2_sqrnorm(T[i]) != 0.0)
				{
					if (m_solver[i]->BackSolve(&X[i][0], &T[i][0]) == false)
						return false;
				}
				else zero(X[i]);
			}
		}

		// backsolve the equations
		if (m_method == JACOBI)
		{
			for (int i = 0; i < NP; ++i)
			{
				if (m_solver[i]->BackSolve(&X[i][0], &T[i][0]) == false)
					return false;
			}
		}

		// combine solution into single solution vector
		neq0 = 0;
		for (int i = 0; i<NP; ++i)
		{
			int neq = m_pA->PartitionEquations(i);
			vector<double>& Xi = X[i];
			for (int j = 0; j<neq; ++j) x[neq0 + j] = Xi[j];
			neq0 += neq;
		}

		// increment iteration count
		m_iter++;

		// calculate residual
		m_pA->mult_vector(&x[0], &res[0]);
		for (int i=0; i<neq0; ++i) res[i] -= b[i];
		norm = l2_norm(res);
		if (m_printLevel == 1) feLog("%d: %lg\n", m_iter, norm);
		if (norm <= norm0*m_tol)
		{
			bconv = true;
			break;	
		}
	}

	if (m_printLevel == 2) feLog("%d: %lg\n", m_iter, norm);

	if ((bconv == false) && (m_failMaxIter == false))
	{
		if (m_printLevel != 0)
			feLogWarning("Iterative solver reached max iterations, but convergence is forced.");
		bconv = true;
	}

	UpdateStats(m_iter);

	return bconv;
}

//-----------------------------------------------------------------------------
//! Clean up
void BlockIterativeSolver::Destroy()
{
	int N = (int) m_solver.size();
	for (int i=0; i<N; ++i) m_solver[i]->Destroy();
}
