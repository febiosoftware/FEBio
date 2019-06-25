/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "BlockPreconditioner.h"
#include "BlockMatrix.h"

BlockPreconditioner::BlockPreconditioner(FEModel* fem) : Preconditioner(fem)
{
	m_method = 0;
}

BlockPreconditioner::~BlockPreconditioner()
{
	for (size_t i = 0; i < m_solver.size(); ++i) delete m_solver[i];
	m_solver.clear();
}

// set the solution method
void BlockPreconditioner::SetSolutionMethod(int method)
{
	m_method = method;
}

// create a preconditioner for a sparse matrix
bool BlockPreconditioner::Create()
{
	// make sure this is a block solver
	BlockMatrix* A = dynamic_cast<BlockMatrix*>(GetSparseMatrix());
	if (A == nullptr) return false;

	// get the partitions
	int partitions = A->Partitions();

	// allocate solvers
	for (int i = 0; i < partitions; ++i)
	{
		PardisoSolver* solver = new PardisoSolver(GetFEModel());
		m_solver.push_back(solver);

		CompactMatrix* Kii = A->Block(i, i).pA;
		if (solver->SetSparseMatrix(Kii) == false) return false;

		if (solver->PreProcess() == false) return false;
		if (solver->Factor() == false) return false;
	}

	return true;
}

// apply to vector P x = y
bool BlockPreconditioner::mult_vector(double* x, double* y)
{
	// make sure this is a block solver
	BlockMatrix* A = dynamic_cast<BlockMatrix*>(GetSparseMatrix());
	if (A == nullptr) return false;

	if (m_method == 0)
	{
		int neq = 0;
		for (int i = 0; i < A->Partitions(); ++i)
		{
			m_solver[i]->BackSolve(y + neq, x + neq);
			neq += A->PartitionEquations(i);
		}
	}
	else if (m_method == 1)
	{
		int n = 0;
		for (int i = 0; i < A->Partitions(); ++i)
		{
			int neq = A->PartitionEquations(i);
			vector<double> tmp(neq, 0.0), rhs(neq, 0.0);
			int m = 0;
			for (int j = 0; j < i; ++j)
			{
				SparseMatrix* Kij = A->Block(i, j).pA;
				Kij->mult_vector(y + m, &tmp[0]);
				m += Kij->Columns();

				for (int k = 0; k < neq; ++k) rhs[k] += tmp[k];
			}

			for (int k = 0; k < neq; ++k) rhs[k] = x[k + n] - rhs[k];

			m_solver[i]->BackSolve(y + n, &rhs[0]);
			n += neq;
		}
	}
	else
	{
		// partitions
		int NP = A->Partitions();

		// offsets into x,y arrays
		vector<int> off(NP+1, 0);
		for (int i = 1; i <= NP; ++i) off[i] = off[i - 1] + A->PartitionEquations(i - 1);

		// loop over partitions (backwards!)
		for (int i = NP-1; i >= 0; --i)
		{
			int neq = A->PartitionEquations(i);
			vector<double> tmp(neq, 0.0), rhs(neq, 0.0);
			int m = 0;
			for (int j = i+1; j < NP; ++j)
			{
				SparseMatrix* Kij = A->Block(i, j).pA;
				Kij->mult_vector(y + off[j], &tmp[0]);
				m += Kij->Columns();

				for (int k = 0; k < neq; ++k) rhs[k] += tmp[k];
			}

			for (int k = 0; k < neq; ++k) rhs[k] = x[k + off[i]] - rhs[k];

			m_solver[i]->BackSolve(y + off[i], &rhs[0]);
		}
	}

	return true;
}


//=============================================================================
FGMRES_Jacobi_Block_Solver::FGMRES_Jacobi_Block_Solver(FEModel* fem) : FGMRESSolver(fem)
{
	SetPreconditioner(m_PC = new BlockPreconditioner(fem));
}

//-----------------------------------------------------------------------------
// set the solution method
void FGMRES_Jacobi_Block_Solver::SetSolutionMethod(int method)
{
	m_PC->SetSolutionMethod(method);
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_Jacobi_Block_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (m_part.empty()) return 0;
	BlockMatrix* pA = new BlockMatrix();
	pA->Partition(m_part, ntype, 1);
	m_PC->SetSparseMatrix(pA);
	FGMRESSolver::SetSparseMatrix(pA);
	return pA;
}

//-----------------------------------------------------------------------------
//! set the sparse matrix
bool FGMRES_Jacobi_Block_Solver::SetSparseMatrix(SparseMatrix* A)
{
	BlockMatrix* pA = dynamic_cast<BlockMatrix*>(A);
	if (pA == nullptr) return false;

	int np = pA->Partitions();
	if (np == 0) return false;

	m_part.resize(np, 0);
	for (int i = 0; i < np; ++i) m_part[i] = pA->PartitionEquations(i);

	m_PC->SetSparseMatrix(pA);
	return FGMRESSolver::SetSparseMatrix(pA);
}

//-----------------------------------------------------------------------------
// this is used to build the preconditioner
bool FGMRES_Jacobi_Block_Solver::Factor()
{
	if (FGMRESSolver::Factor() == false) return false;
	return m_PC->Create();
}
