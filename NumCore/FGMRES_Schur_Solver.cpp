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
#include "FGMRES_Schur_Solver.h"

//=============================================================================
FGMRES_Schur_Solver::FGMRES_Schur_Solver(FEModel* fem) : FGMRESSolver(fem)
{
	SetPreconditioner(m_PC = new SchurPreconditioner(fem));
}

//-----------------------------------------------------------------------------
void FGMRES_Schur_Solver::ZeroDBlock(bool b)
{
	m_PC->ZeroDBlock(b);
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_Schur_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (m_part.size() != 2) return 0;
	BlockMatrix* pA = new BlockMatrix();

	int offset = 1;
	int m_nsolver = m_PC->GetLinearSolver();
	if ((m_nsolver == SchurSolver::A_Solver_HYPRE) ||
		(m_nsolver == SchurSolver::A_Solver_FGMRES_AMG)) offset = 0;

	pA->Partition(m_part, ntype, offset);
	FGMRESSolver::SetSparseMatrix(pA);
	return pA;
}

//-----------------------------------------------------------------------------
// this is used to build the preconditioner
bool FGMRES_Schur_Solver::Factor()
{
	if (FGMRESSolver::Factor() == false) return false;
//	m_PC->SetMaxIterations(GetMaxIterations());
	m_PC->SetSparseMatrix(GetSparseMatrix());
	return m_PC->Create();
}
