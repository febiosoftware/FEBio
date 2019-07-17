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
#include "FGMRES_AMG_Solver.h"
#include "CompactUnSymmMatrix.h"

BoomerAMGPC::BoomerAMGPC(FEModel* fem) : Preconditioner(fem)
{
	m_amg = nullptr;
}

// create a preconditioner for a sparse matrix
bool BoomerAMGPC::Create()
{
	if (m_amg) { delete m_amg; m_amg = nullptr; }
	m_amg = fecore_alloc(BoomerAMGSolver, GetFEModel());
	if (m_amg == nullptr) return false;

	// make sure Jacobi PC is off
	if (m_amg->GetJacobiPC())
	{
		m_amg->SetJacobiPC(false);
		fprintf(stderr, "WARNING: Cannot use Jacobi PC in BoomerAMG when this solver is used as a preconditioner.\n");
		fprintf(stderr, "         This option was turned off.\n\n");
	}

	if (m_amg->SetSparseMatrix(GetSparseMatrix()) == false) return false;
	if (m_amg->PreProcess() == false) return false;
	if (m_amg->Factor() == false) return false;

	return true;
}

// apply to vector P x = y
bool BoomerAMGPC::mult_vector(double* x, double* y)
{
	assert(m_amg);
	m_amg->BackSolve(y, x);

	return true;
}

//! constructor
FGMRES_AMG_Solver::FGMRES_AMG_Solver(FEModel* fem) : FGMRESSolver(fem)
{
	m_PC = new BoomerAMGPC(fem);
	SetPreconditioner(m_PC);
}

//! Return a sparse matrix compatible with this solver
SparseMatrix* FGMRES_AMG_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	// allocate the correct matrix format depending on matrix symmetry type
	SparseMatrix* A = nullptr;
	switch (ntype)
	{
	case REAL_UNSYMMETRIC:
	case REAL_SYMM_STRUCTURE:
		A = new CRSSparseMatrix(0); break;
	default:
		return nullptr;
	}

	m_PC->SetSparseMatrix(A);
	FGMRESSolver::SetSparseMatrix(A);

	return A;
}

// this is used to build the preconditioner
bool FGMRES_AMG_Solver::Factor()
{
	if (FGMRESSolver::Factor() == false) return false;
	return m_PC->Create();
}
