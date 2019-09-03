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
#include "SchurPreconditioner.h"

SchurPreconditioner::SchurPreconditioner(FEModel* fem) : Preconditioner(fem)
{
	m_solver = fecore_alloc(SchurSolver, fem);
	m_nsize = 0;
	m_solver->SetRelativeResidualTolerance(1e-5);
	m_solver->SetMaxIterations(150);
	m_solver->FailOnMaxIterations(true);

	m_solver->SetLinearSolver(SchurSolver::A_Solver_LU);
	m_solver->SetSchurSolver(SchurSolver::Schur_Solver_PC);
	m_solver->SetSchurPreconditioner(SchurSolver::Schur_PC_NONE);
	m_solver->DoJacobiPreconditioning(false);
}

SchurPreconditioner::~SchurPreconditioner()
{
	delete m_solver;
}

void SchurPreconditioner::SetMaxIterations(int n)
{
	m_solver->SetMaxIterations(n);
}

void SchurPreconditioner::DoJacobiPreconditioning(bool b)
{
	m_solver->DoJacobiPreconditioning(b);
}

void SchurPreconditioner::FailOnMaxIterations(bool b)
{
	m_solver->FailOnMaxIterations(b);
}

void SchurPreconditioner::SetPrintLevel(int n)
{
	m_solver->SetPrintLevel(n);
}

void SchurPreconditioner::SetTolerance(double tol)
{
	m_solver->SetRelativeResidualTolerance(tol);
}

void SchurPreconditioner::ZeroDBlock(bool b)
{
	m_solver->ZeroDBlock(b);
}

void SchurPreconditioner::SetLinearSolver(int n)
{
	m_solver->SetLinearSolver(n);
}

void SchurPreconditioner::SetSchurSolver(int n)
{
	m_solver->SetSchurSolver(n);
}

void SchurPreconditioner::SetSchurASolver(int n)
{
	m_solver->SetSchurASolver(n);
}

void SchurPreconditioner::SetSchurBlock(int n)
{
	m_solver->SetSchurBlock(n);
}

void SchurPreconditioner::SetSchurPreconditioner(int n)
{
	m_solver->SetSchurPreconditioner(n);
}

int SchurPreconditioner::GetLinearSolver()
{
	return m_solver->GetLinearSolver();
}


bool SchurPreconditioner::Create()
{
	// clean up old solver
	m_solver->Destroy();

	SparseMatrix* A = GetSparseMatrix();
	m_nsize = A->Rows();
	if (m_solver->SetSparseMatrix(A) == false) return false;
	if (m_solver->PreProcess() == false) return false;
	if (m_solver->Factor() == false) return false;
	return true;
}

// apply to vector P x = y
bool SchurPreconditioner::mult_vector(double* x, double* y)
{
	return m_solver->BackSolve(y, x);
}
