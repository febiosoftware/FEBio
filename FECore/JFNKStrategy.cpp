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
#include "JFNKStrategy.h"
#include "FENewtonSolver.h"
#include "JFNKMatrix.h"
#include "FEException.h"
#include "LinearSolver.h"
#include "log.h"

BEGIN_FECORE_CLASS(JFNKStrategy, FENewtonStrategy)
	ADD_PARAMETER(m_jfnk_eps, "jfnk_eps");
END_FECORE_CLASS();

JFNKStrategy::JFNKStrategy(FEModel* fem) : FENewtonStrategy(fem)
{
//	m_jfnk_eps = 5e-12;
	m_jfnk_eps = 1e-6;
	m_bprecondition = false;
	m_A = nullptr;

	m_plinsolve = nullptr;
	m_neq = 0;
}

//! New initialization method
bool JFNKStrategy::Init()
{
	if (m_pns == nullptr) return false;
	m_plinsolve = m_pns->GetLinearSolver();
	return true;
}

SparseMatrix* JFNKStrategy::CreateSparseMatrix(Matrix_Type mtype)
{
	if (m_A) delete m_A;
	m_A = 0;

	// make sure the linear solver is an iterative linear solver
	IterativeLinearSolver* ls = dynamic_cast<IterativeLinearSolver*>(m_pns->m_plinsolve);
	if (ls)
	{
		// see if this solver has a preconditioner
		m_bprecondition = ls->HasPreconditioner();

		// if the solver has a preconditioner, we still need to create the stiffness matrix
		SparseMatrix* K = 0;
		if (m_bprecondition) 
		{
			K = ls->CreateSparseMatrix(mtype);
			if (K == 0) return 0;
		}

		// Now, override the matrix used
		m_A = new JFNKMatrix(m_pns, K);
		ls->SetSparseMatrix(m_A);

		m_A->SetEpsilon(m_jfnk_eps);

		// If there is no preconditioner we can do the pre-processing here
		if (m_bprecondition == false)
		{
			ls->PreProcess();
		}
		else
		{
			if (ls->GetLeftPreconditioner()) ls->GetLeftPreconditioner()->SetSparseMatrix(K);
			if (ls->GetRightPreconditioner()) ls->GetRightPreconditioner()->SetSparseMatrix(K);
		}
	}

	return m_A;
}

//! perform a quasi-Newton udpate
bool JFNKStrategy::Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
	// nothing to do here
	return true;
}

//! solve the equations
void JFNKStrategy::SolveEquations(vector<double>& x, vector<double>& b)
{
	// perform a backsubstitution
	if (m_plinsolve->BackSolve(x, b) == false)
	{
		throw LinearSolverFailed();
	}
}

bool JFNKStrategy::ReformStiffness()
{
	if (m_bprecondition) return m_pns->ReformStiffness();
	else return true;
}

//! override so we can store a copy of the residual before we add Fd
bool JFNKStrategy::Residual(std::vector<double>& R, bool binit)
{
	// first calculate the residual
	bool b = m_pns->Residual(R);
	if (b == false) return false;

	// store a copy
	m_A->SetReferenceResidual(R);

	// at the first iteration we need to calculate the Fd
	if (binit)
	{
		// get the vector of prescribed displacements
		std::vector<double>& ui = m_pns->m_ui;

		// build m_Fd
		vector<double>& Fd = m_pns->m_Fd;
		m_A->SetPolicy(JFNKMatrix::ZERO_FREE_DOFS);
		m_A->mult_vector(&ui[0], &Fd[0]);
		m_A->SetPolicy(JFNKMatrix::ZERO_PRESCRIBED_DOFS);

		// we need to flip the sign on Fd
		for (size_t i = 0; i < Fd.size(); ++i) Fd[i] = -Fd[i];

		// NOTE: Usually, the residual is called after an Update. However, here, 
		//       the model is updated (via Update2) in m_A->mult_vector. I wonder
		//       if I need to restore the model's state to before the mult_vector call
		//       by calling another Update2 method with a zero vector.
	}

	return true;
}
