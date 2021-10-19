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
#include "StrategySolver.h"
#include "MatrixTools.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(StrategySolver, LinearSolver)
	ADD_PARAMETER(m_strategy, "strategy");
	ADD_PARAMETER(m_ctol, "max_condition_number");
	ADD_PARAMETER(m_print_cn, "print_condition_number");
	ADD_PROPERTY(m_solver1, "solver1");
	ADD_PROPERTY(m_solver2, "solver2");
END_FECORE_CLASS();

StrategySolver::StrategySolver(FEModel* fem) : LinearSolver(fem)
{
	m_strategy = DONT_PERSIST;
	m_print_cn = false;
	m_ctol = 0;

	m_solver1 = nullptr;
	m_solver2 = nullptr;
	m_activeSolver = nullptr;

	m_pA = nullptr;
}

StrategySolver::~StrategySolver()
{
	Destroy();
	delete m_solver1;
	delete m_solver2;
}

//! Preprocess 
bool StrategySolver::PreProcess()
{
	if (m_solver1->PreProcess() == false) return false;
	if (m_solver2->PreProcess() == false) return false;
	return true;
}

//! Factor matrix
bool StrategySolver::Factor()
{
	if ((m_activeSolver == nullptr) || (m_strategy != PERSIST_FOREVER)) m_activeSolver = m_solver1;

	if (m_pA && (m_print_cn || ((m_ctol > 0.0) && (m_activeSolver == m_solver1))))
	{
		double c = NumCore::estimateConditionNumber(m_pA);

		if (m_print_cn) feLog("\nEst. condition number: %lg\n\n", c);
		if ((m_ctol > 0.0) && (c >= m_ctol))
		{
			feLogWarning("Condition number too large: %lg\nSwitching to secondary solver.", c);
			m_activeSolver = m_solver2;
		}
	}

	return m_activeSolver->Factor();
}

//! Backsolve the linear system
bool StrategySolver::BackSolve(double* x, double* b)
{
	assert(m_activeSolver);
	bool success = m_activeSolver->BackSolve(x, b);
	if ((success == false) && (m_activeSolver != m_solver2))
	{
		feLogError("Primary linear solver failed. Trying secondary linear solver ...");
		m_activeSolver = m_solver2;
		m_solver2->Factor();
		success = m_solver2->BackSolve(x, b);

		if (m_strategy == DONT_PERSIST)
		{
			m_solver2->Destroy();
			m_activeSolver = m_solver1;
		}
	}
	return success;
}

//! Clean up
void StrategySolver::Destroy()
{
	m_solver1->Destroy();
	m_solver2->Destroy();
	if (m_strategy != PERSIST_FOREVER) m_activeSolver = nullptr;
}

//! Create a sparse matrix
SparseMatrix* StrategySolver::CreateSparseMatrix(Matrix_Type ntype)
{
	// let solver1 create the sparse system
	if (m_solver1 == nullptr) return nullptr;
	SparseMatrix* A = m_solver1->CreateSparseMatrix(ntype);
	if (A == nullptr) return nullptr;

	// see if this format works for the second solver
	if ((m_solver2 == nullptr) || (m_solver2->SetSparseMatrix(A) == false)) {
		delete A;
		return nullptr;
	}

	m_pA = A;

	return A;
}

//! set the sparse matrix
bool StrategySolver::SetSparseMatrix(SparseMatrix* A)
{
	if ((m_solver1 == nullptr) || (m_solver1->SetSparseMatrix(A) == false)) return false;
	if ((m_solver2 == nullptr) || (m_solver2->SetSparseMatrix(A) == false)) return false;
	return true;
}
