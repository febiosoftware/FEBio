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

#pragma once
#include <FECore/LinearSolver.h>

// This class implements a linear solver that can switch between linear solvers
class StrategySolver : public LinearSolver
{
	// the persist strategy determines when to switch back from solver2 to solver1
	// after solver1 fails. 
	// DONT_PERSIST (0) = switch back immediately.
	// PERSIST          = switch back on the next Destroy()
	// PERSIST_FOREVER  = switch to solver2 permanently
	enum PersistStrategy
	{
		DONT_PERSIST,
		PERSIST,
		PERSIST_FOREVER
	};

public:
	StrategySolver(FEModel* fem);

	~StrategySolver();

public:
	//! Preprocess 
	bool PreProcess() override;

	//! Factor matrix
	bool Factor() override;

	//! Backsolve the linear system
	bool BackSolve(double* x, double* b) override;

	//! Clean up
	void Destroy() override;

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* A) override;

private:
	int				m_strategy;
	double			m_ctol;

	bool			m_print_cn;	// print the condition number
	LinearSolver*	m_solver1;	// the primary solver
	LinearSolver*	m_solver2;	// the seoncdary solver
	LinearSolver*	m_activeSolver;

	SparseMatrix*	m_pA;

	DECLARE_FECORE_CLASS();
};
