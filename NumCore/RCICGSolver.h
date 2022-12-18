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
#include <FECore/CompactSymmMatrix.h>

// This class implements an interface to the RCI CG iterative solver from the MKL math library.
class RCICGSolver : public IterativeLinearSolver
{
public:
	RCICGSolver(FEModel* fem);
	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* b) override;
	void Destroy() override;

public:
	bool HasPreconditioner() const override;

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	bool SetSparseMatrix(SparseMatrix* A) override;

	void SetLeftPreconditioner(LinearSolver* P) override;
	LinearSolver* GetLeftPreconditioner() override;

	void SetMaxIterations(int n) { m_maxiter = n; }
	void SetTolerance(double tol) { m_tol = tol; }
	void SetPrintLevel(int n) override { m_print_level = n; }

protected:
	SparseMatrix*		m_pA;
	LinearSolver*		m_P;

	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level
	bool	m_fail_max_iters;

	DECLARE_FECORE_CLASS();
};
