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
#include <FECore/CompactUnSymmMatrix.h>

//-----------------------------------------------------------------------------
// This class implements the HYPRE GMRES solver
class HypreGMRESsolver : public LinearSolver
{
	class Implementation;

public:
	HypreGMRESsolver(FEModel* fem);
	~HypreGMRESsolver();

	void SetPrintLevel(int n) override;

	void SetMaxIterations(int n);

	void SetConvergencTolerance(double tol);

public:
	// allocate storage
	bool PreProcess() override;

	//! Factor the matrix (for iterative solvers, this can be used for creating pre-conditioner)
	bool Factor() override;

	//! Calculate the solution of RHS b and store solution in x
	bool BackSolve(double* x, double* b) override;

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* A) override;

	//! clean up
	void Destroy() override;

private:
	Implementation*	imp;

	DECLARE_FECORE_CLASS();
};
