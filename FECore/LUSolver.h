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
#include "LinearSolver.h"
#include "DenseMatrix.h"
#include "fecore_api.h"

//-----------------------------------------------------------------------------
//! LU decomposition solver

//! This solver performs an LU decomposition and uses a backsolving algorithm
//! to solve the equations.
//! This solver uses the FullMatrix class and therefore is not the preferred
//! solver. It should only be used for small problems and only when the other
//! solvers are not adequate.

class FECORE_API LUSolver : public LinearSolver
{
public:
	//! constructor
	LUSolver(FEModel* fem = nullptr);

	//! Pre-process data
	bool PreProcess() override;

	//! Factor matrix
	bool Factor() override;

	//! solve using factored matrix
	bool BackSolve(double* x, double* b) override;

	//! Clean-up
	void Destroy() override;

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! Set the matrix
	void SetMatrix(FECore::DenseMatrix* pA);

protected:
	std::vector<int>		indx;	//!< indices
	FECore::DenseMatrix*	m_pA;	//!< sparse matrix
};
