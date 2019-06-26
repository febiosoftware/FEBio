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
#pragma once
#include <FECore/Preconditioner.h>
#include "PardisoSolver.h"
#include "FGMRESSolver.h"

class BlockPreconditioner : public Preconditioner
{
public:
	BlockPreconditioner(FEModel* fem);

	~BlockPreconditioner();

	// create a preconditioner for a sparse matrix
	bool Create() override;

	// apply to vector P x = y
	bool mult_vector(double* x, double* y) override;

	// set the solution method
	void SetSolutionMethod(int method);

private:
	std::vector<PardisoSolver*>	m_solver;
	int m_method; // 0 = Jacobian, 1 = forward Gauss-Seidel, 2 = backward Gauss-Seidel
};


//-----------------------------------------------------------------------------
class FGMRES_Jacobi_Block_Solver : public FGMRESSolver
{
public:
	//! constructor
	FGMRES_Jacobi_Block_Solver(FEModel* fem);

	~FGMRES_Jacobi_Block_Solver();

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* A) override;

	// this is used to build the preconditioner
	bool Factor() override;

	// clean up
	void Destroy() override;

// preconditioner settings
public:

	// set the solution method
	void SetSolutionMethod(int method);

private:
	BlockPreconditioner*	m_PC;		//!< the preconditioner
};
