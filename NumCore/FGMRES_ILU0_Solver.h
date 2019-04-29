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
#include "FGMRESSolver.h"
#include <FECore/Preconditioner.h>
#include "ILU0_Preconditioner.h"

class CRSSparseMatrix;

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver with
//! ILU0 pre-conditioner for nonsymmetric indefinite matrices.
class FGMRES_ILU0_Solver : public FGMRESSolver
{
public:
	//! constructor
	FGMRES_ILU0_Solver(FEModel* fem);

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	// this is used to build the preconditioner
	bool Factor() override;

public: // preconditioner settings

	// do the zero diagonal check during preconditioner
	void DoZeroDiagonalCheck(bool b);

	// Set the zero diagonal tolerance value
	void SetZeroDiagonalTolerance(double tol);

	// set the zero diagonal replacement value
	void SetZeroDiagonalReplacement(double val);

private:
	ILU0_Preconditioner*	m_PC;		//!< the preconditioner
};

//-----------------------------------------------------------------------------
class ILU0_Solver : public LinearSolver
{
public:
	ILU0_Solver(FEModel* fem);
	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* y) override;

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
	bool SetSparseMatrix(SparseMatrix* pA) override;

private:
	Preconditioner*		m_PC;
};
