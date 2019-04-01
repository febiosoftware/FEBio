/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "ILUT_Preconditioner.h"

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver with
//! ILUT pre-conditioner for nonsymmetric indefinite matrices.
class FGMRES_ILUT_Solver : public FGMRESSolver
{
public:
	//! constructor
	FGMRES_ILUT_Solver(FEModel* fem);

	//! Return a sparse matrix compatible with this solver
	//! This is overridden because this solver will only work with nonsymmetric matrices
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	// this is used to build the preconditioner
	bool Factor() override;

public: // Preconditioner properties

	// set the max fill value
	void SetMaxFill(int n);

	// Set the fill tolerance
	void SetFillTolerance(double fillTol);

	// do the zero diagonal check during preconditioner
	void DoZeroDiagonalCheck(bool b);

	// Set the zero diagonal tolerance value
	void SetZeroDiagonalTolerance(double tol);

	// set the zero diagonal replacement value
	void SetZeroDiagonalReplacement(double val);

private:
	ILUT_Preconditioner*	m_PC;		//!< the preconditioner
};
