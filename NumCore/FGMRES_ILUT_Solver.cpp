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

#include "stdafx.h"
#include "FGMRES_ILUT_Solver.h"

//=============================================================================
FGMRES_ILUT_Solver::FGMRES_ILUT_Solver(FEModel* fem) : FGMRESSolver(fem)
{
	// set the preconditioner
	SetPreconditioner(m_PC = new ILUT_Preconditioner(fem));
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_ILUT_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
	// We can only support non-symmetric matrices because of the ILUT preconditioner
	if (ntype != REAL_UNSYMMETRIC) return 0;
	SparseMatrix* A = FGMRESSolver::CreateSparseMatrix(ntype);
	m_PC->SetSparseMatrix(A);
	return A;
}

//-----------------------------------------------------------------------------
// set the max fill value of the preconditioner
void FGMRES_ILUT_Solver::SetMaxFill(int n)
{
	m_PC->m_maxfill = n;
}

//-----------------------------------------------------------------------------
// Set the fill tolerance of the preconditioner
void FGMRES_ILUT_Solver::SetFillTolerance(double fillTol)
{
	m_PC->m_fillTol = fillTol;
}

//-----------------------------------------------------------------------------
// do the zero diagonal check during preconditioner
void FGMRES_ILUT_Solver::DoZeroDiagonalCheck(bool b)
{
	m_PC->m_checkZeroDiagonal = b;
}

//-----------------------------------------------------------------------------
// Set the zero diagonal tolerance value of the preconditioner
void FGMRES_ILUT_Solver::SetZeroDiagonalTolerance(double tol)
{
	m_PC->m_zeroThreshold = tol;
}

//-----------------------------------------------------------------------------
// set the zero diagonal replacement value of the preconditioner
void FGMRES_ILUT_Solver::SetZeroDiagonalReplacement(double val)
{
	m_PC->m_zeroReplace = val;
}

//-----------------------------------------------------------------------------
// this is used to build the preconditioner
bool FGMRES_ILUT_Solver::Factor()
{
	if (FGMRESSolver::Factor() == false) return false;
	return m_PC->Create();
}
