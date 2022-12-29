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
#include "SkylineSolver.h"

//-----------------------------------------------------------------------------
void colsol_factor(int N, double* values, int* pointers);
void colsol_solve(int N, double* values, int* pointers, double* R);

//-----------------------------------------------------------------------------
SkylineSolver::SkylineSolver(FEModel* fem) : LinearSolver(fem), m_pA(0)
{
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* SkylineSolver::CreateSparseMatrix(Matrix_Type ntype)
{ 
	return (m_pA = (ntype == REAL_SYMMETRIC? new SkylineMatrix() : 0)); 
}

//-----------------------------------------------------------------------------
bool SkylineSolver::PreProcess()
{
	// We don't need to do any preprocessing for this solver
	return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool SkylineSolver::Factor()
{
	colsol_factor(m_pA->Rows(), m_pA->values(), m_pA->pointers());
	return true;
}

//-----------------------------------------------------------------------------
bool SkylineSolver::BackSolve(double* x, double* b)
{
	// we need to make a copy of R since colsol overwrites the right hand side vector
	// with the solution
	int neq = m_pA->Rows();
	for (int i=0; i<neq; ++i) x[i] = b[i];
	colsol_solve(m_pA->Rows(), m_pA->values(), m_pA->pointers(), x);

	return true;
}

//-----------------------------------------------------------------------------
void SkylineSolver::Destroy()
{
	// Nothing to destroy
	LinearSolver::Destroy();
}
