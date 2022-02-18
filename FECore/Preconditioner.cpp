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
#include "Preconditioner.h"

//=================================================================================================
DiagonalPreconditioner::DiagonalPreconditioner(FEModel* fem) : Preconditioner(fem)
{
	m_bsqr = false;
}

// take square root of diagonal entries
void DiagonalPreconditioner::CalculateSquareRoot(bool b)
{
	m_bsqr = b;
}

// create a preconditioner for a sparse matrix
bool DiagonalPreconditioner::Factor()
{
	SparseMatrix* A = GetSparseMatrix();
	if (A == nullptr) return false;

	int N = A->Rows();
	if (A->Columns() != N) return false;

	m_D.resize(N);
	for (int i=0; i<N; ++i)
	{
		double dii = A->diag(i);
		if (dii == 0.0) return false;
		if (m_bsqr)
		{
			if (dii < 0) return false;
			dii = sqrt(dii);
		}
		m_D[i] = 1.0 / dii;
	}

	return true;
}

// apply to vector P x = y
bool DiagonalPreconditioner::BackSolve(double* x, double* y)
{
	int N = (int)m_D.size();

#pragma omp parallel for
	for (int i=0; i<N; ++i)
	{
		x[i] = y[i]*m_D[i];
	}

	return true;
}
