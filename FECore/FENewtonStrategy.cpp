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
#include "FENewtonStrategy.h"
#include "FENewtonSolver.h"
#include "LinearSolver.h"

FENewtonStrategy::FENewtonStrategy(FENewtonSolver* pns) : m_pns(pns)
{
	m_maxups = 10;
	m_max_buf_size = 0; // when zero, it should default to m_maxups
	m_cycle_buffer = true;

	m_nups = 0;
}

FENewtonStrategy::~FENewtonStrategy()
{
}

//! initialize the linear system
SparseMatrix* FENewtonStrategy::CreateSparseMatrix(Matrix_Type mtype)
{
	if (m_pns == 0) return 0;

	LinearSolver* plinsolve = m_pns->m_plinsolve;

	SparseMatrix* pS = plinsolve->CreateSparseMatrix(mtype);

	return pS;
}

bool FENewtonStrategy::ReformStiffness()
{
	return m_pns->ReformStiffness();
}

//! calculate the residual
bool FENewtonStrategy::Residual(std::vector<double>& R, bool binit)
{
	return m_pns->Residual(R);
}
