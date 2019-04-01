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
#include "LinearSolver.h"

REGISTER_SUPER_CLASS(LinearSolver, FELINEARSOLVER_ID);

//-----------------------------------------------------------------------------
LinearSolver::LinearSolver(FEModel* fem) : FECoreBase(fem)
{
	ResetStats();
}

//-----------------------------------------------------------------------------
LinearSolver::~LinearSolver()
{ 
}

//-----------------------------------------------------------------------------
bool LinearSolver::PreProcess()
{ 
	return true; 
}

//-----------------------------------------------------------------------------
void LinearSolver::SetPartitions(const vector<int>& part)
{
	m_part = part;
}

//-----------------------------------------------------------------------------
void LinearSolver::SetPartitions(int npart0, int npart1)
{
	m_part.resize(2);
	m_part[0] = npart0;
	m_part[1] = npart1;
}

//-----------------------------------------------------------------------------
const LinearSolverStats& LinearSolver::GetStats() const
{
	return	m_stats;
}

//-----------------------------------------------------------------------------
void LinearSolver::ResetStats()
{
	m_stats.backsolves = 0;
	m_stats.iterations = 0;
}

//-----------------------------------------------------------------------------
void LinearSolver::UpdateStats(int iterations)
{
	m_stats.backsolves++;
	m_stats.iterations += iterations;
}

//-----------------------------------------------------------------------------
void LinearSolver::Destroy()
{

}

//-----------------------------------------------------------------------------
bool LinearSolver::SetSparseMatrix(SparseMatrix* pA)
{
	assert(false);
	return false;
}

//-----------------------------------------------------------------------------
//! convenience function for solving linear systems
bool LinearSolver::Solve(vector<double>& x, vector<double>& y)
{
	if (PreProcess() == false) return false;
	if (Factor() == false) return false;
	if (BackSolve(x, y) == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
//! convenience function for solving linear systems with an iterative solver
bool IterativeLinearSolver::Solve(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b, Preconditioner* pc)
{
	if (SetSparseMatrix(&A) == false) return false;
	SetPreconditioner(pc);
	if (PreProcess() == false) return false;
	if (Factor() == false) return false;
	return BackSolve(x, b);
}
