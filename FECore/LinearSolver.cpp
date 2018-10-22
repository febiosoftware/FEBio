// LinSolver.cpp: implementation of the LinSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "LinearSolver.h"

//-----------------------------------------------------------------------------
LinearSolver::LinearSolver(FEModel* fem) : m_fem(fem)
{

}

//-----------------------------------------------------------------------------
LinearSolver::~LinearSolver()
{ 
	Destroy(); 
}

//-----------------------------------------------------------------------------
FEModel* LinearSolver::GetFEModel() const
{
	return m_fem;
}

//-----------------------------------------------------------------------------
bool LinearSolver::PreProcess()
{ 
	return true; 
}

//-----------------------------------------------------------------------------
void LinearSolver::SetPartition(int nsplit)
{
}

//-----------------------------------------------------------------------------
void LinearSolver::SetPartitions(const vector<int>& part)
{
	
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
