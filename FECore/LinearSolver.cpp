// LinSolver.cpp: implementation of the LinSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "LinearSolver.h"

//-----------------------------------------------------------------------------
LinearSolver::LinearSolver()
{

}

//-----------------------------------------------------------------------------
LinearSolver::~LinearSolver()
{ 
	Destroy(); 
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
void LinearSolver::SetSparseMatrix(SparseMatrix* pA)
{
	assert(false);
}
