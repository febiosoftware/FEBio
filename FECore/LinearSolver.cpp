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
void LinearSolver::Destroy()
{

}
