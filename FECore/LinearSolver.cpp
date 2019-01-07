// LinSolver.cpp: implementation of the LinSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "LinearSolver.h"

REGISTER_SUPER_CLASS(LinearSolver, FELINEARSOLVER_ID);

//-----------------------------------------------------------------------------
LinearSolver::LinearSolver(FEModel* fem) : FECoreBase(fem, FELINEARSOLVER_ID)
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
