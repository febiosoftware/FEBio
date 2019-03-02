// LinSolver.cpp: implementation of the LinSolver class.
//
//////////////////////////////////////////////////////////////////////

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
	m_stats.memsize = 0;
}

//-----------------------------------------------------------------------------
void LinearSolver::UpdateStats(int iterations)
{
	m_stats.backsolves++;
	m_stats.iterations += iterations;
}

//-----------------------------------------------------------------------------
void LinearSolver::SetMemoryUsage(size_t memsize)
{
	m_stats.memsize = memsize;
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
