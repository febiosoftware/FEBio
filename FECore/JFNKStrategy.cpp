#include "stdafx.h"
#include "JFNKStrategy.h"
#include "FENewtonSolver.h"
#include "JFNKMatrix.h"
#include "FEException.h"

JFNKStrategy::JFNKStrategy(FENewtonSolver* pns) : FENewtonStrategy(pns)
{
	m_bprecondition = false;
}

//! New initialization method
void JFNKStrategy::Init(int neq, LinearSolver* pls)
{
	m_plinsolve = pls;
}

SparseMatrix* JFNKStrategy::CreateSparseMatrix(Matrix_Type mtype)
{
	JFNKMatrix* pA = 0;

	// make sure the linear solver is an iterative linear solver
	IterativeLinearSolver* ls = dynamic_cast<IterativeLinearSolver*>(m_pns->m_plinsolve);
	if (ls)
	{
		// see if this solver has a preconditioner
		m_bprecondition = ls->HasPreconditioner();

		// if the solver has a preconditioner, we still need to create the stiffness matrix
		SparseMatrix* K = 0;
		if (m_bprecondition) 
		{
			K = ls->CreateSparseMatrix(mtype);
			if (K == 0) return 0;
		}

		// Now, override the matrix used
		pA = new JFNKMatrix(m_pns, K);
		ls->SetSparseMatrix(pA);

		// If there is no preconditioner we can do the pre-processing here
		if (m_bprecondition == false) ls->PreProcess();
	}

	return pA;
}

//! perform a BFGS udpate
bool JFNKStrategy::Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
	// nothing to do here
	return true;
}

//! solve the equations
void JFNKStrategy::SolveEquations(vector<double>& x, vector<double>& b)
{
	// perform a backsubstitution
	if (m_plinsolve->BackSolve(x, b) == false)
	{
		throw LinearSolverFailed();
	}
}

bool JFNKStrategy::ReformStiffness()
{
	if (m_bprecondition) return m_pns->ReformStiffness();
	else return true;
}
