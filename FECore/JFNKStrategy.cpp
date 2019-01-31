#include "stdafx.h"
#include "JFNKStrategy.h"
#include "FENewtonSolver.h"
#include "JFNKMatrix.h"
#include "FEException.h"
#include "LinearSolver.h"
#include "Preconditioner.h"
#include "FEModel.h"
#include "log.h"

JFNKStrategy::JFNKStrategy(FENewtonSolver* pns) : FENewtonStrategy(pns)
{
	m_bprecondition = false;
	m_A = nullptr;
}

//! New initialization method
void JFNKStrategy::Init(int neq, LinearSolver* pls)
{
	m_plinsolve = pls;

	// we have to turn off the line search for now
	FELineSearch* ls = m_pns->GetLineSearch();
	if (ls && (ls->m_LStol > 0.0))
	{
		felog.printbox("WARNING", "Line search is turned off because the JFNK solver currently does not support it.");
		ls->m_LStol = 0.0;
	}
}

SparseMatrix* JFNKStrategy::CreateSparseMatrix(Matrix_Type mtype)
{
	if (m_A) delete m_A;
	m_A = 0;

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
		m_A = new JFNKMatrix(m_pns, K);
		ls->SetSparseMatrix(m_A);

		// If there is no preconditioner we can do the pre-processing here
		if (m_bprecondition == false)
		{
			ls->PreProcess();
		}
		else
		{
			ls->GetPreconditioner()->SetSparseMatrix(K);
		}
	}

	return m_A;
}

//! perform a quasi-Newton udpate
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

//! override so we can store a copy of the residual before we add Fd
bool JFNKStrategy::Residual(std::vector<double>& R, bool binit)
{
	// first calculate the residual
	bool b = m_pns->Residual(R);
	if (b == false) return false;

	// store a copy
	m_A->SetReferenceResidual(R);

	// at the first iteration we need to calculate the Fd
	if (binit)
	{
		// get the vector of prescribed displacements
		std::vector<double>& ui = m_pns->m_ui;

		// build m_Fd
		vector<double>& Fd = m_pns->m_Fd;
		m_A->SetPolicy(JFNKMatrix::ZERO_FREE_DOFS);
		m_A->mult_vector(&ui[0], &Fd[0]);
		m_A->SetPolicy(JFNKMatrix::ZERO_PRESCRIBED_DOFS);

		// we need to flip the sign on Fd
		for (size_t i = 0; i < Fd.size(); ++i) Fd[i] = -Fd[i];

		// NOTE: Usually, the residual is called after an Update. However, here, 
		//       the model is updated (via Update2) in m_A->mult_vector. I wonder
		//       if I need to restore the model's state to before the mult_vector call
		//       by calling another Update2 method with a zero vector.
	}

	return true;
}
