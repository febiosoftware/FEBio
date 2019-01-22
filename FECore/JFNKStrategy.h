#pragma once
#include "FENewtonStrategy.h"
#include "SparseMatrix.h"

class JFNKMatrix;

//-----------------------------------------------------------------------------
// Implements a Jacobian-Free Newton-Krylov strategy
class JFNKStrategy : public FENewtonStrategy
{
public:
	JFNKStrategy(FENewtonSolver* pns);

	//! New initialization method
	void Init(int neq, LinearSolver* pls) override;

	//! initialize the linear system
	SparseMatrix* CreateSparseMatrix(Matrix_Type mtype) override;

	//! perform a BFGS udpate
	bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1) override;

	//! solve the equations
	void SolveEquations(vector<double>& x, vector<double>& b) override;

	//! Overide reform stiffness because we don't want to do any reformations
	bool ReformStiffness() override;

	//! override so we can store a copy of the residual before we add Fd
	bool Residual(std::vector<double>& R, bool binit) override;

public:
	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;		//!< pointer to linear solver
	int				m_neq;				//!< number of equations
	bool			m_bprecondition;	//!< the solver requires preconditioning

	JFNKMatrix*		m_A;
};
