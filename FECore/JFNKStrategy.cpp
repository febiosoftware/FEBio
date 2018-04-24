#include "stdafx.h"
#include "JFNKStrategy.h"
#include "FENewtonSolver.h"
#include <NumCore/FGMRESSolver.h>
#include <NumCore/FGMRES_ILU0_Solver.h>
#include "FEException.h"

//=================================================================================================
// JFNKMatrix
//=================================================================================================

//-----------------------------------------------------------------------------
// This is a class that just mimics a sparse matrix.
// It is only used by the JFNK strategy. 
// The only function it implements is the mult_vector.
class JFNKMatrix : public CompactUnSymmMatrix
{
public:
	JFNKMatrix(FENewtonSolver* pns) : CompactUnSymmMatrix(1, true), m_pns(pns)
	{
		m_ndim = pns->m_neq;

		// TODO: For contact problems we'll need some mechanism to change the array size
		m_v.resize(m_ndim);
		m_R.resize(m_ndim);
	}

	//! override multiply with vector
	void mult_vector(double* x, double* r) override;

private:
	FENewtonSolver*	m_pns;
	vector<double>	m_v, m_R;
};

void JFNKMatrix::mult_vector(double* x, double* r)
{
	double eps = 0.001;
	int neq = (int)m_pns->m_ui.size();

	for (int i=0; i<neq; ++i) m_v[i] = eps*x[i];

	m_pns->Update(m_v);
	m_pns->Residual(m_R);

	for (int i=0; i<neq; ++i)
	{
		r[i] = (m_pns->m_R0[i] - m_R[i])/eps;
	}
}

//=================================================================================================
// JFNKStrategy
//=================================================================================================

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

	// make sure the solver is of FGMRES type
	FGMRESSolver* fgmres = dynamic_cast<FGMRESSolver*>(m_pns->m_plinsolve);
	if (fgmres)
	{
		pA = new JFNKMatrix(m_pns);

		// set the matrix
		fgmres->SetSparseMatrix(pA);
		fgmres->PreProcess();

		m_bprecondition = false;
	}

	FGMRES_ILU0_Solver* fgmres_ilu0 = dynamic_cast<FGMRES_ILU0_Solver*>(m_pns->m_plinsolve);
	if (fgmres_ilu0)
	{
		pA = new JFNKMatrix(m_pns);

		fgmres_ilu0->SetSparseMatrix(pA);
		m_bprecondition = true;
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
