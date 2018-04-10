#include "stdafx.h"
#include "JFNKStrategy.h"
#include "FENewtonSolver.h"
#include <NumCore/FGMRESSolver.h>
#include "FEException.h"

//=================================================================================================
// JFNKMatrix
//=================================================================================================

//-----------------------------------------------------------------------------
// This is a class that just mimics a sparse matrix.
// It is only used by the JFNK strategy. 
// The only function it implements is the mult_vector.
class JFNKMatrix : public SparseMatrix
{
public:
	JFNKMatrix(FENewtonSolver* pns) : m_pns(pns)
	{
		m_ndim = pns->m_neq;

		// TODO: For contact problems we'll need some mechanism to change the array size
		m_v.resize(m_ndim);
		m_R.resize(m_ndim);
	}

	void zero() override { assert(false); }
	void Create(SparseMatrixProfile& MP) override { assert(false); }
	void Assemble(matrix& ke, std::vector<int>& lm) override { assert(false); }
	void Assemble(matrix& ke, std::vector<int>& lmi, std::vector<int>& lmj) override { assert(false); }
	void set(int i, int j, double v) override { assert(false); }
	void add(int i, int j, double v) override { assert(false); }
	double get(int i, int j) override { assert(false); return 0.0; }
	double diag(int i) override { assert(false); return 0.0; }
	void Clear() {}

	//! multiply with vector
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
}

//! New initialization method
void JFNKStrategy::Init(int neq, LinearSolver* pls)
{
	m_plinsolve = pls;
}

SparseMatrix* JFNKStrategy::CreateSparseMatrix(Matrix_Type mtype)
{
	// make sure the solver is of FGMRES type
	FGMRESSolver* fgmres = dynamic_cast<FGMRESSolver*>(m_pns->m_plinsolve);
	if (fgmres == 0) return 0;

	// create the sparse matrix facade
	JFNKMatrix* pA = new JFNKMatrix(m_pns);

	// set the matrix
	fgmres->SetSparseMatrix(pA);
	fgmres->PreProcess();

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
	// nothing to do here
	return true;
}
