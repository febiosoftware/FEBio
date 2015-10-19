#include "stdafx.h"
#include "FEBroydenStrategy.h"
#include "LinearSolver.h"

//-----------------------------------------------------------------------------
//! constructor
FEBroydenStrategy::FEBroydenStrategy()
{
}

//-----------------------------------------------------------------------------
//! Initialization
void FEBroydenStrategy::Init(int neq, LinearSolver* pls)
{
	// allocate storage for Broyden update vectors
	m_R.resize(m_maxups, neq);
	m_D.resize(m_maxups, neq);
	m_rho.resize(m_maxups);

	m_neq = neq;
	m_nups = 0;

	m_plinsolve = pls;
}

//-----------------------------------------------------------------------------
//! perform a quasi-Newton udpate
bool FEBroydenStrategy::Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
	// calculate q1
	vector<double> q(m_neq);
	m_plinsolve->BackSolve(q, R1);

	// loop over update vectors
	for (int j=0; j<m_nups; ++j)
	{
		double w = 0.0;
		for (int i=0; i<m_neq; ++i) w += (m_D[j][i]*q[i]);

		double g = m_rho[j]*w;
		for (int i=0; i<m_neq; ++i) q[i] += g*(m_D[j][i] - m_R[j][i]);
	}

	// form and store the next update vector
	double rhoi = 0.0;
	for (int i=0; i<m_neq; ++i)
	{
		double ri = q[i] - ui[i];
		double di = -s*ui[i];
		m_R[m_nups][i] = ri;
		m_D[m_nups][i] = di;

		rhoi += di*ri;
	}
	m_rho[m_nups] = 1.0 / (rhoi);

	m_nups++;

	return true;
}

//-----------------------------------------------------------------------------
//! solve the equations
void FEBroydenStrategy::SolveEquations(vector<double>& x, vector<double>& b)
{
	// loop over update vectors
	if (m_nups == 0)
	{
		m_plinsolve->BackSolve(x, b);
	}
	else
	{
		// calculate q1
		vector<double> q(m_neq);
		m_plinsolve->BackSolve(q, b);

		for (int j=0; j<m_nups-1; ++j)
		{
			double w = 0.0;
			for (int i=0; i<m_neq; ++i) w += (m_D[j][i]*q[i]);

			double g = m_rho[j]*w;
			for (int i=0; i<m_neq; ++i) q[i] += g*(m_D[j][i] - m_R[j][i]);
		}

		// calculate solution
		double rho = 0.0;
		for (int i=0; i<m_neq; ++i) rho += m_D[m_nups-1][i]*q[i];
		rho *= m_rho[m_nups-1];
	
		for (int i=0; i<m_neq; ++i)
		{
			x[i] = q[i] + rho*(m_D[m_nups-1][i] - m_R[m_nups-1][i]);
		}
	}
}
