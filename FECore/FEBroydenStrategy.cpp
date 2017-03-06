#include "stdafx.h"
#include "FEBroydenStrategy.h"
#include "LinearSolver.h"
#include "FEException.h"

//-----------------------------------------------------------------------------
//! constructor
FEBroydenStrategy::FEBroydenStrategy()
{
}

//-----------------------------------------------------------------------------
//! Initialization
void FEBroydenStrategy::Init(int neq, LinearSolver* pls)
{
	if (m_max_buf_size <= 0) m_max_buf_size = m_maxups;

	// allocate storage for Broyden update vectors
	m_R.resize(m_max_buf_size, neq);
	m_D.resize(m_max_buf_size, neq);
	m_rho.resize(m_max_buf_size);

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
	if (m_plinsolve->BackSolve(q, R1) == false)
		throw LinearSolverFailed();

	// only update when allowed
	if ((m_nups < m_max_buf_size) || (m_cycle_buffer == true))
	{
		// since the number of updates can be larger than buffer size we need to clamp
		int nups = (m_nups >= m_max_buf_size ? m_max_buf_size - 1: m_nups);

		// get the "0" buffer index
		int n0 = (m_nups >= m_max_buf_size ? (m_nups + 1) % m_max_buf_size : 0);
		int n1 = (m_nups >= m_max_buf_size ? (m_nups    ) % m_max_buf_size : m_nups);

		// loop over update vectors
		for (int j=0; j<nups; ++j)
		{
			int n = (n0 + j) % m_max_buf_size;

			double w = 0.0;
			for (int i=0; i<m_neq; ++i) w += (m_D[n][i]*q[i]);

			double g = m_rho[n]*w;
			for (int i=0; i<m_neq; ++i) q[i] += g*(m_D[n][i] - m_R[n][i]);
		}

		// form and store the next update vector
		double rhoi = 0.0;
		for (int i=0; i<m_neq; ++i)
		{
			double ri = q[i] - ui[i];
			double di = -s*ui[i];
			m_R[n1][i] = ri;
			m_D[n1][i] = di;

			rhoi += di*ri;
		}
		m_rho[n1] = 1.0 / (rhoi);
	}

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
		if (m_plinsolve->BackSolve(x, b) == false)
			throw LinearSolverFailed();
	}
	else
	{
		// calculate q1
		vector<double> q(m_neq);
		if (m_plinsolve->BackSolve(q, b) == false)
			throw LinearSolverFailed();

		// since the number of updates can be larger than buffer size we need to clamp
		int nups = (m_nups > m_max_buf_size ? m_max_buf_size : m_nups);

		// get the first and last buffer index
		int n0 = 0;
		int n1 = nups - 1;
		if ((m_nups > m_max_buf_size) && (m_cycle_buffer == true))
		{
			n0 = m_nups % m_max_buf_size;
			n1 = (m_nups - 1) % m_max_buf_size;
		}

		for (int j=0; j<nups-1; ++j)
		{
			int n = (n0 + j) % m_max_buf_size;

			double w = 0.0;
			for (int i=0; i<m_neq; ++i) w += (m_D[n][i]*q[i]);

			double g = m_rho[n]*w;
			for (int i=0; i<m_neq; ++i) q[i] += g*(m_D[n][i] - m_R[n][i]);
		}

		// calculate solution
		double rho = 0.0;
		for (int i=0; i<m_neq; ++i) rho += m_D[n1][i]*q[i];
		rho *= m_rho[n1];
	
		for (int i=0; i<m_neq; ++i)
		{
			x[i] = q[i] + rho*(m_D[n1][i] - m_R[n1][i]);
		}
	}
}
