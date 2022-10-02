/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "stdafx.h"
#include "FEBroydenStrategy.h"
#include "LinearSolver.h"
#include "FEException.h"
#include "FENewtonSolver.h"
#include "log.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEBroydenStrategy, FENewtonStrategy)
	ADD_PARAMETER(m_maxups, "max_ups");
	ADD_PARAMETER(m_max_buf_size, FE_RANGE_GREATER_OR_EQUAL(0), "max_buffer_size"); 
	ADD_PARAMETER(m_cycle_buffer, "cycle_buffer");
	ADD_PARAMETER(m_cmax, FE_RANGE_GREATER_OR_EQUAL(0.0), "cmax");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEBroydenStrategy::FEBroydenStrategy(FEModel* fem) : FENewtonStrategy(fem)
{
	m_neq = 0;
	m_plinsolve = nullptr;
}

//-----------------------------------------------------------------------------
//! Initialization
bool FEBroydenStrategy::Init()
{
	if (m_pns == nullptr) return false;

	if (m_max_buf_size <= 0) m_max_buf_size = m_maxups;

	int neq = m_pns->m_neq;

	// allocate storage for Broyden update vectors
	m_R.resize(m_max_buf_size, neq);
	m_D.resize(m_max_buf_size, neq);
	m_rho.resize(m_max_buf_size);
	m_q.resize(neq, 0.0);

	m_neq = neq;
	m_nups = 0;

	m_plinsolve = m_pns->GetLinearSolver();

	m_bnewStep = true;

	return true;
}

//! Presolve update
void FEBroydenStrategy::PreSolveUpdate()
{
	m_bnewStep = true;
}

//-----------------------------------------------------------------------------
//! perform a quasi-Newton udpate
bool FEBroydenStrategy::Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
	// for full-Newton, we skip QN update
	if (m_maxups == 0) return false;

	// make sure we didn't reach max updates
	if (m_nups >= m_maxups - 1)
	{
		feLogWarning("Max nr of iterations reached.\nStiffness matrix will now be reformed.");
		return false;
	}

	// calculate q1
	m_q.assign(m_neq, 0.0);
	if (m_plinsolve->BackSolve(m_q, R1) == false)
		throw LinearSolverFailed();

	// only update when allowed
	if ((m_nups < m_max_buf_size) || (m_cycle_buffer == true))
	{
		// since the number of updates can be larger than buffer size we need to clamp
		int nups = (m_nups >= m_max_buf_size ? m_max_buf_size - 1 : m_nups);

		// get the "0" buffer index
		int n0 = (m_nups >= m_max_buf_size ? (m_nups + 1) % m_max_buf_size : 0);
		int n1 = (m_nups >= m_max_buf_size ? (m_nups) % m_max_buf_size : m_nups);

		// loop over update vectors
		for (int j = 0; j<nups; ++j)
		{
			int n = (n0 + j) % m_max_buf_size;

			double w = 0.0;
			for (int i = 0; i<m_neq; ++i) w += (m_D[n][i] * m_q[i]);

			double g = m_rho[n] * w;
			for (int i = 0; i<m_neq; ++i) m_q[i] += g*(m_D[n][i] - m_R[n][i]);
		}

		// form and store the next update vector
		double rhoi = 0.0;
		for (int i = 0; i<m_neq; ++i)
		{
			double ri = m_q[i] - ui[i];
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

		m_bnewStep = false;
	}
	else
	{
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

		if (m_bnewStep)
		{
			m_q = x;
			if (m_plinsolve->BackSolve(m_q, b) == false)
				throw LinearSolverFailed();

			for (int j = 0; j<nups - 1; ++j)
			{
				int n = (n0 + j) % m_max_buf_size;

				double w = 0.0;
				for (int i = 0; i<m_neq; ++i) w += (m_D[n][i] * m_q[i]);

				double g = m_rho[n] * w;
				for (int i = 0; i<m_neq; ++i) m_q[i] += g*(m_D[n][i] - m_R[n][i]);
			}

			m_bnewStep = false;
		}

		// calculate solution
		double rho = 0.0;
		for (int i = 0; i<m_neq; ++i) rho += m_D[n1][i] * m_q[i];
		rho *= m_rho[n1];

		for (int i = 0; i<m_neq; ++i)
		{
			x[i] = m_q[i] + rho*(m_D[n1][i] - m_R[n1][i]);
		}
	}
}

/*
//-----------------------------------------------------------------------------
//! perform a quasi-Newton udpate
bool FEBroydenStrategy::Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
	// calculate q1
	// TODO: q is not initialized. Don't the iterative solvers use this as an initial guess?
	// TODO: Can I use ui as an initial guess for q?
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
*/
