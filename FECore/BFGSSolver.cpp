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
#include "BFGSSolver.h"
#include "FESolver.h"
#include "FEException.h"
#include "FENewtonSolver.h"
#include "log.h"

//-----------------------------------------------------------------------------
// BFGSSolver
//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(BFGSSolver, FENewtonStrategy)
	ADD_PARAMETER(m_maxups, "max_ups");
	ADD_PARAMETER(m_max_buf_size, FE_RANGE_GREATER_OR_EQUAL(0), "max_buffer_size"); 
	ADD_PARAMETER(m_cycle_buffer, "cycle_buffer");
	ADD_PARAMETER(m_cmax, FE_RANGE_GREATER_OR_EQUAL(0.0), "cmax");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
BFGSSolver::BFGSSolver(FEModel* fem) : FENewtonStrategy(fem)
{
	m_neq = 0;

	// pointer to linear solver
	m_plinsolve = 0;
}

//-----------------------------------------------------------------------------
// Initialization method
bool BFGSSolver::Init()
{
	if (m_pns == nullptr) return false;

	if (m_max_buf_size <= 0) m_max_buf_size = m_maxups;

	int neq = m_pns->m_neq;

	// allocate storage for BFGS update vectors
	m_V.resize(m_max_buf_size, neq);
	m_W.resize(m_max_buf_size, neq);

	m_D.resize(neq);
	m_G.resize(neq);
	m_H.resize(neq);

	m_neq = neq;
	m_nups = 0;

	m_plinsolve = m_pns->GetLinearSolver();

	return true;
}

//-----------------------------------------------------------------------------
//! This function performs a BFGS stiffness update.
//! The last line search step is input to this function.
//! This function performs the update assuming the stiffness matrix
//! is positive definite. In the case that this is not the case
//! the function returns false. The function also returns false if 
//! the condition number is too big. A too big condition number might
//! be an indication of an ill-conditioned matrix and the update should
//! not be performed.

bool BFGSSolver::Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
	// for full-Newton, we skip QN update
	if (m_maxups == 0) return false;

	// make sure we didn't reach max updates
	if (m_nups >= m_maxups - 1)
	{
		feLogWarning("Max nr of iterations reached.\nStiffness matrix will now be reformed.");
		return false;
	}

	// calculate the BFGS update vectors
	int neq = m_neq;
	for (int i = 0; i<neq; ++i)
	{
		m_D[i] = s*ui[i];
		m_G[i] = R0[i] - R1[i];
		m_H[i] = R0[i]*s;
	}

	double dg = m_D*m_G;
	double dh = m_D*m_H;
	double dgi = 1.0 / dg;
	double r = dg / dh;

	// check to see if this is still a pos definite update
//	if (r <= 0) 
//	{
//		feLogWarning("The QN update has failed.\nStiffness matrix will now be reformed.");
//		return false;
//	}

	// calculate the condition number
//	c = sqrt(r);
	double c = sqrt(fabs(r));

	// make sure c is less than the the maximum.
	if (c > m_cmax) return false;

	// get the correct buffer to fill
	int n = (m_nups % m_max_buf_size);

	// do the update only when allowed
	if ((m_nups < m_max_buf_size) || (m_cycle_buffer == true))
	{
		double* vn = m_V[n];
		double* wn = m_W[n];

		for (int i=0; i<neq; ++i)	
		{
			vn[i] = -m_H[i]*c - m_G[i];
			wn[i] = m_D[i]*dgi;
		}
	}

	// increment update counter
	++m_nups;

	return true;
}

//-----------------------------------------------------------------------------
// This function solves a system of equations using the BFGS update vectors
// The variable m_nups keeps track of how many updates have been made so far.
// 

void BFGSSolver::SolveEquations(vector<double>& x, vector<double>& b)
{
	// make sure we need to do work
	if (m_neq ==0) return;

	// create temporary storage
	tmp = b;

	// number of updates can be larger than buffer size, so clamp it
	int nups = (m_nups> m_max_buf_size ? m_max_buf_size : m_nups);

	// get the "0" buffer index
	int n0 = 0;
	if ((m_nups > m_max_buf_size) && (m_cycle_buffer == true))
	{
		n0 = m_nups % m_max_buf_size;
	}

	// loop over all update vectors
	for (int i=nups-1; i>=0; --i)
	{
		int n = (n0 + i) % m_max_buf_size;

		double* vi = m_V[n];
		double* wi = m_W[n];

		double wr = 0;
		for (int j = 0; j<m_neq; j++) wr += wi[j] * tmp[j];

		for (int j = 0; j<m_neq; j++) tmp[j] += vi[j] * wr;
	}

	// perform a backsubstitution
	if (m_plinsolve->BackSolve(x, tmp) == false)
	{
		throw LinearSolverFailed();
	}

	// loop again over all update vectors
	for (int i = 0; i<nups; ++i)
	{
		int n = (n0 + i) % m_max_buf_size;

		double* vi = m_V[n];
		double* wi = m_W[n];

		double vr = 0;
		for (int j = 0; j<m_neq; ++j) vr += vi[j] * x[j];

		for (int j = 0; j<m_neq; ++j) x[j] += wi[j] * vr;
	}
}
