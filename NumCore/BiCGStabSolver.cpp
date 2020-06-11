/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "BiCGStabSolver.h"
#include "CompactUnSymmMatrix.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(BiCGStabSolver, IterativeLinearSolver)
	ADD_PARAMETER(m_print_level, "print_level");
	ADD_PARAMETER(m_tol, "tol");
	ADD_PARAMETER(m_maxiter, "max_iter");
	ADD_PARAMETER(m_fail_max_iter, "fail_max_iters");
	ADD_PROPERTY(m_P, "pc_left");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
BiCGStabSolver::BiCGStabSolver(FEModel* fem) : IterativeLinearSolver(fem), m_pA(0), m_P(0)
{
	m_maxiter = 0;
	m_tol = 1e-5;
	m_abstol = 0.0;
	m_print_level = 0;
	m_fail_max_iter = true;
}

//-----------------------------------------------------------------------------
SparseMatrix* BiCGStabSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	// let the preconditioner decide
	m_pA = nullptr;
	if (m_P)
	{
		m_P->SetPartitions(m_part);
		m_pA = m_P->CreateSparseMatrix(ntype);
	}
	else
	{
		if (ntype == REAL_SYMMETRIC) m_pA = new CompactSymmMatrix;
		else m_pA = new CRSSparseMatrix(1);
	}
	return m_pA;
}

//-----------------------------------------------------------------------------
bool BiCGStabSolver::SetSparseMatrix(SparseMatrix* A)
{
	m_pA = A;
	return (m_pA != 0);
}

//-----------------------------------------------------------------------------
void BiCGStabSolver::SetLeftPreconditioner(LinearSolver* P)
{
	m_P = P;
}

//-----------------------------------------------------------------------------
LinearSolver* BiCGStabSolver::GetLeftPreconditioner()
{
	return m_P;
}

//-----------------------------------------------------------------------------
bool BiCGStabSolver::HasPreconditioner() const
{
	return (m_P != nullptr);
}

//-----------------------------------------------------------------------------
bool BiCGStabSolver::PreProcess()
{
	return true;
}

//-----------------------------------------------------------------------------
bool BiCGStabSolver::Factor()
{
	if (m_pA == 0) return false;
	if (m_P)
	{
		if (m_P->PreProcess() == false) return false;
		if (m_P->Factor() == false) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool BiCGStabSolver::BackSolve(double* x, double* b)
{
	if (m_pA == nullptr) return false;

	SparseMatrix& A = *m_pA;
	int neq = A.Rows();

	// assume initial guess is zero
	for (int i = 0; i < neq; ++i) x[i] = 0.0;

	// calculate initial norm
	// r0 = b - A*x0
	vector<double> r_i(neq); double norm0 = 0.0, normi = 0.0;
	for (int j = 0; j < neq; ++j)
	{
		r_i[j] = b[j];
		norm0 += r_i[j] * r_i[j];
	}
	norm0 = sqrt(norm0);

	// if the norm is zero, there is nothing to do
	if (norm0 == 0.0) return true;

//	Choose an arbitrary vector rt such that(rt, r0) != 0, e.g., rt = r0
	vector<double> rt(r_i);

	// initialize some stuff
	double rho_p = 1, alpha = 1, w_p = 1;
	vector<double> v_p(neq, 0.0), p_p(neq, 0.0), p_i(neq), y(neq, 0.0), h(neq), s(neq), z(neq), t(neq), q(neq);

	int iter = 0;
	bool converged = false;
	do
	{
		double rho_i = rt*r_i;

		double beta = (rho_i / rho_p)*(alpha / w_p);

		for (int j = 0; j < neq; ++j) p_i[j] = r_i[j] + beta*(p_p[j] - w_p*v_p[j]);

		// apply preconditioner
		if (m_P)
		{
			m_P->BackSolve(y, p_i);
		}
		else y = p_i;

		A.mult_vector(&y[0], &v_p[0]);

		alpha = rho_i / (rt*v_p);

		for (int j = 0; j < neq; ++j) h[j] = x[j] + alpha*y[j];

		for (int j = 0; j < neq; ++j) s[j] = r_i[j] - alpha*v_p[j];
//		If h is accurate enough then xi = h and quit

		if (m_P)
		{
			m_P->BackSolve(z, s);
		}
		else z = s;

		A.mult_vector(&z[0], &t[0]);

		if (m_P)
		{
			m_P->BackSolve(q, t);
		}
		else q = t;

		w_p = (q*z) / (q*q);

		for (int j = 0; j < neq; ++j) x[j] = h[j] + w_p*z[j];
		
		normi = 0.0;
		for (int j = 0; j < neq; ++j)
		{
			r_i[j] = s[j] - w_p*t[j];
			normi += r_i[j] * r_i[j];
		}
		normi = sqrt(normi);

		// see if we have converged
		double tol = norm0*m_tol + m_abstol;
		if (normi <= tol) converged = true;
		else
		{
			// prepare for next iteration
			rho_p = rho_i;
			p_p = p_i;
		}

		// check max iterations
		iter++;
		if (iter > m_maxiter) break;

		if (m_print_level > 1)
		{
			feLog("%d:%lg, %lg\n", iter, normi, tol);
		}
	}
	while (!converged);

	if (m_print_level == 1)
	{
		feLog("%d:%lg, %lg\n", iter, normi, norm0);
	}

	return (m_fail_max_iter ? converged : true);
}

//-----------------------------------------------------------------------------
void BiCGStabSolver::Destroy()
{
}
