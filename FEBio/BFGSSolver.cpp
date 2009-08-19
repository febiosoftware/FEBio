#include "stdafx.h"
#include "BFGSSolver.h"

//-----------------------------------------------------------------------------
// NonLinearSystem
//-----------------------------------------------------------------------------

NonLinearSystem::NonLinearSystem(int neq, int ndof)
{
	m_neq = neq;
	m_ndof = ndof;
}

NonLinearSystem::~NonLinearSystem()
{

}

//-----------------------------------------------------------------------------
// BFGSSolver
//-----------------------------------------------------------------------------

BFGSSolver::BFGSSolver()
{
	m_pS = 0;
	m_pK = 0;
	m_maxups = 20;
	m_nups = 0;
}

BFGSSolver::~BFGSSolver()
{
}

//-----------------------------------------------------------------------------
//! The solve function calculates the solution to the nonlinear equation system
//! S. On input, x is the initial guess for the solution, where on output it will
//! be the solution of the equations.

bool BFGSSolver::Solve(vector<double>& x, NonLinearSystem& S)
{
	// keep a pointer to the nonlinear system of equations
	m_pS = &S;

	// get the nr of equations and nr of degrees of freedom
	int neq = S.equations();
	int ndof = S.dofs();

	// assert that the neq and ndof are the same
	assert(neq == ndof);

	// TODO: setup a linear solver
	m_pls = 0;

	// allocate a sparse matrix
	m_pK = m_pls->GetMatrix(0);

	// set the solution vector
	// it is assumed that the input parameter x an intial guess
	// contains for the solution.
	m_x = x;

	// allocate storage for the search direction
	vector<double> u(ndof);

	// calculate the initial function evaluation and its jacobian
	S.eval(x, m_F0, m_pK);

	// convergence flag
	bool bconv = false;

	// line search factor
	double s = 1.0;

	// reset BFGS update counter
	m_nups = 0;

	// repeat until converged
	do
	{
		// find a new search direction
		findsd();

		// do a line search
		linesearch();

		// update the solution
		update();
		
		// check convergence
		bconv = converged();
	}
	while (!bconv);

	return bconv;
}

//-----------------------------------------------------------------------------
//! Find a new search direction.

void BFGSSolver::findsd()
{
	// first, we see if we have done any
	// BFGS updates
	if (m_nups == 0)
	{
		// if not, we recalculate the jacobian
		// then we refactor it and do a backsolve.
		m_pS->eval(m_x, m_F0, m_pK);
		m_pls->Factor(*m_pK);
		m_pls->Solve(*m_pK, m_u, m_F0);
	}
	else
	{
		// else, we use the BFGS update vectors
		// to solve for the new search direction
		int i, j;
		double *vi, *wi, vr, wr;

		// get the nr of equations
		int neq = m_pS->dofs();

		// create temporary storage
		static vector<double> tmp;
		tmp = m_F0;

		// loop over all update vectors
		for (i=m_nups-1; i>=0; --i)
		{
			vi = m_V[i];
			wi = m_W[i];

			wr = 0;
			for (j=0; j<neq; j++) wr += wi[j]*tmp[j];
			for (j=0; j<neq; j++) tmp[j] += vi[j]*wr;
		}

		// perform a backsubstitution
		m_pls->Solve(*m_pK, m_u, tmp);

		// loop again over all update vectors
		for (i=0; i<m_nups; ++i)
		{
			vi = m_V[i];
			wi = m_W[i];

			vr = 0;
			for (j=0; j<neq; ++j) vr += vi[j]*m_u[j];
			for (j=0; j<neq; ++j) m_u[j] += wi[j]*vr;
		}
	}
}

//-----------------------------------------------------------------------------
//! do a line search

void BFGSSolver::linesearch()
{
	m_ls = 1.0;
	m_xt = m_x + m_u*m_ls;
	m_pS->eval(m_xt, m_F1);
}

//-----------------------------------------------------------------------------
//! Check convergence of the solution.

bool BFGSSolver::converged()
{
	return true;
}

//-----------------------------------------------------------------------------
//! Update the solution vector with the last increment and line search factor.
//! This also updates the BFGS vectors for the next iteration.

void BFGSSolver::update()
{
	// update the solution
	int N = m_pS->dofs();
	for (int i=0; i<N; ++i) m_x[i] += m_ls*m_u[i];

	// update the BFGS update vectors
	int i;
	double dg, dh,dgi, c, r;
	double *vn, *wn;

	int neq = m_pS->dofs();

	// calculate the BFGS update vectors
	for (i=0; i<neq; ++i)	
	{
		m_D[i] = m_ls*m_u[i];
		m_G[i] = m_F0[i] - m_F1[i];
		m_H[i] = m_F0[i]*m_ls;
	}

	dg = m_D*m_G;
	dh = m_D*m_H;
	dgi = 1.0 / dg;
	r = dg/dh;

	// check to see if this is still a pos definite update
//	if (r <= 0) 
//	{
//		return false;
//	}

	// calculate the condition number
//	c = sqrt(r);
	c = sqrt(fabs(r));

	// make sure c is less than the the maximum.
	if (c < m_cmax)
	{
		vn = m_V[m_nups];
		wn = m_W[m_nups];

		// TODO: There might be a bug here. Check signs!
		for (i=0; i<neq; ++i)	
		{
			vn[i] = -m_H[i]*c - m_G[i];
			wn[i] = m_D[i]*dgi;
		}

		// increase the BFGS update counter
		++m_nups;
	}
}
