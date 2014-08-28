#include "stdafx.h"
#include "BFGSSolver.h"
#include "FESolver.h"

//-----------------------------------------------------------------------------
// BFGSSolver
//-----------------------------------------------------------------------------

BFGSSolver::BFGSSolver()
{
	m_maxups = 10;
	m_maxref = 15;
	m_cmax   = 1e5;

	// default line search parameters
	m_LSmin  = 0.01;
	m_LStol = 0.9;
	m_LSiter = 5;

	// pointer to linear solver
	m_plinsolve = 0;

	// pointer to non-linear system
	m_pNLS = 0;
}

//-----------------------------------------------------------------------------
// Initialization method
void BFGSSolver::Init(int neq, FESolver* pNLS, LinearSolver* pls)
{
	// allocate storage for BFGS update vectors
	m_V.resize(m_maxups, neq);
	m_W.resize(m_maxups, neq);

	m_D.resize(neq);
	m_G.resize(neq);
	m_H.resize(neq);

	m_ui.assign(neq, 0);
	m_R0.assign(neq, 0);
	m_R1.assign(neq, 0);

	m_neq = neq;
	m_nups = 0;

	m_plinsolve = pls;

	assert(pNLS);
	m_pNLS = pNLS;
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
	int i;
	double dg, dh,dgi, c, r;
	double *vn, *wn;

	int neq = ui.size();

	// calculate the BFGS update vectors
	for (i=0; i<neq; ++i)	
	{
		m_D[i] = s*ui[i];
		m_G[i] = R0[i] - R1[i];
		m_H[i] = R0[i]*s;
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
	if (c > m_cmax) return false;

	vn = m_V[m_nups];
	wn = m_W[m_nups];

	// TODO: There might be a bug here. Check signs!
	for (i=0; i<neq; ++i)	
	{
		vn[i] = -m_H[i]*c - m_G[i];
		wn[i] = m_D[i]*dgi;
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
	int i, j;
	double *vi, *wi, vr, wr;

	// get the nr of equations
	int neq = x.size();

	// make sure we need to do work
	if (neq==0) return;

	// create temporary storage
	static vector<double> tmp;
	tmp = b;

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
	m_plinsolve->BackSolve(x, tmp);

	// loop again over all update vectors
	for (i=0; i<m_nups; ++i)
	{
		vi = m_V[i];
		wi = m_W[i];

		vr = 0;
		for (j=0; j<neq; ++j) vr += vi[j]*x[j];

		for (j=0; j<neq; ++j) x[j] += wi[j]*vr;
	}
}

//-----------------------------------------------------------------------------
//! Performs a linesearch on a NR iteration
//! The description of this method can be found in:
//!    "Nonlinear Continuum Mechanics for Finite Element Analysis", 	Bonet & Wood.
//
//! \todo Find a different way to update the deformation based on the ls.
//! For instance, define a di so that ui = s*di. Also, define the 
//! position of the nodes at the previous iteration.

double BFGSSolver::LineSearch(double s)
{
	double smin = s;

	double a, A, B, D;
	double r0, r1, r;

	// max nr of line search iterations
	int nmax = m_LSiter;
	int n = 0;

	// initial energy
	r0 = m_ui*m_R0;

	double rmin = fabs(r0);

/*	// get the logfile
	Logfile& log = GetLogfile();

	if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
	{
		log.printf("\nEntering line search\n");
		log.printf("       STEPSIZE     INITIAL        CURRENT         REQUIRED\n");
	}
*/

	// ul = ls*ui
	vector<double> ul(m_ui.size());
	do
	{
		// Update geometry
		vcopys(ul, m_ui, s);
		m_pNLS->Update(ul);

		// calculate residual at this point
		m_pNLS->Evaluate(m_R1);

		// make sure we are still in a valid range
		if (s < m_LSmin) 
		{
			// it appears that we are not converging
			// I found in the NIKE3D code that when this happens,
			// the line search step is simply set to 0.5.
			// so let's try it here too
			s = 0.5;

			// reupdate  
			vcopys(ul, m_ui, s);
			m_pNLS->Update(ul);

			// recalculate residual at this point
			m_pNLS->Evaluate(m_R1);

			// return and hope for the best
			break;
		}

		// calculate energies
		r1 = m_ui*m_R1;

		if ((n==0) || (fabs(r1) < rmin))
		{
			smin = s;
			rmin = fabs(r1);
		}

		// make sure that r1 does not happen to be really close to zero,
		// since in that case we won't find any better solution.
		if (fabs(r1) < 1.e-20) r = 0;
		else r = fabs(r1/r0);

/*		if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
		{
			log.printf("%15lf%15lg%15lg%15lg\n", s, fabs(r0), fabs(r1), fabs(r0*m_bfgs.m_LStol));
		}
*/
		if (r > m_LStol)
		{
			// calculate the line search step
			a = r0/r1;

			A = 1 + a*(s-1);
			B = a*(s*s);
			D = B*B - 4*A*B;

			// update the line step
			if (D >= 0) 
			{
				s = (B + sqrt(D))/(2*A);
				if (s < 0) s = (B - sqrt(D))/(2*A);
				if (s < 0) s = 0;
			}
			else 
			{
				s = 0.5*B/A;
			}

			++n;
		}
	}
	while ((r > m_LStol) && (n < nmax));

	if (n >= nmax)
	{
		// max nr of iterations reached.
		// we choose the line step that reached the smallest energy
		s = smin;
		vcopys(ul, m_ui, s);
		m_pNLS->Update(ul);
		m_pNLS->Evaluate(m_R1);
	}
/*
	if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
	{
		if ((r < m_bfgs.m_LStol) && (n < nmax))
		{
			log.printf("------------------------------------------------------------\n");
			log.printf("%15lf%15lg%15lg%15lg\n\n", s, fabs(r0), fabs(r1), fabs(r0*m_bfgs.m_LStol));
		}
		else if (n >= nmax)
		{
			log.printf("Line search failed: max iterations reached.\n\n");
		}
	}
*/
	return s;
}


double BFGSSolver::LineSearchCG(double s)
{
	//  s is now passed from the solver routine instead of defaulting to 1.0
	double smin = s;

	double FA, FB, FC, AA, AB, r, r0 ;
	bool failed=false;

	// max nr of line search iterations
	int nmax = m_LSiter;
	int n = 0;

	// initial energy
	FA = m_ui*m_R0;
	AA=0.0;
	r0=FA;

	double rmin = fabs(FA);

	vector<double> ul(m_ui.size());

	// so we can set AA = 0 and FA= r0
	// AB=s and we need to evaluate FB (called r1)
	// n is a count of the number of linesearch attempts

		// calculate residual at this point, reducing s if necessary
		do
		{
			// Update geometry using the initial guess s
			vcopys(ul, m_ui, s);
			failed=false;
			try
				{
					m_pNLS->Update(ul);
					m_pNLS->Evaluate(m_R1);
				}
				catch (...)
				{
//					printf("reducing s at initial evaluation");
					failed=true;
					s=0.1*s;
				}
		}
		while (failed==true);

		// calculate energies
		FB = m_ui*m_R1;
		AB=s;

		if (FB<rmin){
			rmin=FB;
			smin=s;
		}

	do
	{
		// make sure that r1 does not happen to be really close to zero,
		// since in that case we won't find any better solution.
		if (fabs(FB) < 1.e-20) r = 0;  // we've hit the converged solution and don't need to do any more
		else r = fabs(FB/r0);

		if (r > m_LStol)	// we need to search and find a better value of s
		{
            if (FB==FA) s=(AA+AB)*1000;
			else { 
				s=(AA*FB-AB*FA)/(FB-FA);
				s=min(s,100*max(AA,AB));
			}
		// calculate residual at this point, reducing s if necessary
		do
		{
			// Update geometry using the initial guess s
			vcopys(ul, m_ui, s);
			failed=false;
			try
				{
					m_pNLS->Update(ul);
					m_pNLS->Evaluate(m_R1);
				}
				catch (...)
				{
//					printf("reducing s at FC");
					failed=true;
					s=0.1*s;
				}
		}
		while ((failed==true)&&(s>m_LSmin));

			// calculate energies
			FC = m_ui*m_R1;
			r=fabs(FC/r0);
            
			if (fabs(FC)>100*min(fabs(FA),fabs(FB)))  //  it was a bad guess and we need to go back a bit
			{
				s=0.1*s;

		// calculate residual at this point, reducing s if necessary
		do
		{
			// Update geometry using the initial guess s
			vcopys(ul, m_ui, s);
			failed=false;
			try
				{
					m_pNLS->Update(ul);
					m_pNLS->Evaluate(m_R1);
				}
				catch (...)
				{
//					printf("reducing s after bad guess");
					failed=true;
					s=0.1*s;
				}
		}
		while (failed==true);

			// calculate energies
			FC = m_ui*m_R1;
			r=fabs(FC/r0);
			}
            
            if (fabs(FA)<fabs(FB)) // use the new value and the closest of the previous ones
			{
				FB=FC;
				AB=s;
			}
			else
			{
               FA=FC;
			   AA=s;
			}

			++n;
		}
	}
	while ((r > m_LStol) && (n < nmax));


	if (n >= nmax)
	{
		// max nr of iterations reached.
		// we choose the line step that reached the smallest energy
		s = smin;
		// calculate residual at this point, reducing s if necessary
		do
		{
			// Update geometry using the initial guess s
			vcopys(ul, m_ui, s);
			failed=false;
			try
				{
					m_pNLS->Update(ul);
					m_pNLS->Evaluate(m_R1);
				}
				catch (...)
				{
//					printf("reducing s after failed line search");
					failed=true;
					s=0.1*s;
				}
		}
		while (failed==true);
	}

	return s;
}
