#include "stdafx.h"
#include "FENewtonSolver.h"

//-----------------------------------------------------------------------------
FENewtonSolver::FENewtonSolver(FEModel* pfem) : FESolver(pfem)
{
	// default line search parameters
	m_LSmin = 0.01;
	m_LStol = 0.9;
	m_LSiter = 5;
}

//-----------------------------------------------------------------------------
//! Performs a linesearch on a NR iteration
//! The description of this method can be found in:
//!    "Nonlinear Continuum Mechanics for Finite Element Analysis", 	Bonet & Wood.
//
//! \todo Find a different way to update the deformation based on the ls.
//! For instance, define a di so that ui = s*di. Also, define the 
//! position of the nodes at the previous iteration.

double FENewtonSolver::LineSearch(double s)
{
	double smin = s;

	double a, A, B, D;
	double r0, r1, r;

	// max nr of line search iterations
	int nmax = m_LSiter;
	int n = 0;

	// initial energy
	r0 = m_bfgs.m_ui*m_bfgs.m_R0;

	double rmin = fabs(r0);

	// ul = ls*ui
	vector<double> ul(m_bfgs.m_ui.size());
	do
	{
		// Update geometry
		vcopys(ul, m_bfgs.m_ui, s);
		Update(ul);

		// calculate residual at this point
		Evaluate(m_bfgs.m_R1);

		// make sure we are still in a valid range
		if (s < m_LSmin)
		{
			// it appears that we are not converging
			// I found in the NIKE3D code that when this happens,
			// the line search step is simply set to 0.5.
			// so let's try it here too
			s = 0.5;

			// reupdate  
			vcopys(ul, m_bfgs.m_ui, s);
			Update(ul);

			// recalculate residual at this point
			Evaluate(m_bfgs.m_R1);

			// return and hope for the best
			break;
		}

		// calculate energies
		r1 = m_bfgs.m_ui*m_bfgs.m_R1;

		if ((n == 0) || (fabs(r1) < rmin))
		{
			smin = s;
			rmin = fabs(r1);
		}

		// make sure that r1 does not happen to be really close to zero,
		// since in that case we won't find any better solution.
		if (fabs(r1) < 1.e-20) r = 0;
		else r = fabs(r1 / r0);

		if (r > m_LStol)
		{
			// calculate the line search step
			a = r0 / r1;

			A = 1 + a*(s - 1);
			B = a*(s*s);
			D = B*B - 4 * A*B;

			// update the line step
			if (D >= 0)
			{
				s = (B + sqrt(D)) / (2 * A);
				if (s < 0) s = (B - sqrt(D)) / (2 * A);
				if (s < 0) s = 0;
			}
			else
			{
				s = 0.5*B / A;
			}

			++n;
		}
	} while ((r > m_LStol) && (n < nmax));

	if (n >= nmax)
	{
		// max nr of iterations reached.
		// we choose the line step that reached the smallest energy
		s = smin;
		vcopys(ul, m_bfgs.m_ui, s);
		Update(ul);
		Evaluate(m_bfgs.m_R1);
	}
	return s;
}
