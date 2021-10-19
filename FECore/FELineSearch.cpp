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
#include "FELineSearch.h"
#include "FENewtonSolver.h"
#include "DumpStream.h"
#include <vector>
using namespace std;

FELineSearch::FELineSearch(FENewtonSolver* pns) : m_pns(pns)
{
	m_LSmin = 0.01;
	m_LStol = 0.9;
	m_LSiter = 5;
}

// serialization
void FELineSearch::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar & m_LSmin & m_LStol & m_LSiter;
}

//! Performs a linesearch on a NR iteration
//! The description of this method can be found in:
//!    "Nonlinear Continuum Mechanics for Finite Element Analysis", 	Bonet & Wood.
//
//! \todo Find a different way to update the deformation based on the ls.
//! For instance, define a di so that ui = s*di. Also, define the 
//! position of the nodes at the previous iteration.
double FELineSearch::DoLineSearch(double s)
{
	assert(m_pns);

	double smin = s;

	double a, A, B, D;
	double r0, r1, r;

	// max nr of line search iterations
	int nmax = m_LSiter;
	int n = 0;

	// vectors
	vector<double>& ui = m_pns->m_ui;
	vector<double>& R0 = m_pns->m_R0;
	vector<double>& R1 = m_pns->m_R1;

	// initial energy
	r0 = ui*R0;

	double rmin = fabs(r0);

	FENewtonStrategy* ns = m_pns->m_qnstrategy;

	// ul = ls*ui
	vector<double> ul(ui.size());
	do
	{
		// Update geometry
		vcopys(ul, ui, s);
		m_pns->Update(ul);

		// calculate residual at this point
		ns->Residual(R1, false);

		// make sure we are still in a valid range
		if (s < m_LSmin)
		{
			// it appears that we are not converging
			// I found in the NIKE3D code that when this happens,
			// the line search step is simply set to 0.5.
			// so let's try it here too
			s = 0.5;

			// reupdate  
			vcopys(ul, ui, s);
			m_pns->Update(ul);

			// recalculate residual at this point
			ns->Residual(R1, false);

			// return and hope for the best
			break;
		}

		// calculate energies
		r1 = ui*R1;

		if ((n == 0) || (fabs(r1) < rmin))
		{
			smin = s;
			rmin = fabs(r1);
		}

		// make sure that r1 does not happen to be really close to zero,
		// since in that case we won't find any better solution.
		if (fabs(r1) < 1.e-17) r = 0;
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
		vcopys(ul, ui, s);
		m_pns->Update(ul);
		ns->Residual(R1, false);
	}
	return s;
}
