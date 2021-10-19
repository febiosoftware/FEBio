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
#include "FESolidLinearSystem.h"
#include "FESolidSolver.h"
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FEModel.h>

FESolidLinearSystem::FESolidLinearSystem(FESolver* solver, FERigidSolver* rigidSolver, FEGlobalMatrix& K, std::vector<double>& F, std::vector<double>& u, bool bsymm, double alpha, int nreq) : FELinearSystem(solver, K, F, u, bsymm)
{
	m_rigidSolver = rigidSolver;
	m_alpha = alpha;
	m_nreq = nreq;
	m_stiffnessScale = 1.0;
}

// scale factor for stiffness matrix
void FESolidLinearSystem::StiffnessAssemblyScaleFactor(double a)
{
	m_stiffnessScale = a;
}

void FESolidLinearSystem::Assemble(const FEElementMatrix& ke)
{
	// Rigid joints require a different assembly approach in that we can do 
	// a direct assembly as defined by the base class. 
	// Currently, we assume that if the node list of the element matrix is not
	// defined, then we are dealing with rigid joints.
	if (ke.Nodes().empty())
	{
		FELinearSystem::Assemble(ke);
	}
	else
	{
		// assemble into global stiffness matrix
		if (m_stiffnessScale == 1.0)
		{
			m_K.Assemble(ke);
		}
		else
		{
			// NOTE: What is doing here?
			FEElementMatrix kes(ke, m_stiffnessScale);
			m_K.Assemble(kes);
		}

		// get the vector that stores the prescribed BC values
		vector<double>& ui = m_u;

		// adjust for linear constraints
		FEModel* fem = m_solver->GetFEModel();
		FELinearConstraintManager& LCM = fem->GetLinearConstraintManager();
		if (LCM.LinearConstraints() > 0)
		{
			#pragma omp critical 
			LCM.AssembleStiffness(m_K, m_F, m_u, ke.Nodes(), ke.RowIndices(), ke.ColumnsIndices(), ke);
		}

		// adjust stiffness matrix for prescribed degrees of freedom
		// NOTE: I had to comment this if statement out since otherwise
		//       poroelastic DOF's that are set as free-draining in the
		//       sliding2 contact code are skipt and zeroes will appear
		//       on the diagonal of the stiffness matrix.
		//	if (m_fem.m_DC.size() > 0)
		{
			SparseMatrix& K = m_K;

			int N = ke.rows();

			// loop over columns
			const vector<int>& elmi = ke.RowIndices();
			const vector<int>& elmj = ke.ColumnsIndices();
			for (int j = 0; j < N; ++j)
			{
				int J = -elmj[j] - 2;
				if ((J >= 0) && (J < m_nreq))
				{
					// dof j is a prescribed degree of freedom

					// loop over rows
					for (int i = 0; i < N; ++i)
					{
						int I = elmi[i];
						if (I >= 0)
						{
							// dof i is not a prescribed degree of freedom
							#pragma omp atomic
							m_F[I] -= ke[i][j] * ui[J];
						}
					}

					// set the diagonal element of K to 1
					K.set(J, J, 1);
				}
			}
		}

		// see if there are any rigid body dofs here
		#pragma omp critical 
		m_rigidSolver->RigidStiffness(m_K, m_u, m_F, ke, m_alpha);
	}
}
