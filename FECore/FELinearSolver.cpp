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
#include "FELinearSolver.h"
#include "FEModel.h"
#include "LinearSolver.h"
#include "FEGlobalMatrix.h"
#include "log.h"
#include "FENodeReorder.h"
#include "FELinearSystem.h"
#include "FEBoundaryCondition.h"
#include "FEGlobalVector.h"
#include "FEDomain.h"
#include "FENodalLoad.h"
#include "FESurfaceLoad.h"
#include "FEBodyLoad.h"
#include "DumpStream.h"
#include "FELinearConstraintManager.h"

//-----------------------------------------------------------------------------
//! constructor
FELinearSolver::FELinearSolver(FEModel* pfem) : FESolver(pfem)
{
	m_pls = 0;
	m_pK = 0;
	m_neq = 0;
	m_breform = true;
}

//-----------------------------------------------------------------------------
void FELinearSolver::Clean()
{
	if (m_pls) m_pls->Destroy();
}

//-----------------------------------------------------------------------------
//! This function sets the degrees of freedom that will be used by this solver.
//! This is used in the InitEquations method that initializes the equation numbers
//! and the Update method which maps the solution of the linear system to the nodal
//! data.
void FELinearSolver::SetDOF(vector<int>& dof)
{
	m_dof = dof;
}

//-----------------------------------------------------------------------------
int FELinearSolver::NumberOfEquations() const 
{
	return m_neq;
}

//-----------------------------------------------------------------------------
//! Get the linear solver
LinearSolver* FELinearSolver::GetLinearSolver()
{
	return m_pls;
}

//-----------------------------------------------------------------------------
bool FELinearSolver::Init()
{
	if (FESolver::Init() == false) return false;

	// Now that we have determined the equation numbers we can continue
	// with creating the stiffness matrix. First we select the linear solver
	// The stiffness matrix is created in CreateStiffness
	// Note that if a particular solver was requested in the input file
	// then the solver might already be allocated. That's way we need to check it.
	if (m_pls == 0)
	{
		FEModel* fem = GetFEModel();
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		m_pls = fecore.CreateDefaultLinearSolver(fem);
		if (m_pls == 0)
		{
			feLogError("Unknown solver type selected\n");
			return false;
		}

		if (m_part.empty() == false)
		{
			m_pls->SetPartitions(m_part);
		}
	}

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the solver allocate the correct type of matrix format
	Matrix_Type mtype = MatrixType();
	SparseMatrix* pS = m_pls->CreateSparseMatrix(mtype);
	if ((pS == 0) && (m_msymm == REAL_SYMMETRIC))
	{
		// oh, oh, something went wrong. It's probably because the user requested a symmetric matrix for a 
		// solver that wants a non-symmetric. If so, let's force a non-symmetric format.
		pS = m_pls->CreateSparseMatrix(REAL_UNSYMMETRIC);

		if (pS)
		{
			// Problem solved! Let's inform the user.
			m_msymm = REAL_UNSYMMETRIC;
			feLogWarning("The matrix format was changed to non-symmetric since the selected linear solver does not support a symmetric format. \n");
		}
	}

	if (pS == 0)
	{
		feLogError("The selected linear solver does not support the requested matrix format.\nPlease select a different linear solver.\n");
		return false;
	}

	// clean up the stiffness matrix if we have one
	if (m_pK) delete m_pK; m_pK = 0;

	// Create the stiffness matrix.
	// Note that this does not construct the stiffness matrix. This
	// is done later in the StiffnessMatrix routine.
	m_pK = new FEGlobalMatrix(pS);
	if (m_pK == 0)
	{
		feLogError("Failed allocating stiffness matrix\n\n");
		return false;
	}

	// Set the matrix formation flag
	m_breform = true;

	// get number of equations
	int neq = m_neq;

	// allocate data structures
	m_R.resize(neq);
	m_u.resize(neq);

	return true;
}

//-----------------------------------------------------------------------------
//! Solve an analysis step
bool FELinearSolver::SolveStep()
{
	// Make sure we have a linear solver and a stiffness matrix
	if (m_pls == 0) return false;
	if (m_pK == 0) return false;

	// reset counters
	m_niter = 0;
	m_nrhs = 0;
	m_nref = 0;
	m_ntotref = 0;

	FEModel& fem = *GetFEModel();

	// Set up the prescribed dof vector
	// The stiffness matrix assembler uses this to update the RHS vector
	// for prescribed dofs.
	zero(m_u);
	int nbc = fem.BoundaryConditions();
	for (int i=0; i<nbc; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.PrepStep(m_u, false);
	}

	// build the right-hand side
	// (Is done by the derived class)
	zero(m_R);
	vector<double> F(m_neq);
	FEGlobalVector rhs(fem, m_R, F);
	{
		TRACK_TIME(TimerID::Timer_Residual);
		ForceVector(rhs);
	}

	// increase RHS counter
	m_nrhs++;

	// build the stiffness matrix
	ReformStiffness();

	// solve the equations
	vector<double> u(m_neq);
	{
		TRACK_TIME(TimerID::Timer_LinSolve);
		if (m_pls->BackSolve(u, m_R) == false)
			throw LinearSolverFailed();
	}

	// print norms
	double rnorm = sqrt(m_R*m_R);
	double unorm = sqrt(u*u);
	feLog("\trhs norm      :  %15le\n", rnorm);
	feLog("\tsolution norm :  %15le\n", unorm);

	// the prescribed values are stored in m_u, so let's add it
	u += m_u;

	// update solution
	Update(u);

	// increase iteration count
	m_niter++;

	return true;
}

//-----------------------------------------------------------------------------
bool FELinearSolver::ReformStiffness()
{
	// recalculate the shape of the stiffness matrix if necessary
	if (m_breform)
	{
		TRACK_TIME(TimerID::Timer_Reform);

		if (!CreateStiffness()) return false;
		
		// since it's not likely that the matrix form changes
		// in a linear analysis, we don't recalculate the form
		m_breform = false;
	}

	// Make sure it is all set to zero
	m_pK->Zero();

	// calculate the stiffness matrix
	// (This is done by the derived class)
	{
		TRACK_TIME(TimerID::Timer_Stiffness);

		FELinearSystem K(this, *m_pK, m_R, m_u, (m_msymm == REAL_SYMMETRIC));
		if (!StiffnessMatrix(K)) return false;

		// do call back
		FEModel& fem = *GetFEModel();
		fem.DoCallback(CB_MATRIX_REFORM);
	}

	// factorize the stiffness matrix
	{
		TRACK_TIME(TimerID::Timer_LinSolve);
		m_pls->Factor();
	}

	// increase total nr of reformations
	m_nref++;
	m_ntotref++;

	return true;
}

//-----------------------------------------------------------------------------
bool FELinearSolver::CreateStiffness()
{
	// clean up the solver
	if (m_pK->NonZeroes()) m_pls->Destroy();

	// clean up the stiffness matrix
	m_pK->Clear();

	// create the stiffness matrix
	feLog("===== reforming stiffness matrix:\n");
	if (m_pK->Create(GetFEModel(), m_neq, true) == false) 
	{
		feLogError("An error occured while building the stiffness matrix\n\n");
		return false;
	}
	else
	{
		// output some information about the direct linear solver
		int neq = m_pK->Rows();
		int nnz = m_pK->NonZeroes();
		feLog("\tNr of equations ........................... : %d\n", neq);
		feLog("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
	}

	// Do the preprocessing of the solver
	{
		TRACK_TIME(TimerID::Timer_LinSolve);
		if (!m_pls->PreProcess()) throw FatalError();
	}

	// done!
	return true;
}

//-----------------------------------------------------------------------------
FEGlobalMatrix* FELinearSolver::GetStiffnessMatrix()
{
	return m_pK;
}

//-----------------------------------------------------------------------------
//! get the RHS
std::vector<double>	FELinearSolver::GetLoadVector()
{
	return m_R;
}

//-----------------------------------------------------------------------------
void FELinearSolver::Serialize(DumpStream& ar)
{
	FESolver::Serialize(ar);

	ar & m_breform & m_neq;

	if (ar.IsSaving() == false)
	{
		// We need to rebuild the stiffness matrix at some point.
		// Currently this is done during Activation, but we don't
		// call FEAnalysis::Activate after restart so for now,
		// I'll just do it here.
		// TODO: Find a better way.
		FELinearSolver::Init();
	}
}

//-----------------------------------------------------------------------------
//! Evaluate the right-hand side "force" vector
void FELinearSolver::ForceVector(FEGlobalVector& R)
{
	// get the modal and mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FETimeInfo& tp = fem.GetTime();

	// loop over model loads
	int nml = fem.ModelLoads();
	for (int i = 0; i < nml; ++i)
	{
		FEModelLoad& ml = *fem.ModelLoad(i);
		if (ml.IsActive()) ml.LoadVector(R);
	}
}

//-----------------------------------------------------------------------------
//! Evaluate the stiffness matrix
bool FELinearSolver::StiffnessMatrix(FELinearSystem& K) 
{ 
	return false; 
}

//-----------------------------------------------------------------------------
// This function copies the solution back to the nodal variables
// and class the FEDomain::Update to give domains a chance to update
// their local data.
// TODO: Instead of calling Update on all domains, perhaps I should introduce
//       a mechanism for solvers only update the domains that are relevant.
void FELinearSolver::Update(vector<double>& u)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	const FETimeInfo& tp = fem.GetTime();

	// update nodal variables
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<m_dof.size(); ++j)
		{
			int n = node.m_ID[m_dof[j]];
			if (n >= 0) node.set(m_dof[j], u[n]);
			else if (-n-2 >= 0) node.set(m_dof[j], u[-n-2]);
		}
	}

	// make sure linear constraints are satisfied
	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
	if (LCM.LinearConstraints())
	{
		LCM.Update();
	}

	// update the domains
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.Update(tp);
	}

	fem.IncrementUpdateCounter();
}
