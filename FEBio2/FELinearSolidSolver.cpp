#include "stdafx.h"
#include "FELinearSolidSolver.h"
#include "FEBioLib/FELinearSolidDomain.h"
#include "FECore/FENodeReorder.h"
#include "FEBioLib/FEPressureLoad.h"
#include "FEBioLib/FERigid.h"
#include "FEBioLib/log.h"

//-----------------------------------------------------------------------------
//! Class constructor
FELinearSolidSolver::FELinearSolidSolver(FEModel& fem) : FESolver(fem)
{
	m_Dtol = 1e-9;
}

//-----------------------------------------------------------------------------
//! Initialization
//! This function will be called before each analysis step
bool FELinearSolidSolver::Init()
{
	// initialize base class
	if (FESolver::Init() == false) return false;

	// number of equations
	int neq = m_neq;

	// allocate data structures
	m_u.resize(neq);
	m_R.resize(neq);
	m_d.resize(neq);

	// Identify the linear elastic domains
	// TODO: I want this to be done automatically
	//       e.g. while the input file is being read
	// TODO: Do this in the analysis class
	FEMesh& mesh = m_fem.m_mesh;
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	pstep->ClearDomains();
	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FELinearElasticDomain* pd = dynamic_cast<FELinearElasticDomain*>(&mesh.Domain(nd));
		if (pd) pstep->AddDomain(nd);
	}
	assert(pstep->Domains() != 0);

	return true;
}


//-----------------------------------------------------------------------------
//!	This function initializes the equation system.
//! It is assumed that all free dofs up until now have been given an ID >= 0
//! and the fixed or rigid dofs an ID < 0.
//! After this operation the nodal ID array will contain the equation
//! number assigned to the corresponding degree of freedom. To distinguish
//! between free or unconstrained dofs and constrained ones the following rules
//! apply to the ID array:
//!
//!           /
//!          |  >=  0 --> dof j of node i is a free dof
//! ID[i][j] <  == -1 --> dof j of node i is a fixed (no equation assigned too)
//!          |  <  -1 --> dof j of node i is constrained and has equation nr = -ID[i][j]-2
//!           \
//!
bool FELinearSolidSolver::InitEquations()
{
	FEMesh& mesh = m_fem.m_mesh;

	// initialize nr of equations
	int neq = 0;

	// see if we need to optimize the bandwidth
	if (m_fem.m_bwopt)
	{
		// reorder the node numbers
		vector<int> P(mesh.Nodes());
		FENodeReorder mod;
		mod.Apply(mesh, P);

		// set the equation numbers
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(P[i]);
			if (node.m_ID[DOF_X] >= 0) node.m_ID[DOF_X] = neq++;
			if (node.m_ID[DOF_Y] >= 0) node.m_ID[DOF_Y] = neq++;
			if (node.m_ID[DOF_Z] >= 0) node.m_ID[DOF_Z] = neq++;
		}
	}
	else
	{
		// give all free dofs an equation number
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			if (node.m_ID[DOF_X] >= 0) node.m_ID[DOF_X] = neq++;
			if (node.m_ID[DOF_Y] >= 0) node.m_ID[DOF_Y] = neq++;
			if (node.m_ID[DOF_Z] >= 0) node.m_ID[DOF_Z] = neq++;
		}
	}

	// store the number of equations
	m_neq = neq;

	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
//! This function will solve the FE problem
//!
bool FELinearSolidSolver::SolveStep(double time)
{
	// prepare step
	FEMaterialPoint::dt = m_fem.GetCurrentStep()->m_dt;
	FEMaterialPoint::time = m_fem.m_ftime;

	FEMesh& mesh = m_fem.m_mesh;
	for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).InitElements();

	// set-up the prescribed displacements
	zero(m_d);
	vector<double> DT(m_d), DI(m_d);
	for (size_t i=0; i<m_fem.m_DC.size(); ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;

			double D = s*m_fem.GetLoadCurve(lc)->Value();

			FENode& node = m_fem.m_mesh.Node(n);

			if (bc == DOF_X) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_neq) { DT[I] = D; DI[I] = D - (node.m_rt.x - node.m_r0.x); }}
			if (bc == DOF_Y) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_neq) { DT[I] = D; DI[I] = D - (node.m_rt.y - node.m_r0.y); }}
			if (bc == DOF_Z) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_neq) { DT[I] = D; DI[I] = D - (node.m_rt.z - node.m_r0.z); }}
		}
	}

	// start Newton-loop
	bool bconv = false;
	int N = m_u.size();
	vector<double> du(N), Du(N), U(m_u), D(m_d); zero(du); zero(Du);
	const int NMAX = 10;
	int n = 0;
	m_niter = 0;
	m_nrhs = 0;
	m_ntotref = 0;
	do
	{
		// build the residual
		Residual();

		// build the stiffness matrix
		m_d = DI;
		ReformStiffness();

		// solve the equations
		m_plinsolve->BackSolve(du, m_R);

		// update solution
		// NOTE: m_u is not being used in Update!
		Du += du;
		m_u = U + Du;

		m_d = DT;
		Update(m_u);
		zero(DI);

		// check convergence
		double normu = fabs(du*du);
		double normR = fabs(m_R*m_R);

		clog.printf("normu = %lg\n", normu);

		if (normu < m_Dtol) bconv = true;
		n++;

		m_niter++;
		m_nrhs++;
	}
	while (!bconv && (n < NMAX));

	return bconv;
}

//-----------------------------------------------------------------------------
//! update solution
void FELinearSolidSolver::Update(vector<double>& u)
{
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	FEMesh& mesh = m_fem.m_mesh;

	// update nodal positions
	int n;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		n = node.m_ID[DOF_X]; if (n >= 0) node.m_rt.x = node.m_r0.x + u[n]; else if (-n-2 >= 0) node.m_rt.x = node.m_r0.x + m_d[-n-2];
		n = node.m_ID[DOF_Y]; if (n >= 0) node.m_rt.y = node.m_r0.y + u[n]; else if (-n-2 >= 0) node.m_rt.y = node.m_r0.y + m_d[-n-2];
		n = node.m_ID[DOF_Z]; if (n >= 0) node.m_rt.z = node.m_r0.z + u[n]; else if (-n-2 >= 0) node.m_rt.z = node.m_r0.z + m_d[-n-2];
	}

	// update the stresses on all domains
	for (int i=0; i<pstep->Domains(); ++i)
	{
		FELinearElasticDomain& d = dynamic_cast<FELinearElasticDomain&>(*pstep->Domain(i));
		d.UpdateStresses(m_fem);
	}

	// dump all states to the plot file
	// when requested
	if (m_fem.GetCurrentStep()->m_nplot == FE_PLOT_MINOR_ITRS) m_fem.Write();
}

//-----------------------------------------------------------------------------
//! Calculate the residual
void FELinearSolidSolver::Residual()
{
	zero(m_R);

	FEMesh& mesh = m_fem.m_mesh;
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// loop over nodal forces
	int ncnf = m_fem.m_FC.size();
	for (int i=0; i<ncnf; ++i)
	{
		FENodalForce& fc = *m_fem.m_FC[i];
		if (fc.IsActive())
		{
			int id	 = fc.node;	// node ID
			int bc   = fc.bc;	// direction of force
			int lc   = fc.lc;	// loadcurve number
			double s = fc.s;	// force scale factor

			FENode& node = mesh.Node(id);

			double f = s*m_fem.GetLoadCurve(lc)->Value();

			int n = node.m_ID[bc];
			if ((bc == 0) && (n >= 0)) m_R[n] = f;
			if ((bc == 1) && (n >= 0)) m_R[n] = f;
			if ((bc == 2) && (n >= 0)) m_R[n] = f;
		}
	}

	// add contribution from domains
	for (int i=0; i<pstep->Domains(); ++i) 
	{
		FELinearElasticDomain& d = dynamic_cast<FELinearElasticDomain&>(*pstep->Domain(i));
		d.RHS(this, m_R);
	}

	// add contribution of linear surface loads
	for (int i=0; i<(int) m_fem.m_SL.size(); ++i)
	{
		FEPressureLoad* pl = dynamic_cast<FEPressureLoad*>(m_fem.m_SL[i]);
		if (pl && (pl->IsLinear())) pl->Residual(this, m_R);
	}
}

//-----------------------------------------------------------------------------
//! Reform the stiffness matrix. That is, calculate the shape of the stiffness
//! matrix (i.e. figure out the sparsity pattern), fill the stiffness matrix
//! with the element contributions and then factor the matrix. 
//! 
bool FELinearSolidSolver::ReformStiffness()
{
	// recalculate the shape of the stiffness matrix if necessary
	if (!CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	bool bret = StiffnessMatrix();

	if (bret)
	{
		m_SolverTime.start();
		{
			// factorize the stiffness matrix
			m_plinsolve->Factor();
		}
		m_SolverTime.stop();

		// increase total nr of reformations
		m_nref++;
		m_ntotref++;

		// reset bfgs update counter
		m_bfgs.m_nups = 0;
	}

	return bret;
}

//-----------------------------------------------------------------------------
//! Calculate the global stiffness matrix. This function simply calls 
//! FELinearSolidSolver::ElementStiffness() for each element and then assembles the
//! element stiffness matrix into the global matrix.
//!
bool FELinearSolidSolver::StiffnessMatrix()
{
	FEMesh& mesh = m_fem.m_mesh;

	// zero the stiffness matrix
	m_pK->Zero();

	// add contribution from domains
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	for (int i=0; i<pstep->Domains(); ++i)
	{
		FELinearElasticDomain& bd = dynamic_cast<FELinearElasticDomain&>(*pstep->Domain(i));
		bd.StiffnessMatrix(this);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Assembles the element stiffness matrix into the global stiffness matrix. 
//! This function also modifies the residual according to the prescribed
//! degrees of freedom. 
//!
void FELinearSolidSolver::AssembleStiffness(matrix& ke, vector<int>& lm)
{
	// assemble into the global stiffness
	m_pK->Assemble(ke, lm);

	// if there are prescribed bc's we need to adjust the residual
	if (m_fem.m_DC.size() > 0)
	{
		int i, j;
		int I, J;

		SparseMatrix& K = *m_pK;

		int N = ke.rows();

		// loop over columns
		for (j=0; j<N; ++j)
		{
			J = -lm[j]-2;
			if ((J >= 0) && (J<m_neq))
			{
				// dof j is a prescribed degree of freedom

				// loop over rows
				for (i=0; i<N; ++i)
				{
					I = lm[i];
					if (I >= 0)
					{
						// dof i is not a prescribed degree of freedom
						m_R[I] -= ke[i][j]*m_d[J];
					}
				}

				// set the diagonal element of K to 1
				K.set(J,J, 1);			
			}
		}
	}
}

//-----------------------------------------------------------------------------
//!  Assembles the element into the global residual. This function
//!  also checks for rigid dofs and assembles the residual using a condensing
//!  procedure in the case of rigid dofs.

void FELinearSolidSolver::AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R)
{
	// assemble the element residual into the global residual
	int ndof = fe.size();
	for (int i=0; i<ndof; ++i)
	{
		int I = elm[i];
		if ( I >= 0) R[I] += fe[i];
//		else if (-I-2 >= 0) m_Fr[-I-2] -= fe[i];	// reaction forces
	}
}

//-----------------------------------------------------------------------------
//! Store data to restart file
void FELinearSolidSolver::Serialize(DumpFile &ar)
{
	// TODO: implement serialization
}
