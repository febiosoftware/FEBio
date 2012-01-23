#include "stdafx.h"
#include "FEHeatSolver.h"
#include "FEBioLib/FEIsotropicFourier.h"
#include "FEBioLib/FEHeatFlux.h"
#include "FEHeatSolidDomain.h"
#include "FECore/FENodeReorder.h"

//-----------------------------------------------------------------------------
//! constructor for the class
FEHeatSolver::FEHeatSolver(FEModel &fem) : FESolver(fem)
{
}

//-----------------------------------------------------------------------------
//! Do one-time initialization for data
bool FEHeatSolver::Init()
{
	// initialize base class
	if (FESolver::Init() == false) return false;

	// get number of equations
	int neq = m_neq;

	// allocate data structures
	m_R.resize(neq);
	m_T.resize(neq);
	m_u.resize(neq);
	m_Tp.assign(neq, 0);

	// Identify the heat-transfer domains
	// TODO: I want this to be done automatically
	//       e.g. while the input file is being read
	FEMesh& mesh = m_fem.m_mesh;
	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FEHeatSolidDomain* pd = dynamic_cast<FEHeatSolidDomain*>(&mesh.Domain(nd));
		if (pd) m_Dom.push_back(nd);
	}
	assert(m_Dom.empty() == false);

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
bool FEHeatSolver::InitEquations()
{
	FEM& fem = dynamic_cast<FEM&>(m_fem);
	FEMesh& mesh = fem.m_mesh;

	// initialize nr of equations
	int neq = 0;

	// see if we need to optimize the bandwidth
	if (fem.m_bwopt)
	{
		// reorder the node numbers
		vector<int> P(mesh.Nodes());
		FENodeReorder mod;
		mod.Apply(mesh, P);

		// set the equation numbers
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(P[i]);
			if (node.m_ID[DOF_T] >= 0) node.m_ID[DOF_T] = neq++;
		}
	}
	else
	{
		// give all free dofs an equation number
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			if (node.m_ID[DOF_T] >= 0) node.m_ID[DOF_T] = neq++;
		}
	}

	// store the number of equations
	m_neq = neq;

	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
void FEHeatSolver::PrepStep()
{
	zero(m_u);
	for (size_t i=0; i<m_fem.m_DC.size(); ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;	// GAA

			double T = r + s*m_fem.GetLoadCurve(lc)->Value(); // GAA

			FENode& node = m_fem.m_mesh.Node(n);

			if (bc == DOF_T)
			{
				int I = -node.m_ID[bc]-2;
				if (I>=0 && I<m_neq) m_u[I] = T;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! solve a time step
bool FEHeatSolver::SolveStep(double time)
{
	// set up the prescribed temperatures vector
	PrepStep();

	// build the residual
	Residual();

	// build the stiffness matrix
	ReformStiffness();

	// solve the equations
	m_plinsolve->BackSolve(m_T, m_R);

	// update solution
	// NOTE: m_u is not being used in Update!
	Update(m_u);

	return true;
}

//-----------------------------------------------------------------------------
//! update solution
void FEHeatSolver::Update(vector<double>& u)
{
	FEMesh& mesh = m_fem.m_mesh;

	// update temperatures
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int n = node.m_ID[DOF_T];
		if (n >= 0) node.m_T = m_T[n];
		else if (-n-2 >= 0) node.m_T = m_T[-n-2] = m_u[-n-2];
	}

	// TODO: calculate fluxes

	// copy new temperatures to old temperature
	m_Tp = m_T;
}

//-----------------------------------------------------------------------------
//! Calculate the residual
void FEHeatSolver::Residual()
{
	// intialize residual to zero
	zero(m_R);

	// Add nodal flux contributions
	NodalFluxes(m_R);

	// add surface fluxes
	SurfaceFluxes(m_R);
}

//-----------------------------------------------------------------------------
//! Add nodal fluxes to residual
void FEHeatSolver::NodalFluxes(vector<double>& R)
{
	int i, id, bc, lc, n;
	double s, f;

	// get the FE mesh
	FEMesh& mesh = m_fem.m_mesh;

	// loop over nodal force cards
	int ncnf = m_fem.m_FC.size();
	for (i=0; i<ncnf; ++i)
	{
		FENodalForce& fc = *m_fem.m_FC[i];
		if (fc.IsActive())
		{
			id	 = fc.node;	// node ID
			bc   = fc.bc;	// direction of force
			lc   = fc.lc;	// loadcurve number
			s    = fc.s;	// force scale factor

			FENode& node = mesh.Node(id);

			n = node.m_ID[bc];
			if ((n >= 0) && (bc == DOF_T)) 
			{
				f = s*m_fem.GetLoadCurve(lc)->Value();
				R[n] = f;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate heat surface flux contribution to residual.
void FEHeatSolver::SurfaceFluxes(vector<double>& R)
{
	FEM& fem = dynamic_cast<FEM&>(m_fem);
	int nsl = (int) fem.m_SL.size();
	for (int i=0; i<nsl; ++i)
	{
		FEHeatFlux* phf = dynamic_cast<FEHeatFlux*>(fem.m_SL[i]);
		if (phf) phf->Residual(this, R);
	}
}

//-----------------------------------------------------------------------------
//! Reform the stiffness matrix. That is, calculate the shape of the stiffness
//! matrix (i.e. figure out the sparsity pattern), fill the stiffness matrix
//! with the element contributions and then factor the matrix. 
//! 
bool FEHeatSolver::ReformStiffness()
{
	// recalculate the shape of the stiffness matrix if necessary
	if (!CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	if (!StiffnessMatrix()) return false;

	// factorize the stiffness matrix
	m_SolverTime.start();
	{
		m_plinsolve->Factor();
	}
	m_SolverTime.stop();

	// increase total nr of reformations
	m_nref++;
	m_ntotref++;

	// reset bfgs update counter
	m_bfgs.m_nups = 0;

	return true;
}

//-----------------------------------------------------------------------------
//! Calculate the global stiffness matrix. This function simply calls 
//! HeatStiffnessMatrix() for each domain which will calculate the
//! contribution to the global stiffness matrix from each domain.
bool FEHeatSolver::StiffnessMatrix()
{
	// see if this is a dynamic problem
	bool bdyn = (m_fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC);

	// get the time step size
	double dt = m_fem.GetCurrentStep()->m_dt;

	// zero the stiffness matrix
	m_pK->Zero();

	// Add stiffness contribution from all domains
	for (int i=0; i<(int) m_Dom.size(); ++i)
	{
		FEHeatSolidDomain& bd = dynamic_cast<FEHeatSolidDomain&>(*Domain(i));

		// add the conduction stiffness
		bd.ConductionMatrix(this);

		// for a dynamic analysis add the capacitance matrix
		if (bdyn) bd.CapacitanceMatrix(this, dt);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Assembles the element stiffness matrix into the global stiffness matrix. 
//! This function also modifies the residual according to the prescribed
//! degrees of freedom. 
//!
void FEHeatSolver::AssembleStiffness(vector<int>& en, vector<int>& lm, matrix& ke)
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
						m_R[I] -= ke[i][j]*m_u[J];
					}
				}

				// set the diagonal element of K to 1
				K.set(J,J, 1);			
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Serializes data to the archive.
//! Still need to implement this.
//!
void FEHeatSolver::Serialize(DumpFile &ar)
{

}
