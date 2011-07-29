#include "stdafx.h"
#include "FELinearSolidSolver.h"
#include "FELinearSolidDomain.h"

//-----------------------------------------------------------------------------
//! Class constructor
FELinearSolidSolver::FELinearSolidSolver(FEM& fem) : FESolver(fem)
{
	m_niter = 0;
}

//-----------------------------------------------------------------------------
//! Initialization
//! This function will be called before each analysis step
bool FELinearSolidSolver::Init()
{
	// initialize base class
	if (FESolver::Init() == false) return false;

	// get the nr of equations
	int neq = m_fem.m_neq;

	// allocate data structures
	m_u.resize(neq);
	m_R.resize(neq);
	m_d.resize(neq);

	return true;
}

//-----------------------------------------------------------------------------
//! This function will solve the FE problem
//!
bool FELinearSolidSolver::SolveStep(double time)
{
	// set-up the prescribed displacements
	zero(m_d);
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

			if (bc == 0) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_fem.m_neq) m_d[I] = D; }
			if (bc == 1) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_fem.m_neq) m_d[I] = D; }
			if (bc == 2) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_fem.m_neq) m_d[I] = D; }
		}
	}

	// build the residual
	Residual();

	// build the stiffness matrix
	ReformStiffness();

	// solve the equations
	m_plinsolve->BackSolve(m_u, m_R);

	// update solution
	// NOTE: m_u is not being used in Update!
	Update(m_u);

	return true;
}

//-----------------------------------------------------------------------------
//! update solution
void FELinearSolidSolver::Update(vector<double>& u)
{
	FEMesh& mesh = m_fem.m_mesh;

	// update nodal positions
	int n;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		n = node.m_ID[0]; if (n >= 0) node.m_rt.x = node.m_r0.x + u[n]; else if (-n-2 >= 0) node.m_rt.x = node.m_r0.x + m_d[-n-2];
		n = node.m_ID[1]; if (n >= 0) node.m_rt.y = node.m_r0.y + u[n]; else if (-n-2 >= 0) node.m_rt.y = node.m_r0.y + m_d[-n-2];
		n = node.m_ID[2]; if (n >= 0) node.m_rt.z = node.m_r0.z + u[n]; else if (-n-2 >= 0) node.m_rt.z = node.m_r0.z + m_d[-n-2];
	}

	// update the stresses on all domains
	for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).UpdateStresses(m_fem);
}

//-----------------------------------------------------------------------------
//! Calculate the residual
void FELinearSolidSolver::Residual()
{
	zero(m_R);

	// apply nodal fluxes
	int i, id, bc, lc, n;
	double s, f;

	FEMesh& mesh = m_fem.m_mesh;

	// loop over nodal forces
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

			f = s*m_fem.GetLoadCurve(lc)->Value();

			n = node.m_ID[bc];
			if ((bc == 0) && (n >= 0)) m_R[n] = f;
			if ((bc == 1) && (n >= 0)) m_R[n] = f;
			if ((bc == 2) && (n >= 0)) m_R[n] = f;
		}
	}

	// TODO: surface tractions

	// TODO: initial stress
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

	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FELinearSolidDomain& bd = dynamic_cast<FELinearSolidDomain&>(mesh.Domain(nd));
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

		int neq = m_fem.m_neq;

		// loop over columns
		for (j=0; j<N; ++j)
		{
			J = -lm[j]-2;
			if ((J >= 0) && (J<m_fem.m_nreq))
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
//! Store data to restart file
void FELinearSolidSolver::Serialize(DumpFile &ar)
{
	// TODO: implement serialization
}
