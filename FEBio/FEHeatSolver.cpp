#include "stdafx.h"
#include "FEHeatSolver.h"
#include "FEIsotropicFourier.h"

//-----------------------------------------------------------------------------
//! constructor for the class
FEHeatSolver::FEHeatSolver(FEM &fem) : FESolver(fem)
{
	m_niter = 0;
}

//-----------------------------------------------------------------------------
//! Do one-time initialization for data
bool FEHeatSolver::Init()
{
	// initialize base class
	if (FESolver::Init() == false) return false;

	// get the nr of equations
	int neq = m_fem.m_neq;

	// allocate data structures
	m_R.resize(neq);
	m_T.resize(neq);
	m_u.resize(neq);
	m_Tp.assign(neq, 0);

	return true;
}

//-----------------------------------------------------------------------------
//! solve a time step
bool FEHeatSolver::SolveStep(double time)
{
	// set up the prescribed temperatures vector
	zero(m_u);
	for (size_t i=0; i<m_fem.m_DC.size(); ++i)
	{
		FENodalDisplacement& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;

			double T = s*m_fem.GetLoadCurve(lc)->Value();

			FENode& node = m_fem.m_mesh.Node(n);

			if (bc == 10)
			{
				int I = -node.m_ID[bc]-2;
				if (I>=0 && I<m_fem.m_neq) m_u[I] = T;
			}
		}
	}

	// build the residual
	Residual();

	// build the stiffness matrix
	ReformStiffness();

	// solve the equations
	m_plinsolve->Solve(m_T, m_R);

	// update solution
	Update();

	return true;
}

//-----------------------------------------------------------------------------
//! update solution
void FEHeatSolver::Update()
{
	FEMesh& mesh = m_fem.m_mesh;

	// update temperatures
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int n = node.m_ID[10];
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
	zero(m_R);

	// apply nodal fluxes
	int i, j, id, bc, lc, n;
	double s, f;

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
			s    = fc.s;		// force scale factor

			FENode& node = mesh.Node(id);

			n = node.m_ID[bc];
			if ((n >= 0) && (bc == 10)) 
			{
				f = s*m_fem.GetLoadCurve(lc)->Value();
				m_R[n] = f;
			}
		}
	}

	// add surface fluxes
	if (m_fem.m_phflux)
	{
		int nfc = m_fem.m_phflux->Elements();
		for (i=0; i<nfc; ++i)
		{
			FEHeatFlux& hf = m_fem.m_phflux->HeatFlux(i);
			FESurfaceElement& el = m_fem.m_phflux->Element(i);
			m_fem.m_phflux->UnpackElement(el);

			int ne = el.Nodes();
			int ni = el.GaussPoints();

			double g = m_fem.GetLoadCurve(hf.lc)->Value();

			// calculate nodal fluxes
			double qn[4];
			for (j=0; j<el.Nodes(); ++j) qn[j] = g*hf.s[j];

			vector<double> fe(ne);

			// nodal coordinates
			vec3d *rt = el.rt();

			double* Gr, *Gs;
			double* N;
			double* w  = el.GaussWeights();

			// pressure at integration points
			double q;

			vec3d dxr, dxs;

			vector<int> lm(ne);
			for (j=0; j<ne; ++j) lm[j] = (el.LM())[ne*10 + j];

			// force vector
			// repeat over integration points
			zero(fe);
			for (n=0; n<ni; ++n)
			{
				N  = el.H(n);
				Gr = el.Gr(n);
				Gs = el.Gs(n);

				q = 0;
				dxr = dxs = vec3d(0,0,0);
				for (j=0; j<ne; ++j) 
				{
					q += N[j]*qn[j];
					dxr.x += Gr[j]*rt[j].x;
					dxr.y += Gr[j]*rt[j].y;
					dxr.z += Gr[j]*rt[j].z;

					dxs.x += Gs[j]*rt[j].x;
					dxs.y += Gs[j]*rt[j].y;
					dxs.z += Gs[j]*rt[j].z;
				}
		
				double J = (dxr ^ dxs).norm();

				for (j=0; j<ne; ++j) fe[j] += N[j]*q*J*w[n];
			}

			// add element force vector to global force vector
			for (j=0; j<ne; ++j)
			{
				if (lm[j] >= 0) m_R[lm[j]] += fe[j];
			}
		}
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
//! FEHeatSolver::ElementStiffness() for each element and then assembles the
//! element stiffness matrix into the global matrix.
//!
bool FEHeatSolver::StiffnessMatrix()
{
	FEMesh& mesh = m_fem.m_mesh;

	// zero the stiffness matrix
	m_pK->Zero();

	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FEHeatSolidDomain& bd = dynamic_cast<FEHeatSolidDomain&>(mesh.Domain(nd));
		bd.HeatStiffnessMatrix(this);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Assembles the element stiffness matrix into the global stiffness matrix. 
//! This function also modifies the residual according to the prescribed
//! degrees of freedom. 
//!
void FEHeatSolver::AssembleStiffness(matrix& ke, vector<int>& lm)
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
