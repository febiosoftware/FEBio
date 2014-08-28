#include "FEHeatSolver.h"
#include "FEHeatFlux.h"
#include "FEConvectiveHeatFlux.h"
#include "FEHeatTransferMaterial.h"
#include "FEHeatSource.h"
#include "FEHeatStiffnessMatrix.h"
#include "FECore/FENodeReorder.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "NumCore/NumCore.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEHeatSolver, FESolver)
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor for the class
FEHeatSolver::FEHeatSolver(FEModel* pfem) : FESolver(pfem)
{
	m_brhs = false;
	m_ntotref = 0;
	m_niter = 0;
	m_nrhs = 0;

	m_pK = 0;
	m_neq = 0;
	m_plinsolve = 0;
}

//-----------------------------------------------------------------------------
FEHeatSolver::~FEHeatSolver()
{
	delete m_plinsolve;	// clean up linear solver data
	delete m_pK;		// clean up stiffnes matrix data
}

//-----------------------------------------------------------------------------
//! Clean
//! \todo Why can this not be done in destructor?
void FEHeatSolver::Clean()
{
	if (m_plinsolve) m_plinsolve->Destroy();
}

//-----------------------------------------------------------------------------
//! Do one-time initialization for data
bool FEHeatSolver::Init()
{
	// Now that we have determined the equation numbers we can continue
	// with creating the stiffness matrix. First we select the linear solver
	// The stiffness matrix is created in CreateStiffness
	// Note that if a particular solver was requested in the input file
	// then the solver might already be allocated. That's way we need to check it.
	if (m_plinsolve == 0)
	{
		m_plinsolve = NumCore::CreateLinearSolver(m_fem.m_nsolver);
		if (m_plinsolve == 0)
		{
			felog.printbox("FATAL ERROR","Unknown solver type selected\n");
			return false;
		}
	}

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the solver allocate the correct type of matrix format
	SparseMatrix* pS = m_plinsolve->CreateSparseMatrix(m_bsymm? SPARSE_SYMMETRIC : SPARSE_UNSYMMETRIC);
	if (pS == 0)
	{
		felog.printbox("FATAL ERROR", "The selected linear solver does not support the requested\n matrix format.\nPlease select a different linear solver.\n");
		return false;
	}

	// clean up the stiffness matrix if we have one
	if (m_pK) delete m_pK; m_pK = 0;

	// Create the stiffness matrix.
	// Note that this does not construct the stiffness matrix. This
	// is done later in the StiffnessMatrix routine.
	m_pK = new FEHeatStiffnessMatrix(pS);
	if (m_pK == 0)
	{
		felog.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
		return false;
	}

	// get number of equations
	int neq = m_neq;

	// allocate data structures
	m_R.resize(neq);
	m_T.resize(neq);
	m_u.resize(neq);
	m_Tp.assign(neq, 0);

	// set initial temperatures
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int nid = node.m_ID[DOF_T];
		if (nid >= 0) m_Tp[nid] = node.m_T;
	}

	// Identify the heat-transfer domains
	// TODO: I want this to be done automatically
	//       e.g. while the input file is being read
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	pstep->ClearDomains();
	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FEHeatSolidDomain* pd = dynamic_cast<FEHeatSolidDomain*>(&mesh.Domain(nd));
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
bool FEHeatSolver::InitEquations()
{
	FEMesh& mesh = m_fem.GetMesh();

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
	int nbc = m_fem.PrescribedBCs();
	for (int i=0; i<nbc; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;

			double T = r + s*m_fem.GetLoadCurve(lc)->Value();

			FENode& node = m_fem.GetMesh().Node(n);

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
	m_niter = 0;
	m_nrhs = 0;
	m_nref = 0;
	m_ntotref = 0;

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

	// increase iteration count
	m_niter++;

	return true;
}

//-----------------------------------------------------------------------------
//! update solution
void FEHeatSolver::Update(vector<double>& u)
{
	FEMesh& mesh = m_fem.GetMesh();

	// update temperatures
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int n = node.m_ID[DOF_T];
		if (n >= 0) node.m_T = m_T[n];
		else if (-n-2 >= 0) node.m_T = m_T[-n-2] = m_u[-n-2];
	}

	// update heat fluxes
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEHeatSolidDomain* pbd = dynamic_cast<FEHeatSolidDomain*>(&mesh.Domain(i));
		if (pbd)
		{
			FEHeatTransferMaterial* pmat = dynamic_cast<FEHeatTransferMaterial*>(pbd->GetMaterial());
			assert(pmat);

			int NE = pbd->Elements();
			for (int j=0; j<NE; ++j)
			{
				FESolidElement& el = pbd->Element(j);
				int ni = el.GaussPoints();
				int ne = el.Nodes();

				// get the nodal temperatures
				double T[FEElement::MAX_NODES];
				for (int n=0; n<ne; ++n) T[n] = mesh.Node(el.m_node[n]).m_T;

				// calculate heat flux for each integration point
				for (int n=0; n<ni; ++n)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(n);
					FEHeatMaterialPoint* pt = (mp.ExtractData<FEHeatMaterialPoint>());
					assert(pt);

					vec3d gradT = pbd->gradient(el, T, n);
					pt->m_q = pmat->HeatFlux(gradT);
				}
			}
		}
	}

	// copy new temperatures to old temperature
	m_Tp = m_T;
}

//-----------------------------------------------------------------------------
//! Calculate the residual
void FEHeatSolver::Residual()
{
	vector<double> dummy(m_R);

	FEGlobalVector RHS(GetFEModel(), m_R, dummy);

	// intialize residual to zero
	zero(m_R);

	// Add nodal flux contributions
	NodalFluxes(RHS);

	// add surface fluxes
	SurfaceFluxes(RHS);

	// heat sources
	HeatSources(RHS);

	// increase RHS counter
	m_nrhs++;
}

//-----------------------------------------------------------------------------
//! Add nodal fluxes to residual
void FEHeatSolver::NodalFluxes(FEGlobalVector& R)
{
	int i, id, bc, lc, n;
	double s, f;

	// get the FE mesh
	FEMesh& mesh = m_fem.GetMesh();

	// loop over nodal force cards
	int ncnf = m_fem.NodalLoads();
	for (i=0; i<ncnf; ++i)
	{
		FENodalForce& fc = *m_fem.NodalLoad(i);
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
void FEHeatSolver::SurfaceFluxes(FEGlobalVector& R)
{
	int nsl = m_fem.SurfaceLoads();
	for (int i=0; i<nsl; ++i)
	{
		// heat flux
		FEHeatFlux* phf = dynamic_cast<FEHeatFlux*>(m_fem.SurfaceLoad(i));
		if (phf && phf->IsActive()) phf->Residual(R);

		// convective heat flux
		FEConvectiveHeatFlux* pchf = dynamic_cast<FEConvectiveHeatFlux*>(m_fem.SurfaceLoad(i));
		if (pchf && pchf->IsActive()) pchf->Residual(R);
	}
}

//-----------------------------------------------------------------------------
//! Calculate the heat generation from heat sources
void FEHeatSolver::HeatSources(FEGlobalVector& R)
{
	int nbl = m_fem.BodyLoads();
	for (int i=0; i<nbl; ++i)
	{
		FEHeatSource* psh = dynamic_cast<FEHeatSource*>(m_fem.GetBodyLoad(i));
		if (psh) psh->Residual(R);
	}
}

//-----------------------------------------------------------------------------
//!  Creates the global stiffness matrix
//! \todo Can we move this to the FEStiffnessMatrix::Create function?
bool FEHeatSolver::CreateStiffness(bool breset)
{
	// clean up the solver
	if (m_pK->NonZeroes()) m_plinsolve->Destroy();

	// clean up the stiffness matrix
	m_pK->Clear();

	// create the stiffness matrix
	felog.printf("===== reforming stiffness matrix:\n");
	if (m_pK->Create(&GetFEModel(), m_neq, breset) == false) 
	{
		felog.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
		return false;
	}
	else
	{
		// output some information about the direct linear solver
		int neq = m_pK->Rows();
		int nnz = m_pK->NonZeroes();
		felog.printf("\tNr of equations ........................... : %d\n", neq);
		felog.printf("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
		felog.printf("\n");
	}

	// Do the preprocessing of the solver
	m_SolverTime.start();
	{
		if (!m_plinsolve->PreProcess()) throw FatalError();
	}
	m_SolverTime.stop();

	// done!
	return true;
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

	return true;
}

//-----------------------------------------------------------------------------
//! Calculate the global stiffness matrix. This function simply calls 
//! HeatStiffnessMatrix() for each domain which will calculate the
//! contribution to the global stiffness matrix from each domain.
bool FEHeatSolver::StiffnessMatrix()
{
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// see if this is a dynamic problem
	bool bdyn = (pstep->m_nanalysis == FE_DYNAMIC);

	// get the time step size
	double dt = m_fem.GetCurrentStep()->m_dt;

	// zero the stiffness matrix
	m_pK->Zero();

	// Add stiffness contribution from all domains
	for (int i=0; i<pstep->Domains(); ++i)
	{
		FEHeatDomain& bd = dynamic_cast<FEHeatDomain&>(*pstep->Domain(i));

		// add the conduction stiffness
		m_brhs = false;
		bd.ConductionMatrix(this);

		// for a dynamic analysis add the capacitance matrix
		if (bdyn) 
		{
			m_brhs = true;
			bd.CapacitanceMatrix(this, dt);
		}
	}

	// add convective heat fluxes
	m_brhs = false;
	for (int i=0; i<m_fem.SurfaceLoads(); ++i)
	{
		FEConvectiveHeatFlux* pbc = dynamic_cast<FEConvectiveHeatFlux*>(m_fem.SurfaceLoad(i));
		if (pbc && pbc->IsActive()) pbc->StiffnessMatrix(this);
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

	// see if we need to modify the RHS
	// (This is needed for the capacitance matrix)
	if (m_brhs)
	{
		int ne = (int) lm.size();
		for (int j=0; j<ne; ++j)
		{
			if (lm[j] >= 0)
			{
				double q = 0;
				for (int k=0; k<ne; ++k)
				{
					if (lm[k] >= 0) q += ke[j][k]*m_Tp[lm[k]];
					else if (-lm[k]-2 >= 0) q += ke[j][k]*m_Tp[-lm[k]-2];
				}
				m_R[lm[j]] += q;
			}
		}
	}

	// if there are prescribed bc's we need to adjust the residual
	if (m_fem.PrescribedBCs() > 0)
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
	FESolver::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_T << m_Tp << m_R << m_u;
		ar << m_brhs;
	}
	else
	{
		ar >> m_T >> m_Tp >> m_R >> m_u;
		ar >> m_brhs;
	}

}
