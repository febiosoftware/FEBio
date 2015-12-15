#include "FEHeatSolver.h"
#include "FEHeatFlux.h"
#include "FEConvectiveHeatFlux.h"
#include "FEHeatTransferMaterial.h"
#include "FEHeatSource.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEHeatSolver, FELinearSolver)
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor for the class
FEHeatSolver::FEHeatSolver(FEModel* pfem) : FELinearSolver(pfem)
{
	m_brhs = false;
	m_ntotref = 0;
	m_niter = 0;
	m_nrhs = 0;

	// set the active degrees of freedom for this solver
	const int dof_T = pfem->GetDOFS().GetDOF("t");
	vector<int> dof;
	dof.push_back(dof_T);
	SetDOF(dof);
}

//-----------------------------------------------------------------------------
FEHeatSolver::~FEHeatSolver()
{
}

//-----------------------------------------------------------------------------
//! Do one-time initialization for data
bool FEHeatSolver::Init()
{
	// Call base class first
	if (FELinearSolver::Init() == false) return false;

	const int dof_T = m_fem.GetDOFS().GetDOF("t");
	if (dof_T == -1) { assert(false); return false; }

	// allocate data structures
	int neq = NumberOfEquations();
	m_Tp.assign(neq, 0);

	// set initial temperatures
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int nid = node.m_ID[dof_T];
		if (nid >= 0) m_Tp[nid] = node.get(dof_T);
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
//! update solution
void FEHeatSolver::Update(vector<double>& u)
{
	FEMesh& mesh = m_fem.GetMesh();
	const int dof_T = m_fem.GetDOFS().GetDOF("t");
	if (dof_T == -1) { assert(false); return; }

	// update temperatures
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int n = node.m_ID[dof_T];
		if (n >= 0) node.set(dof_T, u[n]);
		else if (-n-2 >= 0) node.set(dof_T, u[-n-2]);
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
				for (int n=0; n<ne; ++n) T[n] = mesh.Node(el.m_node[n]).get(dof_T);

				// calculate heat flux for each integration point
				for (int n=0; n<ni; ++n)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(n);
					FEHeatMaterialPoint* pt = (mp.ExtractData<FEHeatMaterialPoint>());
					assert(pt);

					vec3d gradT = pbd->gradient(el, T, n);
					mat3ds D = pmat->Conductivity(mp);
					pt->m_q = -(D*gradT);
				}
			}
		}
	}

	// copy new temperatures to old temperature
	m_Tp = u;
}

//-----------------------------------------------------------------------------
//! Calculate the residual
void FEHeatSolver::RHSVector(FEGlobalVector& R)
{
	// Add nodal flux contributions
	NodalFluxes(R);

	// add surface fluxes
	SurfaceFluxes(R);

	// heat sources
	HeatSources(R);
}

//-----------------------------------------------------------------------------
//! Add nodal fluxes to residual
void FEHeatSolver::NodalFluxes(FEGlobalVector& R)
{
	// get the FE mesh
	FEMesh& mesh = m_fem.GetMesh();
	const int dof_T = m_fem.GetDOFS().GetDOF("t");
	if (dof_T == -1) { assert(false); return; }

	// loop over nodal loads
	int ncnf = m_fem.NodalLoads();
	for (int i=0; i<ncnf; ++i)
	{
		FENodalLoad& fc = *m_fem.NodalLoad(i);
		if (fc.IsActive())
		{
			int id	 = fc.m_node;	// node ID
			int bc   = fc.m_bc;		// direction of force

			FENode& node = mesh.Node(id);

			int n = node.m_ID[bc];
			if ((n >= 0) && (bc == dof_T)) 
			{
				R[n] = fc.Value();
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
//! This function is modified from the base class for including capacitance matrix
void FEHeatSolver::AssembleStiffness(vector<int>& en, vector<int>& lm, matrix& ke)
{
	// Call base class
	FELinearSolver::AssembleStiffness(en, lm, ke);

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
}

//-----------------------------------------------------------------------------
//! Serializes data to the archive.
//! Still need to implement this.
//!
void FEHeatSolver::Serialize(DumpFile &ar)
{
	FELinearSolver::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_Tp;
		ar << m_brhs;
	}
	else
	{
		ar >> m_Tp;
		ar >> m_brhs;
	}
}
