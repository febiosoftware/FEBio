#include "stdafx.h"
#include "FELinearSolver.h"
#include "FEModel.h"
#include "LinearSolver.h"
#include "FEGlobalMatrix.h"
#include "log.h"
#include "FENodeReorder.h"
#include "FELinearSystem.h"
#include "FEPrescribedBC.h"
#include "FEGlobalVector.h"
#include "FEDomain.h"
#include "FENodalLoad.h"
#include "FESurfaceLoad.h"
#include "FEBodyLoad.h"

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
//! add equations
void FELinearSolver::AddEquations(int neq)
{
	m_neq += neq;
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
		m_pls = fecore.CreateLinearSolver(fem);
		if (m_pls == 0)
		{
			felog.printbox("FATAL ERROR","Unknown solver type selected\n");
			return false;
		}

		if (m_part.empty() == false)
		{
			m_pls->SetPartitions(m_part);
		}
	}

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the solver allocate the correct type of matrix format
	SparseMatrix* pS = m_pls->CreateSparseMatrix(m_bsymm? REAL_SYMMETRIC : REAL_UNSYMMETRIC);
	if ((pS == 0) && m_bsymm)
	{
		// oh, oh, something went wrong. It's probably because the user requested a symmetric matrix for a 
		// solver that wants a non-symmetric. If so, let's force a non-symmetric format.
		pS = m_pls->CreateSparseMatrix(REAL_UNSYMMETRIC);

		if (pS)
		{
			// Problem solved! Let's inform the user.
			m_bsymm = false;
			felog.printbox("WARNING", "The matrix format was changed to non-symmetric since the selected\nlinear solver does not support a symmetric format. \n");
		}
	}

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
	m_pK = new FEGlobalMatrix(pS);
	if (m_pK == 0)
	{
		felog.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
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
	int nbc = fem.PrescribedBCs();
	for (int i=0; i<nbc; ++i)
	{
		FEPrescribedBC& dc = *fem.PrescribedBC(i);
		if (dc.IsActive()) dc.PrepStep(m_u, false);
	}

	// build the right-hand side
	// (Is done by the derived class)
	zero(m_R);
	vector<double> F(m_neq);
	FEGlobalVector rhs(fem, m_R, F);
	{
		TRACK_TIME("residual");
		ForceVector(rhs);
	}

	// increase RHS counter
	m_nrhs++;

	// build the stiffness matrix
	ReformStiffness();

	// solve the equations
	vector<double> u(m_neq);
	{
		TRACK_TIME("solve");
		if (m_pls->BackSolve(u, m_R) == false)
			throw LinearSolverFailed();
	}

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
		TRACK_TIME("reform");

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
		TRACK_TIME("stiffness");

		FELinearSystem K(*m_pK, m_R, m_u);
		if (!StiffnessMatrix(K)) return false;
	}

	// factorize the stiffness matrix
	{
		TRACK_TIME("solve");
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
	felog.printf("===== reforming stiffness matrix:\n");
	if (m_pK->Create(GetFEModel(), m_neq, true) == false) 
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
	{
		TRACK_TIME("solve");
		if (!m_pls->PreProcess()) throw FatalError();
	}

	// done!
	return true;
}

//-----------------------------------------------------------------------------
void FELinearSolver::Serialize(DumpStream& ar)
{
	FESolver::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_breform << m_neq;
	}
	else
	{
		ar >> m_breform >> m_neq;

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
	// Add nodal loads
	NodalLoads(R);
}

//-----------------------------------------------------------------------------
void FELinearSolver::NodalLoads(FEGlobalVector& R)
{
	// get the modal and mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// Get the DOF handler
	DOFS& dofs = fem.GetDOFS();

	// loop over nodal loads
	int ncnf = fem.NodalLoads();
	for (int i = 0; i<ncnf; ++i)
	{
		FENodalLoad& fc = *fem.NodalLoad(i);
		if (fc.IsActive())
		{
			// make sure the dof is valid
			int dof = fc.GetDOF();
			if ((dof >= 0) && (dof < dofs.GetTotalDOFS()))
			{
				int N = fc.Nodes();
				for (int j = 0; j < N; ++j)
				{
					int nid = fc.NodeID(j);
					FENode& node = mesh.Node(nid);
					int n = node.m_ID[dof];
					if (n >= 0)
					{
						R[n] = fc.NodeValue(j);
					}
				}
			}
			else assert(false);
		}
	}
}

//-----------------------------------------------------------------------------
// add surface loads to the RHS vector
void FELinearSolver::SurfaceLoads(FEGlobalVector& R)
{
	// get the time information
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	int nsl = fem.SurfaceLoads();
	for (int i = 0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = fem.SurfaceLoad(i);
		if (psl && psl->IsActive())
		{
			psl->Residual(tp, R);
		}
	}
}

//-----------------------------------------------------------------------------
// add body loads to RHS vector
void FELinearSolver::BodyLoads(FEGlobalVector& R)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	int nbl = fem.BodyLoads();
	for (int i = 0; i<nbl; ++i)
	{
		FEBodyLoad* pbl = fem.GetBodyLoad(i);
		if (pbl && pbl->IsActive())
		{
			pbl->ForceVector(R);
		}
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

	// update the domains
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.Update(tp);
	}
}
