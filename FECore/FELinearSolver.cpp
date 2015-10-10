#include "stdafx.h"
#include "FELinearSolver.h"
#include "FEModel.h"
#include "LinearSolver.h"
#include "FEGlobalMatrix.h"
#include "log.h"
#include "FENodeReorder.h"

//-----------------------------------------------------------------------------
//! constructor
FELinearSolver::FELinearSolver(FEModel* pfem) : FESolver(pfem)
{
	m_pls = 0;
	m_pK = 0;
	m_neq = 0;
}

//-----------------------------------------------------------------------------
void FELinearSolver::Clean()
{
	if (m_pls) m_pls->Destroy();
}

//-----------------------------------------------------------------------------
void FELinearSolver::SetDOF(vector<int>& dof)
{
	m_dof = dof;
}

//-----------------------------------------------------------------------------
bool FELinearSolver::Init()
{
	// Now that we have determined the equation numbers we can continue
	// with creating the stiffness matrix. First we select the linear solver
	// The stiffness matrix is created in CreateStiffness
	// Note that if a particular solver was requested in the input file
	// then the solver might already be allocated. That's way we need to check it.
	if (m_pls == 0)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		m_pls = fecore.CreateLinearSolver(m_fem.m_nsolver);
		if (m_pls == 0)
		{
			felog.printbox("FATAL ERROR","Unknown solver type selected\n");
			return false;
		}
	}

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the solver allocate the correct type of matrix format
	SparseMatrix* pS = m_pls->CreateSparseMatrix(m_bsymm? SPARSE_SYMMETRIC : SPARSE_UNSYMMETRIC);
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

	// get number of equations
	int neq = m_neq;

	// allocate data structures
	m_R.resize(neq);
	m_u.resize(neq);

	return true;
}

//-----------------------------------------------------------------------------
// Initialize linear equation system
bool FELinearSolver::InitEquations()
{
	FEMesh& mesh = m_fem.GetMesh();

	// initialize nr of equations
	int neq = 0;

	// degrees of freedom
	int ndof = m_dof.size();
	if (ndof == 0) return false;

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
			for (int j=0; j<ndof; ++j)
			{
				int dofj = m_dof[j];
				if      (node.m_ID[dofj] == DOF_FIXED     ) { node.m_ID[dofj] = -1; }
				else if (node.m_ID[dofj] == DOF_OPEN      ) { node.m_ID[dofj] =  neq++; }
				else if (node.m_ID[dofj] == DOF_PRESCRIBED) { node.m_ID[dofj] = -neq-2; neq++; }
				else { assert(false); return false; }
			}
		}
	}
	else
	{
		// give all free dofs an equation number
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			for (int j=0; j<ndof; ++j)
			{
				int dofj = m_dof[j];
				if      (node.m_ID[dofj] == DOF_FIXED     ) { node.m_ID[dofj] = -1; }
				else if (node.m_ID[dofj] == DOF_OPEN      ) { node.m_ID[dofj] =  neq++; }
				else if (node.m_ID[dofj] == DOF_PRESCRIBED) { node.m_ID[dofj] = -neq-2; neq++; }
				else { assert(false); return false; }
			}
		}
	}

	// store the number of equations
	m_neq = neq;

	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
//! Solve an analysis step
bool FELinearSolver::SolveStep(double time)
{
	// Make sure we have a linear solver and a stiffness matrix
	if (m_pls == 0) return false;
	if (m_pK == 0) return false;

	// reset counters
	m_niter = 0;
	m_nrhs = 0;
	m_nref = 0;
	m_ntotref = 0;

	// Set up the prescribed dof vector
	// The stiffness matrix assembler uses this to update the RHS vector
	// for prescribed dofs.
	zero(m_u);
	int nbc = m_fem.PrescribedBCs();
	for (int i=0; i<nbc; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive()) dc.PrepStep(m_u, false);
	}

	// build the right-hand side
	// (Is done by the derived class)
	zero(m_R);
	vector<double> F(m_neq);
	FEModel& fem = GetFEModel();
	FEGlobalVector rhs(fem, m_R, F);
	RHSVector(rhs);

	// build the stiffness matrix
	ReformStiffness();

	// solve the equations
	vector<double> u(m_neq);
	m_pls->BackSolve(u, m_R);

	// update solution
	Update(u);

	// increase iteration count
	m_niter++;

	return true;
}

//-----------------------------------------------------------------------------
//! assemble global stiffness matrix
void FELinearSolver::AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{
	// assemble into the global stiffness
	m_pK->Assemble(ke, elm);

	// if there are prescribed bc's we need to adjust the rhs vector
	if (m_fem.PrescribedBCs() > 0)
	{
		int i, j;
		int I, J;

		SparseMatrix& K = *m_pK;

		int N = ke.rows();

		// loop over columns
		for (j=0; j<N; ++j)
		{
			J = -elm[j]-2;
			if ((J >= 0) && (J<m_neq))
			{
				// dof j is a prescribed degree of freedom

				// loop over rows
				for (i=0; i<N; ++i)
				{
					I = elm[i];
					if (I >= 0)
					{
						// dof i is not a prescribed degree of freedom
						m_R[I] -= ke[i][j]*m_u[J];
					}
				}

				// set the diagonal element of K to 1
				K.set(J,J, 1);

				// set the rhs vector to the prescribed value
				// that way the solution vector will contain the prescribed value
				m_R[J] = m_u[J];
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FELinearSolver::ReformStiffness()
{
	// recalculate the shape of the stiffness matrix if necessary
	if (!CreateStiffness()) return false;

	// Make sure it is all set to zero
	m_pK->Zero();

	// calculate the stiffness matrix
	// (This is done by the derived class)
	if (!StiffnessMatrix()) return false;

	// factorize the stiffness matrix
	m_SolverTime.start();
	{
		m_pls->Factor();
	}
	m_SolverTime.stop();

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
	if (m_pK->Create(&GetFEModel(), m_neq, true) == false) 
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
		if (!m_pls->PreProcess()) throw FatalError();
	}
	m_SolverTime.stop();

	// done!
	return true;
}
