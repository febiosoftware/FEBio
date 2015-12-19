#include "stdafx.h"
#include "FELinearSolidSolver.h"
#include "FELinearSolidDomain.h"
#include "FECore/FENodeReorder.h"
#include "FEPressureLoad.h"
#include "FERigidMaterial.h"
#include "FECore/log.h"
#include "FEStiffnessMatrix.h"
#include "NumCore/NumCore.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FELinearSolidSolver, FESolver)
	ADD_PARAMETER(m_Dtol         , FE_PARAM_DOUBLE, "dtol"    );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Class constructor
FELinearSolidSolver::FELinearSolidSolver(FEModel* pfem) : FENewtonSolver(pfem)
{
	m_Dtol = 1e-9;

	m_pK = 0;
	m_neq = 0;
	m_plinsolve = 0;

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int varD = dofs.AddVariable("displacement");
	dofs.AddDOF(varD, "x");
	dofs.AddDOF(varD, "y");
	dofs.AddDOF(varD, "z");

	// get DOF indices
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
}

//-----------------------------------------------------------------------------
FELinearSolidSolver::~FELinearSolidSolver()
{
}

//-----------------------------------------------------------------------------
//! This function will solve the FE problem
//!
bool FELinearSolidSolver::Quasin(double time)
{
	// prepare step
	FEMaterialPoint::dt = m_fem.GetCurrentStep()->m_dt;
	FEMaterialPoint::time = m_fem.m_ftime;

	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).InitElements();

	// set-up the prescribed displacements
	zero(m_d);
	vector<double> DT(m_d), DI(m_d);
	int nbc = m_fem.PrescribedBCs();
	for (int i=0; i<nbc; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive())
		{
			int bc = dc.GetDOF();
			for (size_t j = 0; j<dc.Items(); ++j)
			{
				double D = dc.NodeValue(j);
				int n = dc.NodeID(j);

				FENode& node = mesh.Node(n);

				if (bc == m_dofX) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_neq) { DT[I] = D; DI[I] = D - node.get(m_dofX); }}
				if (bc == m_dofY) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_neq) { DT[I] = D; DI[I] = D - node.get(m_dofY); }}
				if (bc == m_dofZ) { int I = -node.m_ID[bc]-2; if (I>=0 && I<m_neq) { DT[I] = D; DI[I] = D - node.get(m_dofZ); }}
			}
		}
	}

	FETimePoint tp = m_fem.GetTime();

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
		ReformStiffness(tp);

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

		felog.printf("normu = %lg\n", normu);

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
	FEMesh& mesh = m_fem.GetMesh();

	// update nodal positions
	int n;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		n = node.m_ID[m_dofX]; if (n >= 0) node.m_rt.x = node.m_r0.x + u[n]; else if (-n-2 >= 0) node.m_rt.x = node.m_r0.x + m_d[-n-2];
		n = node.m_ID[m_dofY]; if (n >= 0) node.m_rt.y = node.m_r0.y + u[n]; else if (-n-2 >= 0) node.m_rt.y = node.m_r0.y + m_d[-n-2];
		n = node.m_ID[m_dofZ]; if (n >= 0) node.m_rt.z = node.m_r0.z + u[n]; else if (-n-2 >= 0) node.m_rt.z = node.m_r0.z + m_d[-n-2];
	}

	// update the stresses on all domains
	for (int i=0; i<pstep->Domains(); ++i)
	{
		FELinearElasticDomain& d = dynamic_cast<FELinearElasticDomain&>(*pstep->Domain(i));
		d.UpdateStresses(m_fem);
	}

	// output modified state data
	m_fem.Write(FE_UNCONVERGED);
}

//-----------------------------------------------------------------------------
//! Calculate the residual
void FELinearSolidSolver::Residual()
{
	zero(m_R);

	vector<double> dummy(m_R);

	FEGlobalVector RHS(GetFEModel(), m_R, dummy);

	FEMesh& mesh = m_fem.GetMesh();
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// loop over nodal forces
	int ncnf = m_fem.NodalLoads();
	for (int i=0; i<ncnf; ++i)
	{
		FENodalLoad& fc = *m_fem.NodalLoad(i);
		if (fc.IsActive())
		{
			int id	 = fc.m_node;	// node ID
			int bc   = fc.m_bc;	// direction of force

			FENode& node = mesh.Node(id);

			double f = fc.Value();

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
		d.RHS(RHS);
	}

	// add contribution of linear surface loads
	int nsl = m_fem.SurfaceLoads();
	for (int i=0; i<nsl; ++i)
	{
		FEPressureLoad* pl = dynamic_cast<FEPressureLoad*>(m_fem.SurfaceLoad(i));
		if (pl && (pl->IsLinear())) pl->Residual(RHS);
	}
}

//-----------------------------------------------------------------------------
//! Calculate the global stiffness matrix. This function simply calls 
//! FELinearSolidSolver::ElementStiffness() for each element and then assembles the
//! element stiffness matrix into the global matrix.
//!
bool FELinearSolidSolver::StiffnessMatrix(const FETimePoint& tp)
{
	FEMesh& mesh = m_fem.GetMesh();

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
void FELinearSolidSolver::AssembleStiffness(vector<int>& en, vector<int>& lm, matrix& ke)
{
	// assemble into the global stiffness
	m_pK->Assemble(ke, lm);

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
//! \todo Implement serialization
void FELinearSolidSolver::Serialize(DumpFile &ar)
{
	
}
