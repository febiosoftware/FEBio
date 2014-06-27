#include "stdafx.h"
#include "FESolidSolver2.h"
#include "FERigidMaterial.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FEPointBodyForce.h"
#include "FEPressureLoad.h"
#include "FEResidualVector.h"
#include "FECore/FENodeReorder.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "NumCore/NumCore.h"

#ifdef WIN32
	#include <float.h>
	#define ISNAN(x) _isnan(x)
#endif

#ifdef LINUX
	#include <math.h>
	#define ISNAN(x) isnan(x)
#endif

#ifdef __APPLE__
#include <math.h>
#define ISNAN(x) isnan(x)
#endif

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FESolidSolver2, FESolver)
	ADD_PARAMETER(m_Dtol         , FE_PARAM_DOUBLE, "dtol"        );
	ADD_PARAMETER(m_Etol         , FE_PARAM_DOUBLE, "etol"        );
	ADD_PARAMETER(m_Rtol         , FE_PARAM_DOUBLE, "rtol"        );
	ADD_PARAMETER(m_Rmin         , FE_PARAM_DOUBLE, "min_residual");
	ADD_PARAMETER(m_bfgs.m_LStol , FE_PARAM_DOUBLE, "lstol"       );
	ADD_PARAMETER(m_bfgs.m_LSmin , FE_PARAM_DOUBLE, "lsmin"       );
	ADD_PARAMETER(m_bfgs.m_LSiter, FE_PARAM_INT   , "lsiter"      );
	ADD_PARAMETER(m_bfgs.m_maxref, FE_PARAM_INT   , "max_refs"    );
	ADD_PARAMETER(m_bfgs.m_maxups, FE_PARAM_INT   , "max_ups"     );
	ADD_PARAMETER(m_bfgs.m_cmax  , FE_PARAM_DOUBLE, "cmax"        );
    ADD_PARAMETER(m_alpha        , FE_PARAM_DOUBLE, "alpha"       );
	ADD_PARAMETER(m_beta         , FE_PARAM_DOUBLE, "beta"        );
	ADD_PARAMETER(m_gamma        , FE_PARAM_DOUBLE, "gamma"       );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FESolidSolver Construction
//
FESolidSolver2::FESolidSolver2(FEModel* pfem) : FESolver(pfem)
{
	// default values
	m_Rtol = 0;	// deactivate residual convergence 
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_Rmin = 1.0e-20;

	m_niter = 0;
	m_nreq = 0;

	m_pK = 0;
	m_neq = 0;
	m_plinsolve = 0;


	// default Newmark parameters for unconditionally stable time integration
	m_beta = 0.25;
	m_gamma = 0.5;
	m_alpha = 1.;
}

//-----------------------------------------------------------------------------
FESolidSolver2::~FESolidSolver2()
{
	delete m_plinsolve;	// clean up linear solver data
	delete m_pK;		// clean up stiffnes matrix data
}

//-----------------------------------------------------------------------------
//! Clean
//! \todo Why can this not be done in destructor?
void FESolidSolver2::Clean()
{
	if (m_plinsolve) m_plinsolve->Destroy();
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FESolidSolver
//
bool FESolidSolver2::Init()
{
	// check parameters
	if (m_Dtol <  0.0) { felog.printf("Error: dtol must be nonnegative.\n"   ); return false; }
	if (m_Etol <  0.0) { felog.printf("Error: etol must be nonnegative.\n"); return false; }
	if (m_Rtol <  0.0) { felog.printf("Error: rtol must be nonnegative.\n"); return false; }
	if (m_Rmin <  0.0) { felog.printf("Error: min_residual must be nonnegative.\n"  ); return false; }
	if (m_bfgs.m_LStol  < 0.0) { felog.printf("Error: lstol must be nonnegative.\n" ); return false; }
	if (m_bfgs.m_LSmin  < 0.0) { felog.printf("Error: lsmin must be nonnegative.\n" ); return false; }
	if (m_bfgs.m_LSiter < 0) { felog.printf("Error: lsiter must be nonnegative.\n"  ); return false; }
	if (m_bfgs.m_maxref < 0) { felog.printf("Error: max_refs must be nonnegative.\n"); return false; }
	if (m_bfgs.m_maxups < 0) { felog.printf("Error: max_ups must be nonnegative.\n" ); return false; }
	if (m_bfgs.m_cmax   < 0) { felog.printf("Error: cmax must be nonnegative.\n"    ); return false; }

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
	m_pK = new FEStiffnessMatrix(pS);
	if (m_pK == 0)
	{
		felog.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
		return false;
	}

	// get nr of equations
	int neq = m_neq;

	// allocate vectors
	m_Fn.assign(neq, 0);
	m_Fd.assign(neq, 0);
	m_Fr.assign(neq, 0);
	m_Ui.assign(neq, 0);
	m_Ut.assign(neq, 0);

	int i, n;

	// we need to fill the total displacement vector m_Ut
	// TODO: I need to find an easier way to do this
	FEMesh& mesh = m_fem.GetMesh();
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// displacement dofs
		n = node.m_ID[DOF_X]; if (n >= 0) m_Ut[n] = node.m_rt.x - node.m_r0.x;
		n = node.m_ID[DOF_Y]; if (n >= 0) m_Ut[n] = node.m_rt.y - node.m_r0.y;
		n = node.m_ID[DOF_Z]; if (n >= 0) m_Ut[n] = node.m_rt.z - node.m_r0.z;

		// rotational dofs
		n = node.m_ID[DOF_U]; if (n >= 0) m_Ut[n] = node.m_Dt.x - node.m_D0.x;
		n = node.m_ID[DOF_V]; if (n >= 0) m_Ut[n] = node.m_Dt.y - node.m_D0.y;
		n = node.m_ID[DOF_W]; if (n >= 0) m_Ut[n] = node.m_Dt.z - node.m_D0.z;
	}

	// initialize BFGS data
	m_bfgs.Init(neq, this, m_plinsolve);

	// set the create stiffness matrix flag
	m_breshape = true;

	return true;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FESolidSolver2::Serialize(DumpFile& ar)
{
	// Serialize parameters
	FESolver::Serialize(ar);
	
	if (ar.IsSaving())
	{
		ar << m_Dtol << m_Etol << m_Rtol << m_Rmin;
		ar << m_bsymm;
		ar << m_nrhs;
		ar << m_niter;
		ar << m_nref << m_ntotref;
		ar << m_naug;
		ar << m_neq << m_nreq;

		ar << m_bfgs.m_LStol << m_bfgs.m_LSiter << m_bfgs.m_LSmin;
		ar << m_bfgs.m_maxups;
		ar << m_bfgs.m_maxref;
		ar << m_bfgs.m_cmax;
		ar << m_bfgs.m_nups;
	}
	else
	{
		ar >> m_Dtol >> m_Etol >> m_Rtol >> m_Rmin;
		ar >> m_bsymm;
		ar >> m_nrhs;
		ar >> m_niter;
		ar >> m_nref >> m_ntotref;
		ar >> m_naug;
		ar >> m_neq >> m_nreq;

		ar >> m_bfgs.m_LStol >> m_bfgs.m_LSiter >> m_bfgs.m_LSmin;
		ar >> m_bfgs.m_maxups;
		ar >> m_bfgs.m_maxref;
		ar >> m_bfgs.m_cmax;
		ar >> m_bfgs.m_nups;
	}
}

//-----------------------------------------------------------------------------
//! Determine the number of linear equations and assign equation numbers
//!

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
bool FESolidSolver2::InitEquations()
{
	int i, j, n;

	// get the mesh
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
		for (i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(P[i]);
			for (j=0; j<(int)node.m_ID.size(); ++j)
				if (node.m_ID[j] >= 0) node.m_ID[j] = neq++;
		}
	}
	else
	{
		// give all free dofs an equation number
		for (i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			for (j=0; j<(int)node.m_ID.size(); ++j)
				if (node.m_ID[j] >= 0) node.m_ID[j] = neq++;
		}
	}

	// Next, we assign equation numbers to the rigid body degrees of freedom
	m_nreq = neq;
	int nrb = m_fem.Objects();
	for (i=0; i<nrb; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (j=0; j<6; ++j)
			if (RB.m_BC[j] >= 0)
			{
				RB.m_LM[j] = neq++;
			}
			else
			{
				RB.m_LM[j] = -1;
			}
	}

	// store the number of equations
	m_neq = neq;

	// we assign the rigid body equation number to
	// Also make sure that the nodes are NOT constrained!
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0)
		{
			FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(node.m_rid));
			node.m_ID[DOF_X] = -RB.m_LM[0]-2;
			node.m_ID[DOF_Y] = -RB.m_LM[1]-2;
			node.m_ID[DOF_Z] = -RB.m_LM[2]-2;
			node.m_ID[DOF_RU] = -RB.m_LM[3]-2;
			node.m_ID[DOF_RV] = -RB.m_LM[4]-2;
			node.m_ID[DOF_RW] = -RB.m_LM[5]-2;
		}
	}

	// adjust the rigid dofs that are prescribed
	for (i=0; i<nrb; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (j=0; j<6; ++j)
		{
			n = RB.m_LM[j];
			if (RB.m_BC[j] > 0) RB.m_LM[j] = -n-2;
		}
	}

	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
//!  Creates the global stiffness matrix
//! \todo Can we move this to the FEStiffnessMatrix::Create function?
bool FESolidSolver2::CreateStiffness(bool breset)
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
//!  This functions performs the Lagrange augmentations
//!  It returns true if all the augmentation have converged, 
//!	otherwise it returns false
//
//! \todo There is an inherent problem with this approach. Since
//!	      Lagrangian multipliers are inherited from previous timesteps
//!       they might not be zero in case a node-surface contact breaks. 
//!       The node's gap value needs to become negative to a certain value
//!       before the Lagr. multipliers dissapears. 
//
bool FESolidSolver2::Augment()
{
	// Assume we will pass (can't hurt to be optimistic)
	bool bconv = true;

	// Do contact augmentations
	if (m_fem.SurfacePairInteractions() > 0)
	{
		// loop over all contact interfaces
		for (int i=0; i<m_fem.SurfacePairInteractions(); ++i)
		{
			FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairInteraction(i));
			if (pci->IsActive()) bconv = (pci->Augment(m_naug) && bconv);
		}
	}

	// do nonlinear constraint augmentations
	int n = m_fem.NonlinearConstraints();
	for (int i=0; i<n; ++i) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		if (plc->IsActive()) bconv = plc->Augment(m_naug) && bconv;
	}

	// do incompressibility multipliers for 3Field domains
	FEMesh& mesh = m_fem.GetMesh();
	int ND = mesh.Domains();
	for (int i=0; i<ND; ++i)
	{
		FE3FieldElasticSolidDomain* pd = dynamic_cast<FE3FieldElasticSolidDomain*>(&mesh.Domain(i));
		if (pd) bconv = (pd->Augment() && bconv);
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FESolidSolver2::UpdateKinematics(vector<double>& ui)
{
	int i, n;

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// update rigid bodies
	UpdateRigidBodies(ui);

	// update flexible nodes
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// displacement dofs
		// current position = initial + total at prev conv step + total increment so far + current increment  
		if ((n = node.m_ID[DOF_X]) >= 0) node.m_rt.x = node.m_r0.x + m_Ut[n] + m_Ui[n] + ui[n];
		if ((n = node.m_ID[DOF_Y]) >= 0) node.m_rt.y = node.m_r0.y + m_Ut[n] + m_Ui[n] + ui[n];
		if ((n = node.m_ID[DOF_Z]) >= 0) node.m_rt.z = node.m_r0.z + m_Ut[n] + m_Ui[n] + ui[n];

		// rotational dofs
		if ((n = node.m_ID[DOF_U]) >= 0) node.m_Dt.x = node.m_D0.x + m_Ut[n] + m_Ui[n] + ui[n];
		if ((n = node.m_ID[DOF_V]) >= 0) node.m_Dt.y = node.m_D0.y + m_Ut[n] + m_Ui[n] + ui[n];
		if ((n = node.m_ID[DOF_W]) >= 0) node.m_Dt.z = node.m_D0.z + m_Ut[n] + m_Ui[n] + ui[n];
	}

	// make sure the prescribed displacements are fullfilled
	int ndis = m_fem.PrescribedBCs();
	for (i=0; i<ndis; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;

			FENode& node = mesh.Node(n);

			double g = r + s*m_fem.GetLoadCurve(lc)->Value();

			switch (bc)
			{
			case 0:
				node.m_rt.x = node.m_r0.x + g;
				break;
			case 1:
				node.m_rt.y = node.m_r0.y + g;
				break;
			case 2:
				node.m_rt.z = node.m_r0.z + g;
				break;
			case 20:
				{
					vec3d dr = node.m_r0;
					dr.x = 0; dr.unit(); dr *= g;

					node.m_rt.y = node.m_r0.y + dr.y;
					node.m_rt.z = node.m_r0.z + dr.z;
				}
				break;
			}
		}
	}

	// enforce the linear constraints
	// TODO: do we really have to do this? Shouldn't the algorithm
	// already guarantee that the linear constraints are satisfied?
	if (m_fem.m_LinC.size() > 0)
	{
		int nlin = m_fem.m_LinC.size();
		list<FELinearConstraint>::iterator it = m_fem.m_LinC.begin();
		double d;
		for (int n=0; n<nlin; ++n, ++it)
		{
			FELinearConstraint& lc = *it;
			FENode& node = mesh.Node(lc.master.node);

			d = 0;
			int ns = lc.slave.size();
			list<FELinearConstraint::SlaveDOF>::iterator si = lc.slave.begin();
			for (int i=0; i<ns; ++i, ++si)
			{
				FENode& node = mesh.Node(si->node);
				switch (si->bc)
				{
				case 0: d += si->val*(node.m_rt.x - node.m_r0.x); break;
				case 1: d += si->val*(node.m_rt.y - node.m_r0.y); break;
				case 2: d += si->val*(node.m_rt.z - node.m_r0.z); break;
				}
			}

			switch (lc.master.bc)
			{
			case 0: node.m_rt.x = node.m_r0.x + d; break;
			case 1: node.m_rt.y = node.m_r0.y + d; break;
			case 2: node.m_rt.z = node.m_r0.z + d; break;
			}
		}
	}


	// update velocity and accelerations
	// for dynamic simulations
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	if (pstep->m_nanalysis == FE_DYNAMIC)
	{
		int N = mesh.Nodes();
		double dt = pstep->m_dt;
		double a = 1.0 / (m_beta*dt);
		double b = a / dt;
		double c = 1.0 - 0.5/m_beta;
		for (i=0; i<N; ++i)
		{
			FENode& n = mesh.Node(i);
			n.m_at = (n.m_rt - n.m_rp)*b - n.m_vp*a + n.m_ap*c;
			n.m_vt = n.m_vp + (n.m_ap*(1.0 - m_gamma) + n.m_at*m_gamma)*dt;
		}
	}
}

//-----------------------------------------------------------------------------
//! Update DOF increments
void FESolidSolver2::UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap)
{
	int i, n;
    
	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();
    
	// update rigid bodies
	int nrb = m_fem.Objects();
	for (int i=0; i<nrb; ++i)
	{
		// get the rigid body
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.Object(i));
		if (RB.IsActive()) RB.UpdateIncrement(Ui, ui, emap);
	}
    
    
	// update flexible nodes
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
        
		// displacement dofs
		// current position = initial + total at prev conv step + total increment so far + current increment
		if ((n = node.m_ID[DOF_X]) >= 0) Ui[n] += ui[n];
		if ((n = node.m_ID[DOF_Y]) >= 0) Ui[n] += ui[n];
		if ((n = node.m_ID[DOF_Z]) >= 0) Ui[n] += ui[n];
        
		// rotational dofs
        // TODO: modify rotational update to account for exponential map or Cayley transform
		if ((n = node.m_ID[DOF_U]) >= 0) Ui[n] += ui[n];
		if ((n = node.m_ID[DOF_V]) >= 0) Ui[n] += ui[n];
		if ((n = node.m_ID[DOF_W]) >= 0) Ui[n] += ui[n];
	}
    
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FESolidSolver2::Update(vector<double>& ui)
{
	// update kinematics
	UpdateKinematics(ui);

	// update contact
	if (m_fem.SurfacePairInteractions() > 0) UpdateContact();

	// update element stresses
	UpdateStresses();

	// update other stuff that may depend on the deformation
	int NBL = m_fem.BodyLoads();
	for (int i=0; i<NBL; ++i)
	{
		FEPointBodyForce* pbf = dynamic_cast<FEPointBodyForce*>(m_fem.GetBodyLoad(i));
		if (pbf) pbf->Update();
	}

	// dump all states to the plot file when requested
	if (m_fem.GetCurrentStep()->GetPlotLevel() == FE_PLOT_MINOR_ITRS) m_fem.Write();
}

//-----------------------------------------------------------------------------
//! Updates the rigid body data
void FESolidSolver2::UpdateRigidBodies(vector<double>& ui)
{
	FEMesh& mesh = m_fem.GetMesh();

	// update rigid bodies
	int nrb = m_fem.Objects();
	for (int i=0; i<nrb; ++i)
	{
		// get the rigid body
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		if (RB.IsActive()) RB.Update(m_Ui, ui);
	}

	// update rigid joints
	int NC = m_fem.NonlinearConstraints();
	for (int i=0; i<NC; ++i)
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->Update();
	}
}

//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FESolidSolver2::UpdateStresses()
{
	FEMesh& mesh = m_fem.GetMesh();

	// update the stresses on all domains
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.UpdateStresses(m_fem);
	}
}

//-----------------------------------------------------------------------------
//! Update contact interfaces.
void FESolidSolver2::UpdateContact()
{
	// Update all contact interfaces
	for (int i=0; i<m_fem.SurfacePairInteractions(); ++i) 
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairInteraction(i));
		if (pci->IsActive()) pci->Update(m_niter);
	}
}

//-----------------------------------------------------------------------------
//! Update nonlinear constraints
void FESolidSolver2::UpdateConstraints()
{
	// Update all nonlinear constraints
	for (int i=0; i<m_fem.NonlinearConstraints(); ++i) 
	{
		FENLConstraint* pci = m_fem.NonlinearConstraint(i);
		if (pci->IsActive()) pci->Update();
	}
}

//-----------------------------------------------------------------------------
//!  This function mainly calls the Quasin routine 
//!  and deals with exceptions that require the immediate termination of
//!	quasi-Newton iterations.
bool FESolidSolver2::SolveStep(double time)
{
	bool bret;

	try
	{
		// let's try to call Quasin
		bret = Quasin(time);
	}
	catch (NegativeJacobian e)
	{
		// A negative jacobian was detected
		felog.printbox("ERROR","Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng+1, e.m_vol);
		if (m_fem.GetDebugFlag()) m_fem.Write();
		return false;
	}
	catch (MaxStiffnessReformations)
	{
		// max nr of reformations is reached
		felog.printbox("ERROR", "Max nr of reformations reached.");
		return false;
	}
	catch (ForceConversion)
	{
		// user forced conversion of problem
		felog.printbox("WARNING", "User forced conversion.\nSolution might not be stable.");
		return true;
	}
	catch (IterationFailure)
	{
		// user caused a forced iteration failure
		felog.printbox("WARNING", "User forced iteration failure.");
		return false;
	}
	catch (ZeroLinestepSize)
	{
		// a zero line step size was detected
		felog.printbox("ERROR", "Zero line step size.");
		return false;
	}
	catch (EnergyDiverging)
	{
		// problem was diverging after stiffness reformation
		felog.printbox("ERROR", "Problem diverging uncontrollably.");
		return false;
	}
	catch (FEMultiScaleException)
	{
		// the RVE problem didn't solve
		felog.printbox("ERROR", "The RVE problem has failed. Aborting macro run.");
		return false;
	}
	catch (DoRunningRestart)
	{
		// a request to fail the iteration and restart the time step
		if (m_fem.GetDebugFlag()) m_fem.Write();
		return false;
	}

	return bret;
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration. 
void FESolidSolver2::PrepStep(double time)
{
	// initialize counters
	m_niter = 0;	// nr of iterations
	m_nrhs  = 0;	// nr of RHS evaluations
	m_nref  = 0;	// nr of stiffness reformations
	m_ntotref = 0;
	m_bfgs.m_nups	= 0;	// nr of stiffness updates between reformations
	m_naug  = 0;	// nr of augmentations

	// zero total displacements
	zero(m_Ui);

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		ni.m_rp = ni.m_rt;
		ni.m_vp = ni.m_vt;
		ni.m_ap = ni.m_at;
		// ---> TODO: move to the FEPoroSoluteSolver
		for (int k=0; k<(int)ni.m_cp.size(); ++k) ni.m_cp[k] = ni.m_ct[k];
	}

	// apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
	NodalForces(m_Fn);

	// apply prescribed displacements
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_bfgs.m_ui;
	zero(ui);
	int neq = m_neq;
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

			double dq = r + s*m_fem.GetLoadCurve(lc)->Value();

			int I;

			FENode& node = m_fem.GetMesh().Node(n);

			switch (bc)
			{
				case DOF_X: 
					I = -node.m_ID[bc]-2;
					if (I>=0 && I<neq) 
						ui[I] = dq - (node.m_rt.x - node.m_r0.x);
					break;
				case DOF_Y: 
					I = -node.m_ID[bc]-2;
					if (I>=0 && I<neq) 
						ui[I] = dq - (node.m_rt.y - node.m_r0.y); 
					break;
				case DOF_Z: 
					I = -node.m_ID[bc]-2;
					if (I>=0 && I<neq) 
						ui[I] = dq - (node.m_rt.z - node.m_r0.z); 
					break;
					// ---> TODO: move to the FEPoroSolidSolver
				case DOF_P: 
					I = -node.m_ID[bc]-2;
					if (I>=0 && I<neq) 
						ui[I] = dq - node.m_pt; 
					break;
/*				case DOF_C: 
					I = -node.m_ID[bc]-2;
					if (I>=0 && I<neq) 
						ui[I] = dq - node.m_ct[0]; 
					break;
				case DOF_C+1: 
					I = -node.m_ID[bc]-2;
					if (I>=0 && I<neq) 
						ui[I] = dq - node.m_ct[1]; 
					break;*/
					// ---> TODO: change bc=20 to something else
				case 20:
				{
					vec3d dr = node.m_r0;
					dr.x = 0; dr.unit(); dr *= dq;
					
					I = -node.m_ID[1]-2;
					if (I>=0 && I<neq) 
						ui[I] = dr.y - (node.m_rt.y - node.m_r0.y); 
					I = -node.m_ID[2]-2;
					if (I>=0 && I<neq) 
						ui[I] = dr.z - (node.m_rt.z - node.m_r0.z); 
				}
					break;
				default:
					if ((bc >= DOF_C) && (bc < (int)node.m_ID.size())) {
						I = -node.m_ID[bc]-2;
						if (I>=0 && I<neq) 
							ui[I] = dq - node.m_ct[bc - DOF_C]; 
					}
			}
		}
	}

	// initialize rigid bodies
	int NO = m_fem.Objects();
	for (int i=0; i<NO; ++i) m_fem.Object(i)->Init();

	// calculate local rigid displacements
	for (int i=0; i<(int) m_fem.m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *m_fem.m_RDC[i];
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(DC.id));
		if (RB.IsActive() && DC.IsActive())
		{
			int I = DC.bc;
			int lc = DC.lc;
			if (lc >= 0)
			{
				RB.m_dul[I] = DC.ref + DC.sf*m_fem.GetLoadCurve(lc)->Value() - RB.m_Ut[DC.bc];
			}
		}
	}

	// calculate global rigid displacements
	for (int i=0; i<NO; ++i)
	{
		FERigidBody* prb = dynamic_cast<FERigidBody*>(m_fem.Object(i));
		if (prb)
		{
			FERigidBody& RB = *prb;
			if (RB.m_prb == 0)
			{
				for (int j=0; j<6; ++j) RB.m_du[j] = RB.m_dul[j];
			}
			else
			{
				double* dul = RB.m_dul;
				vec3d dr = vec3d(dul[0], dul[1], dul[2]);
				
				vec3d v = vec3d(dul[3], dul[4], dul[5]);
				double w = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
				quatd dq = quatd(w, v);

				FERigidBody* pprb = RB.m_prb;

				vec3d r0 = RB.m_rt;
				quatd Q0 = RB.m_qt;

				dr = Q0*dr;
				dq = Q0*dq*Q0.Inverse();

				while (pprb)
				{
					vec3d r1 = pprb->m_rt;
					dul = pprb->m_dul;

					quatd Q1 = pprb->m_qt;
					
					dr = r0 + dr - r1;

					// grab the parent's local displacements
					vec3d dR = vec3d(dul[0], dul[1], dul[2]);
					v = vec3d(dul[3], dul[4], dul[5]);
					w = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
					quatd dQ = quatd(w, v);

					dQ = Q1*dQ*Q1.Inverse();

					// update global displacements
					quatd Qi = Q1.Inverse();
					dr = dR + r1 + dQ*dr - r0;
					dq = dQ*dq;

					// move up in the chain
					pprb = pprb->m_prb;
					Q0 = Q1;
				}

				// set global displacements
				double* du = RB.m_du;

				du[0] = dr.x;
				du[1] = dr.y;
				du[2] = dr.z;

				v = dq.GetVector();
				w = dq.GetAngle();
				du[3] = w*v.x;
				du[4] = w*v.y;
				du[5] = w*v.z;
			}
		}
	}

	// store rigid displacements in Ui vector
	for (int i=0; i<NO; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (int j=0; j<6; ++j)
		{
			int I = -RB.m_LM[j]-2;
			if (I >= 0) ui[I] = RB.m_du[j];
		}
	}

	// apply prescribed rigid body forces
	// TODO: I don't think this does anything since
	//       the reaction forces are zeroed in the FESolidSolver::Residual function
	for (int i=0; i<(int) m_fem.m_RFC.size(); ++i)
	{
		FERigidBodyForce& FC = *m_fem.m_RFC[i];
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(FC.id));
		if (RB.IsActive() && FC.IsActive())
		{
			int lc = FC.lc;
			int I  = RB.m_LM[FC.bc];
			if ((I>=0) && (lc >= 0))
			{
				double f = m_fem.GetLoadCurve(lc)->Value()*FC.sf;
				m_Fn[I] += f;

				switch (FC.bc)
				{
				case 0: RB.m_Fr.x += f; break;
				case 1: RB.m_Fr.y += f; break;
				case 2: RB.m_Fr.z += f; break;
				case 3: RB.m_Mr.x += f; break;
				case 4: RB.m_Mr.y += f; break;
				case 5: RB.m_Mr.z += f; break;
				}
			}
		}
	}

	// initialize contact
	if (m_fem.SurfacePairInteractions() > 0) UpdateContact();

	// initialize nonlinear constraints
	if (m_fem.NonlinearConstraints() > 0) UpdateConstraints();

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	FEMaterialPoint::dt = m_fem.GetCurrentStep()->m_dt;
	FEMaterialPoint::time = m_fem.m_ftime;

	for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).InitElements();

	// update stresses
	UpdateStresses();
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
//! The details of this implementation of the BFGS method can be found in:
//!   "Finite Element Procedures", K.J. Bathe, p759 and following
bool FESolidSolver2::Quasin(double time)
{
	int i;

	vector<double> u0(m_neq);
	vector<double> Rold(m_neq);

	// convergence norms
	double	normR1;		// residual norm
	double	normE1;		// energy norm
	double	normU;		// displacement norm
	double	normu;		// displacement increment norm
	double	normRi;		// initial residual norm
	double	normEi;		// initial energy norm
	double	normEm;		// max energy norm
	double	normUi;		// initial displacement norm

	// initialize flags
	bool bconv = false;		// convergence flag
	bool breform = false;	// reformation flag
	bool sdflag = true;		// flag for steepest descent iterations in NLCG

	// Get the current step
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// prepare for the first iteration
	PrepStep(time);

	// do minor iterations callbacks
	m_fem.DoCallback(CB_MINOR_ITERS);

	// calculate initial stiffness matrix
	if (ReformStiffness() == false) return false;

	// calculate initial residual
	if (Residual(m_bfgs.m_R0) == false) return false;

	m_bfgs.m_R0 += m_Fd;

	// TODO: I can check here if the residual is zero.
	// If it is than there is probably no force acting on the system
	// if (m_R0*m_R0 < eps) bconv = true;

//	double r0 = m_R0*m_R0;

	felog.printf("\n===== beginning time step %d : %lg =====\n", pstep->m_ntimesteps+1, m_fem.m_ftime);

	// set the initial step length estimates to 1.0
	double s, olds, oldolds;  // line search step lengths from the current iteration and the two previous ones
	s=1; olds=1; oldolds=1;

	// loop until converged or when max nr of reformations reached
	do
	{
		Logfile::MODE oldmode = felog.GetMode();
		if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::FILE_ONLY);

		felog.printf(" %d\n", m_niter+1);
		felog.SetMode(oldmode);

		// assume we'll converge. 
		bconv = true;
		// solve the equations
		m_SolverTime.start();
		{
			m_bfgs.SolveEquations(m_bfgs.m_ui, m_bfgs.m_R0);
		}
		m_SolverTime.stop();

		// check for nans
		if (m_fem.GetDebugFlag())
		{
			double du = m_bfgs.m_ui*m_bfgs.m_ui;
			if (ISNAN(du)) throw NANDetected();
		}

		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_bfgs.m_R0*m_bfgs.m_R0);
			normEi = fabs(m_bfgs.m_ui*m_bfgs.m_R0);
			normUi = fabs(m_bfgs.m_ui*m_bfgs.m_ui);
			normEm = normEi;
		}

		// perform a linesearch
		// the geometry is also updated in the line search
		if (m_bfgs.m_LStol > 0) s = m_bfgs.LineSearch(1.0);
		else
		{
			s = 1;

			// Update geometry
			Update(m_bfgs.m_ui);

			// calculate residual at this point
			Residual(m_bfgs.m_R1);
		}

		// update total displacements
		int neq = (int)m_Ui.size();
        vector<double> ui(m_bfgs.m_ui);
        for (i=0; i<neq; ++i) ui[i] *= s;
        UpdateIncrements(m_Ui, ui, false);

		// calculate norms
		normR1 = m_bfgs.m_R1*m_bfgs.m_R1;
		normu  = (m_bfgs.m_ui*m_bfgs.m_ui)*(s*s);
		normU  = m_Ui*m_Ui;
		normE1 = s*fabs(m_bfgs.m_ui*m_bfgs.m_R1);

		// check residual norm
		if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;	

		// check displacement norm
		if ((m_Dtol > 0) && (normu  > (m_Dtol*m_Dtol)*normU )) bconv = false;

		// check energy norm
		if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;

		// check linestep size
		if ((m_bfgs.m_LStol > 0) && (s < m_bfgs.m_LSmin)) bconv = false;

		// check energy divergence
		if (normE1 > normEm) bconv = false;

		// print convergence summary
		oldmode = felog.GetMode();
		if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::FILE_ONLY);

		felog.printf(" Nonlinear solution status: time= %lg\n", time); 
		felog.printf("\tstiffness updates             = %d\n", m_bfgs.m_nups);
		felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_bfgs.m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
		felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		felog.printf("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
		felog.printf("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
		felog.printf("\t   displacement     %15le %15le %15le \n", normUi, normu ,(m_Dtol*m_Dtol)*normU );

		felog.SetMode(oldmode);

		// see if we may have a small residual
		if ((bconv == false) && (normR1 < m_Rmin))
		{
			// check for almost zero-residual on the first iteration
			// this might be an indication that there is no force on the system
			felog.printbox("WARNING", "No force acting on the system.");
			bconv = true;
		}

		// check if we have converged. 
		// If not, calculate the BFGS update vectors
		if (bconv == false)
		{
			if (s < m_bfgs.m_LSmin)
			{
				// check for zero linestep size
				felog.printbox("WARNING", "Zero linestep size. Stiffness matrix will now be reformed");
				breform = true;
			}
			else if (normE1 > normEm)
			{
				// check for diverging
				felog.printbox("WARNING", "Problem is diverging. Stiffness matrix will now be reformed");
				normEm = normE1;
				normEi = normE1;
				normRi = normR1;
				breform = true;
			}
			else
			{
				// If we havn't reached max nr of BFGS updates
				// do an update
				if (!breform)
				{
					if (m_bfgs.m_nups < m_bfgs.m_maxups-1)
					{
						if (m_bfgs.Update(s, m_bfgs.m_ui, m_bfgs.m_R0, m_bfgs.m_R1) == false)
						{
							// Stiffness update has failed.
							// this might be due a too large condition number
							// or the update was no longer positive definite.
							felog.printbox("WARNING", "The BFGS update has failed.\nStiffness matrix will now be reformed.");
							breform = true;
						}
					}
					else
					{
						// we've reached the max nr of BFGS updates, so
						// we need to do a stiffness reformation
						breform = true;

						// print a warning only if the user did not intent full-Newton
						if (m_bfgs.m_maxups > 0)
							felog.printbox("WARNING", "Max nr of iterations reached.\nStiffness matrix will now be reformed.");

					}
				}
			}	

			// zero displacement increments
			// we must set this to zero before the reformation
			// because we assume that the prescribed displacements are stored 
			// in the m_ui vector.
			zero(m_bfgs.m_ui);

			// reform stiffness matrices if necessary
			if (breform)
			{
				felog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);

				// reform the matrix
				if (ReformStiffness() == false) break;
	
				// reset reformation flag
				breform = false;
			}

			// copy last calculated residual
			m_bfgs.m_R0 = m_bfgs.m_R1;
		}
		else if (pstep->m_baugment)
		{
			// we have converged, so let's see if the augmentations have converged as well

			felog.printf("\n........................ augmentation # %d\n", m_naug+1);

			// do the augmentations
			bconv = Augment();

			// update counter
			++m_naug;

			// we reset the reformations counter
			m_nref = 0;
	
			// If we havn't converged we prepare for the next iteration
			if (!bconv) 
			{
				// Since the Lagrange multipliers have changed, we can't just copy 
				// the last residual but have to recalculate the residual
				// we also recalculate the stresses in case we are doing augmentations
				// for incompressible materials
				UpdateStresses();
				Residual(m_bfgs.m_R0);

				// reform the matrix if we are using full-Newton
				if (m_bfgs.m_maxups == 0)
				{
					felog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);
					if (ReformStiffness() == false) break;
				}
			}
		}
	
		// increase iteration number
		m_niter++;

		// let's flush the logfile to make sure the last output will not get lost
		felog.flush();

		// do minor iterations callbacks
		m_fem.DoCallback(CB_MINOR_ITERS);
	}
	while (bconv == false);

	// when converged, 
	// print a convergence summary to the felog file
	if (bconv)
	{
		Logfile::MODE mode = felog.SetMode(Logfile::FILE_ONLY);
		if (mode != Logfile::NEVER)
		{
			felog.printf("\nconvergence summary\n");
			felog.printf("    number of iterations   : %d\n", m_niter);
			felog.printf("    number of reformations : %d\n", m_nref);
		}
		felog.SetMode(mode);
	}

	// if converged we update the total displacements
	if (bconv)
	{
        UpdateIncrements(m_Ut, m_Ui, true);
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Reforms a stiffness matrix and factorizes it
bool FESolidSolver2::ReformStiffness()
{
	// first, let's make sure we have not reached the max nr of reformations allowed
	if (m_nref >= m_bfgs.m_maxref) throw MaxStiffnessReformations();

	// recalculate the shape of the stiffness matrix if necessary
	if (m_breshape)
	{
		// TODO: I don't think I need to update here
//		if (m_fem.m_bcontact) UpdateContact();

		// reshape the stiffness matrix
		if (!CreateStiffness(m_niter == 0)) return false;

		// reset reshape flag, except for contact
		m_breshape = (m_fem.SurfacePairInteractions() > 0? true : false);
	}

	// calculate the global stiffness matrix
	FETimePoint tp = m_fem.GetCurrentTime();
	bool bret = StiffnessMatrix(tp);

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
//! Calculates global stiffness matrix.

bool FESolidSolver2::StiffnessMatrix(const FETimePoint& tp)
{
	// get the stiffness matrix
	SparseMatrix& K = *m_pK;

	// zero stiffness matrix
	K.zero();

	// zero the residual adjustment vector
	zero(m_Fd);

	// nodal degrees of freedom
	int i, j, I;

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// calculate the stiffness matrix for each domain
	for (i=0; i<mesh.Domains(); ++i) 
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.StiffnessMatrix(this);
	}

	// calculate the body force stiffness matrix for each domain
	for (i=0; i<mesh.Domains(); ++i) 
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		int NBL = m_fem.BodyLoads();
		for (int j=0; j<NBL; ++j)
		{
			FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
			if (pbf) dom.BodyForceStiffness(this, *pbf);
		}
	}

	// Add mass matrix for dynamic problems
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	if (pstep->m_nanalysis == FE_DYNAMIC)
	{
		// scale factor
		double dt = tp.dt;
		double a = 1.0 / (m_beta*dt*dt);

		// loop over all domains
		for (int i=0; i<mesh.Domains(); ++i) 
		{
			FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
			dom.MassMatrix(this, a);
		}

		// loop over all rigid bodies
		for (int i=0; i<m_fem.Objects(); ++i)
		{
			FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
			RigidMassMatrix(RB);
		}
	}

	// calculate contact stiffness
	if (m_fem.SurfacePairInteractions() > 0) 
	{
		ContactStiffness();
	}

	// calculate stiffness matrices for surface loads
	int nsl = m_fem.SurfaceLoads();
	for (i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = m_fem.SurfaceLoad(i);

		// respect the pressure stiffness flag
		if ((dynamic_cast<FEPressureLoad*>(psl) == 0) || (m_fem.GetCurrentStep()->m_istiffpr != 0)) psl->StiffnessMatrix(this); 
	}

	// calculate nonlinear constraint stiffness
	// note that this is the contribution of the 
	// constrainst enforced with augmented lagrangian
	NonLinearConstraintStiffness();

	// point constraints
//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->StiffnessMatrix(this);

	// we still need to set the diagonal elements to 1
	// for the prescribed rigid body dofs.
	int NRB = m_fem.Objects();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& rb = static_cast<FERigidBody&>(*m_fem.Object(i));
		for (j=0; j<6; ++j)
			if (rb.m_LM[j] < -1)
			{
				I = -rb.m_LM[j]-2;
				K.set(I,I, 1);
			}
	}

	// let's check the stiffness matrix for zero diagonal elements
	if (m_fem.GetDebugFlag())
	{
		vector<int> zd;
		int neq = K.Size();
		for (i=0; i<neq; ++i)
		{
//			if (K.diag(i) == 0) zd.push_back(i);
			if (fabs(K.diag(i)) < 1e-15) zd.push_back(i);
		}

//		if (zd.empty() == false) throw ZeroDiagonal(zd, m_fem);
		if (zd.empty() == false) throw ZeroDiagonal(-1, -1);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FESolidSolver2::NonLinearConstraintStiffness()
{
	int N = m_fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->StiffnessMatrix(this);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact stiffness matrix

void FESolidSolver2::ContactStiffness()
{
	for (int i=0; i<m_fem.SurfacePairInteractions(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairInteraction(i));
		if (pci->IsActive()) pci->ContactStiffness(this);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contribution to the mass matrix from the rigid bodies
void FESolidSolver2::RigidMassMatrix(FERigidBody& RB)
{
	// element stiffness matrix
	vector<int> lm;
	matrix ke;
        
    // 6 dofs per rigid body
    ke.resize(6, 6);
    ke.zero();
        
    // Newmark integration rule
    double dt = m_fem.GetCurrentStep()->m_dt;
    double beta = m_beta;
    double gamma = m_gamma;
    double a = 1./(beta*dt*dt);
        
    // mass matrix
    double M = RB.m_mass*a;
        
    ke[0][0] = M;
    ke[1][1] = M;
    ke[2][2] = M;
        
    // evaluate mass moment of inertia at t
    mat3d Rt = RB.m_qt.RotationMatrix();
    mat3ds Jt = (Rt*RB.m_moi*Rt.transpose()).sym();
        
    // incremental rotation in spatial frame
    quatd q = RB.m_qt*RB.m_qp.Inverse();
    q.MakeUnit();                           // clean-up roundoff errors
    double theta = 2*tan(q.GetAngle()/2);   // get theta from Cayley transform
    vec3d e = q.GetVector();
        
    // skew-symmetric tensor whose axial vector is the incremental rotation
    mat3d qhat;
    qhat.skew(e*theta);
        
    // generate tensor T(theta)
    mat3d T = mat3dd(1) + qhat/2 + dyad(e*theta)/4;
        
    // skew-symmetric of angular momentum
    mat3d Jw;
    Jw.skew(Jt*RB.m_wt);
        
    // rotational inertia stiffness
    mat3d K = (Jt*T)*a*gamma - Jw/dt;
        
    ke[3][3] = K(0,0); ke[3][4] = K(0,1); ke[3][5] = K(0,2);
    ke[4][3] = K(1,0); ke[4][4] = K(1,1); ke[4][5] = K(1,2);
    ke[5][3] = K(2,0); ke[5][4] = K(2,1); ke[5][5] = K(2,2);
        
    lm.assign(RB.m_LM, RB.m_LM+6);
        
    AssembleStiffness(lm, ke);
}

//-----------------------------------------------------------------------------
//! This function calculates the rigid stiffness matrices

void FESolidSolver2::RigidStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{
	int i, j, k, l, n = en.size();
    
    // get nodal DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_NDOFS = fedofs.GetNDOFS();
    
	double Ri[3][3] = {0}, Rj[3][3] = {0};
	vector< vector<double> > kij; kij.assign(MAX_NDOFS, vector<double>(MAX_NDOFS));
    
	vector< vector<double> > KF; KF.assign(MAX_NDOFS, vector<double>(6));
	double KR[6][6];
    
	int *lmi, *lmj;
	int I, J;
    
	SparseMatrix& K = *m_pK;
    
	vec3d ai, aj;
    
	int ndof = ke.columns() / n;
    
	vector<double>& ui = m_bfgs.m_ui;
	FEMesh& mesh = m_fem.GetMesh();
    
	// loop over columns
	for (j=0; j<n; ++j)
	{
		FENode& nodej = mesh.Node(en[j]);
		if (nodej.m_rid >= 0)
		{
			// this is a rigid interface node
			// get the rigid body this node is attached to
			FERigidBody& RBj = dynamic_cast<FERigidBody&>(*m_fem.Object(nodej.m_rid));
            
			// get the rigid body equation nrs.
			lmj = RBj.m_LM;
            
			// get the relative distance to the center of mass
			aj = nodej.m_rt - RBj.m_rt;
            
			Rj[0][1] = aj.z; Rj[0][2] =-aj.y;
			Rj[1][0] =-aj.z; Rj[1][2] = aj.x;
			Rj[2][0] = aj.y; Rj[2][1] =-aj.x;
            
			// loop over rows
			for (i=0; i<n; ++i)
			{
				// get the element sub-matrix
				for (k=0; k<ndof; ++k)
					for (l=0; l<ndof; ++l)
						kij[k][l] = ke[ndof*i+k][ndof*j+l];
                
				FENode& nodei = mesh.Node(en[i]);
                
				if (nodei.m_rid>=0)
				{
					// node i is also a rigid body node
					// get the rigid body this node is attached to
					FERigidBody& RBi = dynamic_cast<FERigidBody&>(*m_fem.Object(nodei.m_rid));
                    
					lmi = RBi.m_LM;
					
					// get the relative distance (use alpha rule)
					ai = (nodei.m_rt - RBi.m_rt)*m_alpha + (nodei.m_rp - RBi.m_rp)*(1-m_alpha);
                    
					Ri[0][1] = ai.z; Ri[0][2] =-ai.y;
					Ri[1][0] =-ai.z; Ri[1][2] = ai.x;
					Ri[2][0] = ai.y; Ri[2][1] =-ai.x;
                    
					// Kij
					KR[0][0] = kij[0][0]; KR[0][1] = kij[0][1]; KR[0][2] = kij[0][2];
					KR[1][0] = kij[1][0]; KR[1][1] = kij[1][1]; KR[1][2] = kij[1][2];
					KR[2][0] = kij[2][0]; KR[2][1] = kij[2][1]; KR[2][2] = kij[2][2];
                    
                    
					//Kij*Rj
					KR[0][3] = kij[0][0]*Rj[0][0]+kij[0][1]*Rj[1][0]+kij[0][2]*Rj[2][0];
					KR[0][4] = kij[0][0]*Rj[0][1]+kij[0][1]*Rj[1][1]+kij[0][2]*Rj[2][1];
					KR[0][5] = kij[0][0]*Rj[0][2]+kij[0][1]*Rj[1][2]+kij[0][2]*Rj[2][2];
                    
					KR[1][3] = kij[1][0]*Rj[0][0]+kij[1][1]*Rj[1][0]+kij[1][2]*Rj[2][0];
					KR[1][4] = kij[1][0]*Rj[0][1]+kij[1][1]*Rj[1][1]+kij[1][2]*Rj[2][1];
					KR[1][5] = kij[1][0]*Rj[0][2]+kij[1][1]*Rj[1][2]+kij[1][2]*Rj[2][2];
                    
					KR[2][3] = kij[2][0]*Rj[0][0]+kij[2][1]*Rj[1][0]+kij[2][2]*Rj[2][0];
					KR[2][4] = kij[2][0]*Rj[0][1]+kij[2][1]*Rj[1][1]+kij[2][2]*Rj[2][1];
					KR[2][5] = kij[2][0]*Rj[0][2]+kij[2][1]*Rj[1][2]+kij[2][2]*Rj[2][2];
                    
                    
					// Ri^T*Kij
					KR[3][0] = Ri[0][0]*kij[0][0]+Ri[1][0]*kij[1][0]+Ri[2][0]*kij[2][0];
					KR[3][1] = Ri[0][0]*kij[0][1]+Ri[1][0]*kij[1][1]+Ri[2][0]*kij[2][1];
					KR[3][2] = Ri[0][0]*kij[0][2]+Ri[1][0]*kij[1][2]+Ri[2][0]*kij[2][2];
                    
					KR[4][0] = Ri[0][1]*kij[0][0]+Ri[1][1]*kij[1][0]+Ri[2][1]*kij[2][0];
					KR[4][1] = Ri[0][1]*kij[0][1]+Ri[1][1]*kij[1][1]+Ri[2][1]*kij[2][1];
					KR[4][2] = Ri[0][1]*kij[0][2]+Ri[1][1]*kij[1][2]+Ri[2][1]*kij[2][2];
                    
					KR[5][0] = Ri[0][2]*kij[0][0]+Ri[1][2]*kij[1][0]+Ri[2][2]*kij[2][0];
					KR[5][1] = Ri[0][2]*kij[0][1]+Ri[1][2]*kij[1][1]+Ri[2][2]*kij[2][1];
					KR[5][2] = Ri[0][2]*kij[0][2]+Ri[1][2]*kij[1][2]+Ri[2][2]*kij[2][2];
                    
                    
                    
					// Ri^T*Kij*Rj
					KR[3][3] = Ri[0][0]*KR[0][3]+Ri[1][0]*KR[1][3]+Ri[2][0]*KR[2][3];
					KR[3][4] = Ri[0][0]*KR[0][4]+Ri[1][0]*KR[1][4]+Ri[2][0]*KR[2][4];
					KR[3][5] = Ri[0][0]*KR[0][5]+Ri[1][0]*KR[1][5]+Ri[2][0]*KR[2][5];
                    
					KR[4][3] = Ri[0][1]*KR[0][3]+Ri[1][1]*KR[1][3]+Ri[2][1]*KR[2][3];
					KR[4][4] = Ri[0][1]*KR[0][4]+Ri[1][1]*KR[1][4]+Ri[2][1]*KR[2][4];
					KR[4][5] = Ri[0][1]*KR[0][5]+Ri[1][1]*KR[1][5]+Ri[2][1]*KR[2][5];
                    
					KR[5][3] = Ri[0][2]*KR[0][3]+Ri[1][2]*KR[1][3]+Ri[2][2]*KR[2][3];
					KR[5][4] = Ri[0][2]*KR[0][4]+Ri[1][2]*KR[1][4]+Ri[2][2]*KR[2][4];
					KR[5][5] = Ri[0][2]*KR[0][5]+Ri[1][2]*KR[1][5]+Ri[2][2]*KR[2][5];
                    
					// add the stiffness components to the Krr matrix
					for (k=0; k<6; ++k)
						for (l=0; l<6; ++l)
						{
							J = lmj[k];
							I = lmi[l];
                            
							if (I >= 0)
							{
                                // multiply KR by alpha for alpha rule
								if (J < -1) m_Fd[I] -= KR[l][k]*m_alpha*ui[-J-2];
								else if (J >= 0) K.add(I,J, KR[l][k]*m_alpha);
							}
						}
                    
					// we still need to couple the non-rigid degrees of node i to the
					// rigid dofs of node j
					for (k=3; k<ndof; ++k)
						for (l=0; l<3; ++l)
						{
							KF[k][l] = kij[k][l];
							KF[k][3+l] = kij[k][0]*Rj[0][l] + kij[k][1]*Rj[1][l] + kij[k][2]*Rj[2][l];
						}
                    
					for (k=0; k<6; ++k)
						for (l=3; l<ndof; ++l)
						{
							J = lmj[k];
							I = elm[ndof*i+l];
                            
							if (I >= 0)
							{
                                // multiply KF by alpha for alpha rule
								if (J < -1) m_Fd[I] -= KF[l][k]*m_alpha*ui[-J-2];
								else if (J >= 0) K.add(I,J, KF[l][k]*m_alpha);
							}
						}
                    
                    // now the transpose location
					for (k=0; k<3; ++k)
						for (l=3; l<ndof; ++l)
						{
							KF[l][k] = kij[k][l];
							KF[l][3+k] = kij[0][l]*Rj[0][k] + kij[1][l]*Rj[1][k] + kij[2][l]*Rj[2][k];
						}
                    
					for (k=0; k<6; ++k)
						for (l=3; l<ndof; ++l)
						{
							J = elm[ndof*j+l];
							I = lmi[k];
                            
							if (I >= 0)
							{
								if (J < -1) m_Fd[I] -= KF[l][k]*ui[-J-2];
								else if (J >= 0) K.add(I,J, KF[l][k]);
							}
						}
                    
				}
				else
				{
					// node i is not a rigid body node
					// add the stiffness components to the Kfr matrix
                    
					// Kij
					for (k=0; k<ndof; ++k)
						for (l=0; l<3; ++l)
						{
							KF[k][l] = kij[k][l];
							KF[k][3+l] = kij[k][0]*Rj[0][l] + kij[k][1]*Rj[1][l] + kij[k][2]*Rj[2][l];
						}
                    
					for (k=0; k<6; ++k)
						for (l=0; l<ndof; ++l)
						{
							J = lmj[k];
							I = elm[ndof*i+l];
                            
							if (I >= 0)
							{
                                // multiply KF by alpha for alpha rule
								if (J < -1) m_Fd[I] -= KF[l][k]*m_alpha*ui[-J-2];
								else if (J >= 0) K.add(I,J, KF[l][k]*m_alpha);
							}
						}
				}
			}
		}
		else
		{
			// loop over rows
			for (i=0; i<n; ++i)
			{
				FENode& nodei = mesh.Node(en[i]);
				if (nodei.m_rid>=0)
				{
					// node i is a rigid body
					// get the rigid body this node is attached to
					FERigidBody& RBi = dynamic_cast<FERigidBody&>(*m_fem.Object(nodei.m_rid));
                    
					// get the rigid body equation nrs.
					lmi = RBi.m_LM;
                    
					// get the relative distance (use alpha rule)
					ai = (nodei.m_rt - RBi.m_rt)*m_alpha + (nodei.m_rp - RBi.m_rp)*(1-m_alpha);
                    
					Ri[0][1] = ai.z; Ri[0][2] =-ai.y;
					Ri[1][0] =-ai.z; Ri[1][2] = ai.x;
					Ri[2][0] = ai.y; Ri[2][1] =-ai.x;
                    
					// get the element sub-matrix
					for (k=0; k<ndof; ++k)
						for (l=0; l<ndof; ++l)
							kij[k][l] = ke[ndof*i+k][ndof*j+l];
                    
					// add the stiffness components to the Krf matrix
                    
					// Kij
					for (k=0; k<ndof; ++k)
						for (l=0; l<3; ++l)
						{
							KF[k][l] = kij[l][k];
							KF[k][3+l] = Ri[0][l]*kij[0][k] + Ri[1][l]*kij[1][k] + Ri[2][l]*kij[2][k];
						}
                    
					for (k=0; k<6; ++k)
						for (l=0; l<ndof; ++l)
						{
							I = lmi[k];
							J = elm[ndof*j+l];
                            
							if (I >= 0)
							{
								if (J < -1) m_Fd[I] -= KF[l][k]*ui[-J-2];
								else if (J >= 0) K.add(I,J, KF[l][k]);
							}
						}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! \todo This function is only used for rigid joints. I need to figure out if
//!       I can use the other assembly function.
void FESolidSolver2::AssembleStiffness(std::vector<int>& lm, matrix& ke)
{
	m_pK->Assemble(ke, lm);
}

//-----------------------------------------------------------------------------
//!  Assembles the element stiffness matrix into the global stiffness matrix.
//!  Also adjusts the global stiffness matrix and residual to take the 
//!  prescribed displacements into account.

//! \todo In stead of changing the global stiffness matrix to accomodate for 
//!       the rigid bodies and linear constraints, can I modify the element stiffness
//!       matrix prior to assembly? I might have to change the elm vector as well as 
//!       the element matrix size.

void FESolidSolver2::AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{
    // get nodal DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_NDOFS = fedofs.GetNDOFS();

	// assemble into global stiffness matrix
	m_pK->Assemble(ke, elm);

	vector<double>& ui = m_bfgs.m_ui;

	// adjust for linear constraints
	if (m_fem.m_LinC.size() > 0)
	{
		int i, j, l;
		int nlin = m_fem.m_LinC.size();

		int ndof = ke.rows();
		int ndn = ndof / en.size();

		SparseMatrix& K = *m_pK;

		// loop over all stiffness components 
		// and correct for linear constraints
		int ni, nj, li, lj, I, J, k;
		double kij;
		for (i=0; i<ndof; ++i)
		{
			ni = MAX_NDOFS*(en[i/ndn]) + i%ndn;
			li = m_fem.m_LCT[ni];
			for (j=0; j<ndof; ++j)
			{
				nj = MAX_NDOFS*(en[j/ndn]) + j%ndn;
				lj = m_fem.m_LCT[nj];

				if ((li >= 0) && (lj < 0))
				{
					// dof i is constrained
					FELinearConstraint& Li = *m_fem.m_LCA[li];

					assert(elm[i] == -1);

					list<FELinearConstraint::SlaveDOF>::iterator is = Li.slave.begin();
					for (k=0; k<(int)Li.slave.size(); ++k, ++is)
					{
						I = is->neq;
						J = elm[j];
						kij = is->val*ke[i][j];
						if ((J>=I) && (I >=0)) K.add(I,J, kij);
						else
						{
							// adjust for prescribed dofs
							J = -J-2;
							if ((J>=0) && (J<m_nreq) && (I>=0)) m_Fd[I] -= kij*ui[J];
						}
					}
				}
				else if ((lj >= 0) && (li < 0))
				{
					// dof j is constrained
					FELinearConstraint& Lj = *m_fem.m_LCA[lj];

					assert(elm[j] == -1);

					list<FELinearConstraint::SlaveDOF>::iterator js = Lj.slave.begin();

					for (k=0; k<(int)Lj.slave.size(); ++k, ++js)
					{
						I = elm[i];
						J = js->neq;
						kij = js->val*ke[i][j];
						if ((J>=I) && (I >=0)) K.add(I,J, kij);
						else
						{
							// adjust for prescribed dofs
							J = -J-2;
							if ((J>=0) && (J<m_nreq) && (I>=0)) m_Fd[I] -= kij*ui[J];
						}
					}
				}
				else if ((li >= 0) && (lj >= 0))
				{
					// both dof i and j are constrained
					FELinearConstraint& Li = *m_fem.m_LCA[li];
					FELinearConstraint& Lj = *m_fem.m_LCA[lj];

					list<FELinearConstraint::SlaveDOF>::iterator is = Li.slave.begin();
					list<FELinearConstraint::SlaveDOF>::iterator js = Lj.slave.begin();

					assert(elm[i] == -1);
					assert(elm[j] == -1);

					for (k=0; k<(int)Li.slave.size(); ++k, ++is)
					{
						js = Lj.slave.begin();
						for  (l=0; l<(int)Lj.slave.size(); ++l, ++js)
						{
							I = is->neq;
							J = js->neq;
							kij = ke[i][j]*is->val*js->val;

							if ((J>=I) && (I >=0)) K.add(I,J, kij);
							else
							{
								// adjust for prescribed dofs
								J = -J-2;
								if ((J>=0) && (J<m_nreq) && (I>=0)) m_Fd[I] -= kij*ui[J];
							}
						}
					}
				}
			}
		}
	}

	// adjust stiffness matrix for prescribed degrees of freedom
	// NOTE: I had to comment this if statement out since otherwise
	//       poroelastic DOF's that are set as free-draining in the
	//       sliding2 contact code are skipt and zeroes will appear
	//       on the diagonal of the stiffness matrix.
//	if (m_fem.m_DC.size() > 0)
	{
		int i, j;
		int I, J;

		SparseMatrix& K = *m_pK;

		int N = ke.rows();

		// loop over columns
		for (j=0; j<N; ++j)
		{
			J = -elm[j]-2;
			if ((J >= 0) && (J<m_nreq))
			{
				// dof j is a prescribed degree of freedom

				// loop over rows
				for (i=0; i<N; ++i)
				{
					I = elm[i];
					if (I >= 0)
					{
						// dof i is not a prescribed degree of freedom
						m_Fd[I] -= ke[i][j]*ui[J];
					}
				}

				// set the diagonal element of K to 1
				K.set(J,J, 1);			
			}
		}
	}

	// see if there are any rigid body dofs here
	if (m_fem.Objects()) RigidStiffness(en, elm, ke);
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FESolidSolver2::ContactForces(FEGlobalVector& R)
{
	for (int i=0; i<m_fem.SurfacePairInteractions(); ++i) 
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairInteraction(i));
		if (pci->IsActive()) pci->ContactForces(R);
	}
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FESolidSolver2::Residual(vector<double>& R)
{
	int i;
	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// setup the global vector
	FEResidualVector RHS(GetFEModel(), R, m_Fr);

	// zero rigid body reaction forces
	int NRB = m_fem.Objects();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);
	}

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// calculate the internal (stress) forces
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// update body forces
	for (i=0; i<m_fem.BodyLoads(); ++i)
	{
		// TODO: I don't like this but for now I'll hard-code the modification of the
		//       force center position
		FEPointBodyForce* pbf = dynamic_cast<FEPointBodyForce*>(m_fem.GetBodyLoad(i));
		if (pbf)
		{
			if (pbf->m_rlc[0] >= 0) pbf->m_rc.x = m_fem.GetLoadCurve(pbf->m_rlc[0])->Value();
			if (pbf->m_rlc[1] >= 0) pbf->m_rc.y = m_fem.GetLoadCurve(pbf->m_rlc[1])->Value();
			if (pbf->m_rlc[2] >= 0) pbf->m_rc.z = m_fem.GetLoadCurve(pbf->m_rlc[2])->Value();
		}
	}

	// calculate the body forces
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		for (int j=0; j<m_fem.BodyLoads(); ++j)
		{
			FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
			dom.BodyForce(RHS, *pbf);
		}
	}

	// calculate inertial forces for dynamic problems
	if (m_fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC) InertialForces(RHS);

	// calculate forces due to surface loads
	int nsl = m_fem.SurfaceLoads();
	for (i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->Residual(RHS);
	}

	// calculate contact forces
	if (m_fem.SurfacePairInteractions() > 0)
	{
		ContactForces(RHS);
	}

	// calculate nonlinear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	NonLinearConstraintForces(RHS);

	// forces due to point constraints
//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->Residual(this, R);

	// set the nodal reaction forces
	// TODO: Is this a good place to do this?
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		node.m_Fr = vec3d(0,0,0);

		int n;
		if ((n = -node.m_ID[DOF_X]-2) >= 0) node.m_Fr.x = -m_Fr[n];
		if ((n = -node.m_ID[DOF_Y]-2) >= 0) node.m_Fr.y = -m_Fr[n];
		if ((n = -node.m_ID[DOF_Z]-2) >= 0) node.m_Fr.z = -m_Fr[n];
	}

	// increase RHS counter
	m_nrhs++;

	return true;
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces 
void FESolidSolver2::NonLinearConstraintForces(FEGlobalVector& R)
{
	int N = m_fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->Residual(R);
	}
}

//-----------------------------------------------------------------------------
//! calculates the concentrated nodal forces

void FESolidSolver2::NodalForces(vector<double>& F)
{
	int i, id, bc, lc, n;
	double s, f;
	vec3d a;
	int* lm;

	// zero nodal force vector
	zero(F);

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
			s    = fc.s;		// force scale factor

			FENode& node = mesh.Node(id);

			n = node.m_ID[bc];
		
			f = s*m_fem.GetLoadCurve(lc)->Value();
			
			// For pressure and concentration loads, multiply by dt
			// for consistency with evaluation of residual and stiffness matrix
			if ((bc == DOF_P) || (bc >= DOF_C))
				f *= m_fem.GetCurrentStep()->m_dt;

			if (n >= 0) F[n] = f;
			else if (node.m_rid >=0)
			{
				// this is a rigid body node
				FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(node.m_rid));

				// get the relative position
				a = node.m_rt - RB.m_rt;

				lm = RB.m_LM;
				switch (bc)
				{
				case 0:
					if (lm[0] >= 0) F[lm[0]] +=  f;
					if (lm[4] >= 0) F[lm[4]] +=  a.z*f;
					if (lm[5] >= 0) F[lm[5]] += -a.y*f;
					break;
				case 1:
					if (lm[1] >= 0) F[lm[1]] +=  f;
					if (lm[3] >= 0) F[lm[3]] += -a.z*f;
					if (lm[5] >= 0) F[lm[5]] +=  a.x*f;
					break;
				case 2:
					if (lm[2] >= 0) F[lm[2]] +=  f;
					if (lm[3] >= 0) F[lm[3]] +=  a.y*f;
					if (lm[4] >= 0) F[lm[4]] += -a.x*f;
					break;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the inertial forces for dynamic problems

void FESolidSolver2::InertialForces(FEGlobalVector& R)
{
	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();
    FEModel& fem = GetFEModel();
    
	// allocate F
	vector<double> F(3*mesh.Nodes());
	zero(F);
    
	// calculate F
	double dt = m_fem.GetCurrentStep()->m_dt;
    // Newmark rule
    double a = 1.0 / (m_beta*dt);
    double b = a / dt;
    double c = 1.0 - 0.5/m_beta;
    // update kinematics since acceleration is needed for inertial forces
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
        node.m_at = (node.m_rt - node.m_rp)*b - node.m_vp*a + node.m_ap*c;
        node.m_vt = node.m_vp + (node.m_ap*(1.0 - m_gamma) + node.m_at*m_gamma)*dt;
	}
    
	// update rigid bodies
	int nrb = m_fem.Objects();
	for (int i=0; i<nrb; ++i)
	{
		// get the rigid body
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		if (RB.IsActive())
        {
            // acceleration and velocity of center of mass
            RB.m_at = (RB.m_rt - RB.m_rp)*b - RB.m_vp*a + RB.m_ap*c;
            RB.m_vt = RB.m_vp + (RB.m_ap*(1.0 - m_gamma) + RB.m_at*m_gamma)*dt;
            // angular acceleration and velocity of rigid body
            quatd q = RB.m_qt*RB.m_qp.Inverse();
            q.MakeUnit();
            vec3d vq = q.GetVector()*(2*tan(q.GetAngle()/2));  // Cayley transform
            RB.m_wt = vq*(a*m_gamma) - RB.m_wp + (RB.m_wp + RB.m_alp*dt/2.)*(2-m_gamma/m_beta);
            q.RotateVector(RB.m_wt);
            RB.m_alt = vq*b - RB.m_wp*a + RB.m_alp*c;
            q.RotateVector(RB.m_alt);
        }
	}
    
	// calculate the interial forces for all elastic domains (except rigid domains)
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(nd));
		dom.InertialForces(R, F);
	}

	// calculate rigid body inertial forces
	for (int i=0; i<m_fem.Objects(); ++i)
	{
		FERigidBody& RB = static_cast<FERigidBody&>(*m_fem.Object(i));
		if (RB.IsActive())
        {
			RigidInertialForces(RB, R);
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the inertial force contribution from a rigid body and assemble
//! it in the global vector R.
void FESolidSolver2::RigidInertialForces(FERigidBody& RB, FEGlobalVector& R)
{
    // 6 dofs per rigid body
    double fe[6];
        
    // rate of change of linear momentum = mass*acceleration
    vec3d F = RB.m_at*RB.m_mass;
        
    fe[0] = -F.x;
    fe[1] = -F.y;
    fe[2] = -F.z;
        
    double dt = m_fem.GetCurrentStep()->m_dt;
        
    // evaluate mass moment of inertia at t and tp
    mat3d Rt = RB.m_qt.RotationMatrix();
    mat3ds Jt = (Rt*RB.m_moi*Rt.transpose()).sym();
    mat3d Rp = RB.m_qp.RotationMatrix();
    mat3ds Jp = (Rp*RB.m_moi*Rp.transpose()).sym();
        
    // evaluate rate of change of angular momentum
    vec3d M = (Jt*RB.m_wt - Jp*RB.m_wp)/dt;
        
    fe[3] = -M.x;
    fe[4] = -M.y;
    fe[5] = -M.z;
        
    R.Assemble(RB.m_LM, fe, 6);
        
    // add to rigid body force
    RB.m_Fr += F;
        
    // add to rigid body torque
    RB.m_Mr += M;
}
