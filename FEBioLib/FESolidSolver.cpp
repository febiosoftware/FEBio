#include "stdafx.h"
#include "FESolidSolver.h"
#include "FECore/FENodeReorder.h"
#include "FERigid.h"
#include "FERigidBody.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FEPointBodyForce.h"
#include "FEPressureLoad.h"
#include "log.h"

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
//! FESolidSolver Construction
//
FESolidSolver::FESolidSolver(FEModel& fem) : FESolver(fem)
{
	// default values
	m_Rtol = 0;	// deactivate residual convergence 
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_Rmin = 1.0e-20;

	m_niter = 0;
	m_nreq = 0;
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FESolidSolver
//
bool FESolidSolver::Init()
{
	// initialize base class
	if (FESolver::Init() == false) return false;

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
	FEMesh& mesh = m_fem.m_mesh;
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

	// calculate the inverse mass vector for explicit analysis
	if (m_solvertype==2)
	{
		matrix ke;
		int i, j, iel, n;
		int nint, neln;
		double *H, kab;
		vector <int> lm;
		vector <double> el_lumped_mass;
		for (int nd = 0; nd < mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
			if (pbd)
			{
				for (iel=0; iel<pbd->Elements(); ++iel)
				{
					FESolidElement& el = pbd->Element(iel);
					pbd->UnpackLM(el, lm);

					FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(m_fem.GetMaterial(el.GetMatID()));

					double d = pme->Density();

					nint = el.GaussPoints();
					neln = el.Nodes();

					ke.resize(3*neln, 3*neln);
					ke.zero();
					el_lumped_mass.resize(3*neln);
					for (i=0; i<3*neln; ++i) el_lumped_mass[i]=0.0;

					// create the element mass matrix
					for (n=0; n<nint; ++n)
					{
						double detJ0 = pbd->detJ0(el, n)*el.GaussWeights()[n];

						H = el.H(n);
						for (i=0; i<neln; ++i)
							for (j=0; j<neln; ++j)
							{
								kab = H[i]*H[j]*detJ0*d;
								ke[3*i  ][3*j  ] += kab;
								ke[3*i+1][3*j+1] += kab;
								ke[3*i+2][3*j+2] += kab;
							}	
					}
					// reduce to a lumped mass vector
					for (i=0; i<3*neln; ++i)
							for (j=0; j<neln; ++j)
							{
								el_lumped_mass[i]+=ke[i][j];
							}	
					
					// assemble element matrix into inv_mass vector
					AssembleResidual(el.m_node, lm, el_lumped_mass, m_inv_mass);
				}
			}
		}
	}
	// set the create stiffness matrix flag
	m_breshape = true;

	return true;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FESolidSolver::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_Dtol << m_Etol << m_Rtol << m_Rmin;
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
bool FESolidSolver::InitEquations()
{
	int i, j, n;

	// get the mesh
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
		for (i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(P[i]);
			for (j=0; j<MAX_NDOFS; ++j)
				if (node.m_ID[j] >= 0) node.m_ID[j] = neq++;
		}
	}
	else
	{
		// give all free dofs an equation number
		for (i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			for (j=0; j<MAX_NDOFS; ++j)
				if (node.m_ID[j] >= 0) node.m_ID[j] = neq++;
		}
	}

	// Next, we assign equation numbers to the rigid body degrees of freedom
	m_nreq = neq;
	int nrb = m_fem.m_Obj.size();
	for (i=0; i<nrb; ++i)
	{
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[i]);
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_fem.GetMaterial(RB.m_mat));
		assert(pm);
		for (j=0; j<6; ++j)
			if (pm->m_bc[j] >= 0)
			{
				RB.m_LM[j] = neq++;
			}
			else 
				RB.m_LM[j] = -1;
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
			FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[node.m_rid]);
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
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[i]);
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_fem.GetMaterial(RB.m_mat));
		for (j=0; j<6; ++j)
		{
			n = RB.m_LM[j];
			if (pm->m_bc[j] > 0) RB.m_LM[j] = -n-2;
		}
	}

	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
//!  This functions performs the Lagrange augmentations
//!  It returns true if all the augmentation have converged, 
//!	otherwise it returns false
//
// TODO: There is an inherent problem with this approach. Since
//  Lagrangian multipliers are inherited from previous timesteps
//  they might not be zero in case a node-surface contact breaks. 
//  The node's gap value needs to become negative to a certain value
//  before the Lagr. multipliers dissapears. 
//
bool FESolidSolver::Augment()
{
	// Assume we will pass (can't hurt to be optimistic)
	bool bconv = true;

	// Do contact augmentations
	if (m_fem.ContactInterfaces() > 0)
	{
		// loop over all contact interfaces
		for (int i=0; i<m_fem.ContactInterfaces(); ++i) bconv = m_fem.ContactInterface(i)->Augment(m_naug) && bconv;
	}

	// do nonlinear constraint augmentations
	int n = m_fem.NonlinearConstraints();
	for (int i=0; i<n; ++i) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		bconv = plc->Augment(m_naug) && bconv;
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
void FESolidSolver::UpdateKinematics(vector<double>& ui)
{
	int i, n;

	// get the mesh
	FEMesh& mesh = m_fem.m_mesh;

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
	int ndis = m_fem.m_DC.size();
	for (i=0; i<ndis; ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;	// GAA

			FENode& node = mesh.Node(n);

			double g = r + s*m_fem.GetLoadCurve(lc)->Value(); // GAA

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
	if ((pstep->m_nanalysis == FE_DYNAMIC)&&(m_solvertype==0))
	{
		int N = mesh.Nodes();
		double dt = pstep->m_dt;
		double a = 4.0 / dt;
		double b = a / dt;
		for (i=0; i<N; ++i)
		{
			FENode& n = mesh.Node(i);
			n.m_at = (n.m_rt - n.m_rp)*b - n.m_vp*a - n.m_ap;
			n.m_vt = n.m_vp + (n.m_ap + n.m_at)*dt*0.5;
		}
	}
	//	// update velocity - acceleration is calculated in Quasin
	//{
	//	int N = mesh.Nodes();
	//	double dt = m_fem.m_pStep->m_dt;
	//	for (i=0; i<N; ++i)
	//	{
	//		FENode& n = mesh.Node(i);
	//	}
	//}
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FESolidSolver::Update(vector<double>& ui)
{
	// update kinematics
	UpdateKinematics(ui);

	// update contact
	if (m_fem.ContactInterfaces() > 0) UpdateContact();

	// update element stresses
	UpdateStresses();

	// update other stuff that may depend on the deformation
	int NBF = m_fem.BodyForces();
	for (int i=0; i<NBF; ++i)
	{
		FEPointBodyForce* pbf = dynamic_cast<FEPointBodyForce*>(m_fem.GetBodyForce(i));
		if (pbf) pbf->Update();
	}

	// dump all states to the plot file when requested
	if (m_fem.GetCurrentStep()->m_nplot == FE_PLOT_MINOR_ITRS) m_fem.Write();
}

//-----------------------------------------------------------------------------
//! Updates the rigid body data
void FESolidSolver::UpdateRigidBodies(vector<double>& ui)
{
	FEMesh& mesh = m_fem.m_mesh;

	// update rigid bodies
	int nrb = m_fem.m_Obj.size();
	for (int i=0; i<nrb; ++i)
	{
		// get the rigid body
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[i]);
		if (RB.IsActive()) RB.Update(m_Ui, ui);
	}

	// update rigid joints
	int NC = m_fem.NonlinearConstraints();
	for (int i=0; i<NC; ++i)
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		plc->Update();
	}
}

//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FESolidSolver::UpdateStresses()
{
	FEMesh& mesh = m_fem.m_mesh;

	// update the stresses on all domains
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.UpdateStresses(m_fem);
	}
}

//-----------------------------------------------------------------------------
//! Update contact interfaces.
void FESolidSolver::UpdateContact()
{
	// Update all contact interfaces
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.ContactInterface(i)->Update(m_niter);
}

//-----------------------------------------------------------------------------
//!  This function mainly calls the Quasin routine 
//!  and deals with exceptions that require the immediate termination of
//!	quasi-Newton iterations.
bool FESolidSolver::SolveStep(double time)
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
		clog.printbox("ERROR","Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng, e.m_vol);
		if (m_fem.GetDebugFlag()) m_fem.Write();
		return false;
	}
	catch (MaxStiffnessReformations)
	{
		// max nr of reformations is reached
		clog.printbox("ERROR", "Max nr of reformations reached.");
		return false;
	}
	catch (ForceConversion)
	{
		// user forced conversion of problem
		clog.printbox("WARNING", "User forced conversion.\nSolution might not be stable.");
		return true;
	}
	catch (IterationFailure)
	{
		// user caused a forced iteration failure
		clog.printbox("WARNING", "User forced iteration failure.");
		return false;
	}
	catch (ZeroLinestepSize)
	{
		// a zero line step size was detected
		clog.printbox("ERROR", "Zero line step size.");
		return false;
	}
	catch (EnergyDiverging)
	{
		// problem was diverging after stiffness reformation
		clog.printbox("ERROR", "Problem diverging uncontrollably.");
		return false;
	}
	catch (FEMultiScaleException)
	{
		// the RVE problem didn't solve
		clog.printbox("ERROR", "The RVE problem has failed. Aborting macro run.");
		return false;
	}

	return bret;
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration. 
void FESolidSolver::PrepStep(double time)
{
	int i, j;

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
	for (i=0; i<m_fem.m_mesh.Nodes(); ++i)
	{
		m_fem.m_mesh.Node(i).m_rp = m_fem.m_mesh.Node(i).m_rt;
		m_fem.m_mesh.Node(i).m_vp = m_fem.m_mesh.Node(i).m_vt;
		m_fem.m_mesh.Node(i).m_ap = m_fem.m_mesh.Node(i).m_at;
		// ---> TODO: move to the FEPoroSoluteSolver
		for (int k=0; k<MAX_CDOFS; ++k)
			m_fem.m_mesh.Node(i).m_cp[k] = m_fem.m_mesh.Node(i).m_ct[k];
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
	for (i=0; i<(int) m_fem.m_DC.size(); ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;	// GAA

			double dq = r + s*m_fem.GetLoadCurve(lc)->Value();	// GAA

			int I;

			FENode& node = m_fem.m_mesh.Node(n);

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
					if ((bc >= DOF_C) && (bc < MAX_NDOFS)) {
						I = -node.m_ID[bc]-2;
						if (I>=0 && I<neq) 
							ui[I] = dq - node.m_ct[bc - DOF_C]; 
					}
			}
		}
	}

	// initialize rigid bodies
	int NO = m_fem.m_Obj.size();
	for (i=0; i<NO; ++i)
	{
		FERigidBody* prb = dynamic_cast<FERigidBody*>(m_fem.m_Obj[i]);
		if (prb)
		{
			FERigidBody& RB = *prb;

			// clear reaction forces
			RB.m_Fr = RB.m_Mr = vec3d(0,0,0);

			// store previous state
			RB.m_rp = RB.m_rt;
			RB.m_qp = RB.m_qt;
			RB.m_Up[0] = RB.m_Ut[0];
			RB.m_Up[1] = RB.m_Ut[1];
			RB.m_Up[2] = RB.m_Ut[2];
			RB.m_Up[3] = RB.m_Ut[3];
			RB.m_Up[4] = RB.m_Ut[4];
			RB.m_Up[5] = RB.m_Ut[5];

			RB.m_du[0] = RB.m_dul[0] = 0.0;
			RB.m_du[1] = RB.m_dul[1] = 0.0;
			RB.m_du[2] = RB.m_dul[2] = 0.0;
			RB.m_du[3] = RB.m_dul[3] = 0.0;
			RB.m_du[4] = RB.m_dul[4] = 0.0;
			RB.m_du[5] = RB.m_dul[5] = 0.0;
		}
	}

	// calculate local rigid displacements
	for (i=0; i<(int) m_fem.m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *m_fem.m_RDC[i];
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[DC.id]);
		if (RB.IsActive() && DC.IsActive())
		{
			int I = DC.bc;
			int lc = DC.lc;
			if (lc >= 0)
			{
				RB.m_dul[I] = DC.sf*m_fem.GetLoadCurve(lc)->Value() - RB.m_Ut[DC.bc];
			}
		}
	}

	// calculate global rigid displacements
	for (i=0; i<NO; ++i)
	{
		FERigidBody* prb = dynamic_cast<FERigidBody*>(m_fem.m_Obj[i]);
		if (prb)
		{
			FERigidBody& RB = *prb;
			if (RB.m_prb == 0)
			{
				for (j=0; j<6; ++j) RB.m_du[j] = RB.m_dul[j];
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
	for (i=0; i<NO; ++i)
	{
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[i]);
		for (j=0; j<6; ++j)
		{
			int I = -RB.m_LM[j]-2;
			if (I >= 0) ui[I] = RB.m_du[j];
		}
	}

	// apply prescribed rigid body forces
	// TODO: I don't think this does anything since
	//       the reaction forces are zeroed in the FESolidSolver::Residual function
	for (i=0; i<(int) m_fem.m_RFC.size(); ++i)
	{
		FERigidBodyForce& FC = *m_fem.m_RFC[i];
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[FC.id]);
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
	if (m_fem.ContactInterfaces() > 0) UpdateContact();

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	FEMaterialPoint::dt = m_fem.GetCurrentStep()->m_dt;
	FEMaterialPoint::time = m_fem.m_ftime;

	FEMesh& mesh = m_fem.m_mesh;
	for (i=0; i<mesh.Domains(); ++i) mesh.Domain(i).InitElements();

	// intialize the stresses
	// TODO: is this a good place to update the stresses?
	// Perhaps I should place this back in the residual routine?
	UpdateStresses();
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
//! The details of this implementation of the BFGS method can be found in:
//!   "Finite Element Procedures", K.J. Bathe, p759 and following
bool FESolidSolver::Quasin(double time)
{
	int i, n;

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

	// check for CTRL+C interruption before we do any work
	m_fem.CheckInterruption();

	// calculate initial stiffness matrix
	if (m_solvertype == 0)
	{
		if (ReformStiffness() == false) return false;
	}

	// calculate initial residual
	if (Residual(m_bfgs.m_R0) == false) return false;

	m_bfgs.m_R0 += m_Fd;

	// TODO: I can check here if the residual is zero.
	// If it is than there is probably no force acting on the system
	// if (m_R0*m_R0 < eps) bconv = true;

//	double r0 = m_R0*m_R0;

	clog.printf("\n===== beginning time step %d : %lg =====\n", pstep->m_ntimesteps+1, m_fem.m_ftime);

	// set the initial step length estimates to 1.0
	double s, olds, oldolds;  // line search step lengths from the current iteration and the two previous ones
	s=1; olds=1; oldolds=1;

	// loop until converged or when max nr of reformations reached
	do
	{
		Logfile::MODE oldmode = clog.GetMode();
		if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(pstep->GetPrintLevel() != FE_PRINT_NEVER)) clog.SetMode(Logfile::FILE_ONLY);

		clog.printf(" %d\n", m_niter+1);
		clog.SetMode(oldmode);

		// assume we'll converge. 
		bconv = true;
		if (m_solvertype==1) 	//  calculate new direction using Hager & Zhang conjugate gradient method
		{
			if ((m_niter>0)&&(breform==0))  // no need to restart CG
			{ 
				// calculate Hager- Zhang direction
        		double moddU=sqrt(u0*u0);  // needed later for the step length calculation

				// calculate yk
				vector<double> RR(m_neq);
				RR=m_bfgs.m_R1-Rold;

				// calculate dk.yk
				double bdiv=u0*RR;
				double betapcg;
				if (bdiv==0.0) 
				{
					betapcg=0.0;
					sdflag=true;
				}
    			else {
					double RR2=RR*RR;	// yk^2
               			// use m_ui as a temporary vector
					for (i=0; i<m_neq; ++i) {
						m_bfgs.m_ui[i] = RR[i]-2.0*u0[i]*RR2/bdiv;	// yk-2*dk*yk^2/(dk.yk)
						}
					betapcg=m_bfgs.m_ui*m_bfgs.m_R1;	// m_ui*gk+1
					betapcg=-betapcg/bdiv;   
          			double modR=sqrt(m_bfgs.m_R0*m_bfgs.m_R0);
          			double etak=-1.0/(moddU*min(0.01,modR));
          			betapcg=max(etak,betapcg);
					// try Fletcher - Reeves instead
					// betapcg=(m_R0*m_R0)/(m_Rold*m_Rold);
					// betapcg=0.0;
					sdflag=false;
				}

				for (i=0; i<m_neq; ++i) 
				{
            		m_bfgs.m_ui[i]=m_bfgs.m_R1[i]+betapcg*u0[i];
				}
			}
			else 
			{
				// use steepest descent for first iteration or when a restart is needed
            	m_bfgs.m_ui=m_bfgs.m_R0;
				breform=false;
				sdflag=true;
        	 }
			Rold=m_bfgs.m_R1;		// store residual for use next time
			u0=m_bfgs.m_ui;		// store direction for use on the next iteration
		}
		else if (m_solvertype==0)	// we are using the BFGS solver and need to solve the equations
		{
			// solve the equations
			m_SolverTime.start();
			{
				m_bfgs.SolveEquations(m_bfgs.m_ui, m_bfgs.m_R0);
			}
			m_SolverTime.stop();
		}
		else if (m_solvertype==2)	// we are doing an explicit analysis and need to calculate the acceleration and velocity
		{
			// get the mesh
			FEMesh& mesh = m_fem.m_mesh;
			int N = mesh.Nodes();
			double dt=m_fem.GetCurrentStep()->m_dt;
			for (i=0; i<N; ++i)
			{
				FENode& node = mesh.Node(i);
				//  calculate acceleration using F=ma and update
				if ((n = node.m_ID[DOF_X]) >= 0) node.m_at.x = m_inv_mass[n]*m_bfgs.m_R1[n];
				if ((n = node.m_ID[DOF_Y]) >= 0) node.m_at.y = m_inv_mass[n]*m_bfgs.m_R1[n];
				if ((n = node.m_ID[DOF_Z]) >= 0) node.m_at.z = m_inv_mass[n]*m_bfgs.m_R1[n];
				// and update the velocities using the accelerations
				node.m_vt = node.m_vp + node.m_at*dt;	//  update velocity using acceleration m_at
				node.m_vt.x*=m_dyn_damping;
				node.m_vt.y*=m_dyn_damping;
				node.m_vt.z*=m_dyn_damping;
				//	calculate incremental displacement using velocity
				double vel=node.m_vt.x;
				if ((n = node.m_ID[DOF_X]) >= 0) m_bfgs.m_ui[n] = node.m_vt.x*dt;
				if ((n = node.m_ID[DOF_Y]) >= 0) m_bfgs.m_ui[n] = node.m_vt.y*dt;
				if ((n = node.m_ID[DOF_Z]) >= 0) m_bfgs.m_ui[n] = node.m_vt.z*dt;
			}
		}

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
		if (m_solvertype==0) // we are using the BFGS solver
		{
			if (m_bfgs.m_LStol > 0) s = m_bfgs.LineSearch(1.0);
			else
			{
				s = 1;

				// Update geometry
				Update(m_bfgs.m_ui);

				// calculate residual at this point
				Residual(m_bfgs.m_R1);
			}
		}
		else if (m_solvertype==1)	// we are using the Hager - Zhang solver, which starts with a guess for the step length 
		{
			// use the step length from two steps previously as the initial guess
			// note that it has its own linesearch, different from the BFGS one
			s = m_bfgs.LineSearchCG(oldolds);
			// update the old step lengths for use as an initial guess in two iterations' time
			if (m_niter<1) oldolds=s;	// if this is the first iteration, use current step length
			else oldolds=olds;	// otherwise use the previous one
			olds=s;  // and store the current step to be used for the iteration after next
		}
		else if (m_solvertype==2) // need to update everything for the explicit solver
		{
			// Update geometry
			Update(m_bfgs.m_ui);

			// calculate residual at this point
			Residual(m_bfgs.m_R1);
		}

		// update total displacements
		int neq = m_Ui.size();
		if (m_solvertype==2) s=1.0;  // we have not done a line search in the explicit solver 
		for (i=0; i<neq; ++i) m_Ui[i] += s*m_bfgs.m_ui[i];

		// calculate norms
		if (m_solvertype<2) 
		{
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
			oldmode = clog.GetMode();
			if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
				(pstep->GetPrintLevel() != FE_PRINT_NEVER)) clog.SetMode(Logfile::FILE_ONLY);

			clog.printf(" Nonlinear solution status: time= %lg\n", time); 
			clog.printf("\tstiffness updates             = %d\n", m_bfgs.m_nups);
			clog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
			clog.printf("\tstiffness matrix reformations = %d\n", m_nref);
			if (m_bfgs.m_LStol > 0) clog.printf("\tstep from line search         = %lf\n", s);
			clog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
			clog.printf("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
			clog.printf("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
			clog.printf("\t   displacement     %15le %15le %15le \n", normUi, normu ,(m_Dtol*m_Dtol)*normU );

			clog.SetMode(oldmode);
		}

		// check if we have converged. 
		// If not, calculate the BFGS update vectors
		if (bconv == false)
		{
			if ((normR1 < m_Rmin))
			{
				// check for almost zero-residual on the first iteration
				// this might be an indication that there is no force on the system
				clog.printbox("WARNING", "No force acting on the system.");
				bconv = true;
			}
			else if (s < m_bfgs.m_LSmin)
			{
				// check for zero linestep size
				clog.printbox("WARNING", "Zero linestep size. Stiffness matrix will now be reformed");
				breform = true;
			}
			else if (normE1 > normEm)
			{
				// check for diverging
				if (m_solvertype==0) clog.printbox("WARNING", "Problem is diverging. Stiffness matrix will now be reformed");
				normEm = normE1;
				normEi = normE1;
				normRi = normR1;
				breform = true;
			}
			else
			{
				// If we havn't reached max nr of BFGS updates
				// do an update
				if (!breform && (m_solvertype==0))
				{
					if (m_bfgs.m_nups < m_bfgs.m_maxups-1)
					{
						if (m_bfgs.Update(s, m_bfgs.m_ui, m_bfgs.m_R0, m_bfgs.m_R1) == false)
						{
							// Stiffness update has failed.
							// this might be due a too large condition number
							// or the update was no longer positive definite.
							clog.printbox("WARNING", "The BFGS update has failed.\nStiffness matrix will now be reformed.");
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
							clog.printbox("WARNING", "Max nr of iterations reached.\nStiffness matrix will now be reformed.");

					}
				}
			}	

			// zero displacement increments
			// we must set this to zero before the reformation
			// because we assume that the prescribed displacements are stored 
			// in the m_ui vector.
			zero(m_bfgs.m_ui);

			// reform stiffness matrices if necessary
			if (breform && (m_solvertype == 0))
			{
				clog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);

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

			clog.printf("\n........................ augmentation # %d\n", m_naug+1);

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
				if (pstep->m_psolver->m_bfgs.m_maxups == 0)
				{
					clog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);
					if (ReformStiffness() == false) break;
				}
			}
		}
	
		// increase iteration number
		m_niter++;

		// let's flush the logfile to make sure the last output will not get lost
		clog.flush();

		// check for CTRL+C interruption
		m_fem.CheckInterruption();
	}
	while (bconv == false);

	// when converged, 
	// print a convergence summary to the clog file
	if (bconv)
	{
		Logfile::MODE mode = clog.SetMode(Logfile::FILE_ONLY);
		if (mode != Logfile::NEVER)
		{
			clog.printf("\nconvergence summary\n");
			clog.printf("    number of iterations   : %d\n", m_niter);
			clog.printf("    number of reformations : %d\n", m_nref);
		}
		clog.SetMode(mode);
	}

	// if converged we update the total displacements
	if (bconv)
	{
		m_Ut += m_Ui;
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Reforms a stiffness matrix and factorizes it
bool FESolidSolver::ReformStiffness()
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
		m_breshape = (m_fem.ContactInterfaces() > 0? true : false);
	}

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
//! Calculates global stiffness matrix.

bool FESolidSolver::StiffnessMatrix()
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
	FEMesh& mesh = m_fem.m_mesh;

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
		int NF = m_fem.BodyForces();
		for (int j=0; j<NF; ++j)
		{
			FEBodyForce& BF = *m_fem.GetBodyForce(j);
			dom.BodyForceStiffness(this, BF);
		}
	}

	// Add inertial stiffness for dynamic problems
	if (m_fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC)
	{
		for (i=0; i<mesh.Domains(); ++i) 
		{
			FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
			dom.InertialStiffness(this);
		}
	}

	// calculate contact stiffness
	if (m_fem.ContactInterfaces() > 0) 
	{
		ContactStiffness();
	}

	// calculate stiffness matrices for surface loads
	int nsl = (int) m_fem.m_SL.size();
	for (i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = m_fem.m_SL[i];

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
	int NRB = m_fem.m_Obj.size();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[i]);
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
			if (K.diag(i) == 0) zd.push_back(i);
		}

//		if (zd.empty() == false) throw ZeroDiagonal(zd, m_fem);
		if (zd.empty() == false) throw ZeroDiagonal(-1, -1);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FESolidSolver::NonLinearConstraintStiffness()
{
	int N = m_fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		plc->StiffnessMatrix(this);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact stiffness matrix

void FESolidSolver::ContactStiffness()
{
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.ContactInterface(i)->ContactStiffness(this);
}

//-----------------------------------------------------------------------------
//! This function calculates the rigid stiffness matrices

void FESolidSolver::RigidStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{
	int i, j, k, l, n = en.size();
	double kij[MAX_NDOFS][MAX_NDOFS], Ri[3][3] = {0}, Rj[3][3] = {0};

	double KF[MAX_NDOFS][6];
	double KR[6][6];

	int *lmi, *lmj;
	int I, J;

	SparseMatrix& K = *m_pK;

	vec3d ai, aj;

	int ndof = ke.columns() / n;

	vector<double>& ui = m_bfgs.m_ui;
	FEMesh& mesh = m_fem.m_mesh;

	// loop over columns
	for (j=0; j<n; ++j)
	{
		FENode& nodej = mesh.Node(en[j]);
		if (nodej.m_rid >= 0)
		{
			// this is a rigid interface node
			// get the rigid body this node is attached to
			FERigidBody& RBj = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[nodej.m_rid]);

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
					FERigidBody& RBi = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[nodei.m_rid]);

					lmi = RBi.m_LM;
					
					// get the relative distance
					ai = nodei.m_rt - RBi.m_rt;
	
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
								if (J < -1) m_Fd[I] -= KR[l][k]*ui[-J-2];
								else if (J >= 0) K.add(I,J, KR[l][k]);
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
								if (J < -1) m_Fd[I] -= KF[l][k]*ui[-J-2];
								else if (J >= 0) K.add(I,J, KF[l][k]);
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
					FERigidBody& RBi = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[nodei.m_rid]);

					// get the rigid body equation nrs.
					lmi = RBi.m_LM;

					// get the relative distance to the center of mass
					ai = nodei.m_rt - RBi.m_rt;

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
//!  Assembles the element stiffness matrix into the global stiffness matrix.
//!  Also adjusts the global stiffness matrix and residual to take the 
//!  prescribed displacements into account.

// TODO: In stead of changing the global stiffness matrix to accomodate for 
// the rigid bodies and linear constraints, can I modify the element stiffness
// matrix prior to assembly? I might have to change the elm vector as well as 
// the element matrix size.

void FESolidSolver::AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{
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
	if (m_fem.m_Obj.empty() == false) RigidStiffness(en, elm, ke);
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FESolidSolver::ContactForces(vector<double>& R)
{
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.ContactInterface(i)->ContactForces(R, this);
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FESolidSolver::Residual(vector<double>& R)
{
	int i;
	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// zero rigid body reaction forces
	int NRB = m_fem.m_Obj.size();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[i]);
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);
	}

	// get the mesh
	FEMesh& mesh = m_fem.m_mesh;

	// calculate the internal (stress) forces
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(this, R);
	}

	// update body forces
	for (i=0; i<m_fem.BodyForces(); ++i)
	{
		// TODO: I don't like this but for now I'll hard-code the modification of the
		//       force center position
		FEPointBodyForce* pbf = dynamic_cast<FEPointBodyForce*>(m_fem.GetBodyForce(i));
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
		for (int j=0; j<m_fem.BodyForces(); ++j)
		{
			FEBodyForce& BF = *m_fem.GetBodyForce(j);
			dom.BodyForce(this, BF, R);
		}
	}

	// calculate inertial forces for dynamic problems
	if (m_fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC) InertialForces(R);

	// calculate forces due to surface loads
	int nsl = (int) m_fem.m_SL.size();
	for (i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = m_fem.m_SL[i];
		if (psl->IsActive()) psl->Residual(this, R);
	}

	// calculate contact forces
	if (m_fem.ContactInterfaces() > 0)
	{
		ContactForces(R);
	}

	// calculate nonlinear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	NonLinearConstraintForces(R);

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
void FESolidSolver::NonLinearConstraintForces(vector<double> &R)
{
	int N = m_fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		plc->Residual(this, R);
	}
}

//-----------------------------------------------------------------------------
//!  Assembles the element into the global residual. This function
//!  also checks for rigid dofs and assembles the residual using a condensing
//!  procedure in the case of rigid dofs.

void FESolidSolver::AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R)
{
	int i, j, I, n, l;
	vec3d a, d;

	// assemble the element residual into the global residual
	int ndof = fe.size();
	for (i=0; i<ndof; ++i)
	{
		I = elm[i];
		if ( I >= 0) R[I] += fe[i];
		else if (-I-2 >= 0) m_Fr[-I-2] -= fe[i];
	}

	int ndn = ndof / en.size();

	// if there are linear constraints we need to apply them
	if (m_fem.m_LinC.size() > 0)
	{
		// loop over all degrees of freedom of this element
		for (i=0; i<ndof; ++i)
		{
			// see if this dof belongs to a linear constraint
			n = MAX_NDOFS*(en[i/ndn]) + i%ndn;
			l = m_fem.m_LCT[n];
			if (l >= 0)
			{
				// if so, get the linear constraint
				FELinearConstraint& lc = *m_fem.m_LCA[l];
				assert(elm[i] == -1);
	
				// now loop over all "slave" nodes and
				// add the contribution to the residual
				int ns = lc.slave.size();
				list<FELinearConstraint::SlaveDOF>::iterator is = lc.slave.begin();
				for (j=0; j<ns; ++j, ++is)
				{
					I = is->neq;
					if (I >= 0)
					{
						double A = is->val;
						R[I] += A*fe[i];
					}
				}
			}
		}
	}

	// If there are rigid bodies we need to look for rigid dofs
	if (m_fem.m_Obj.empty() == false)
	{
		int *lm;

		for (i=0; i<ndof; i+=ndn)
		{
			FENode& node = m_fem.m_mesh.Node(en[i/ndn]);
			if (node.m_rid >= 0)
			{
				vec3d F(fe[i], fe[i+1], fe[i+2]);

				// this is an interface dof
				// get the rigid body this node is connected to
				FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[node.m_rid]);
				lm = RB.m_LM;

				// add to total torque of this body
				a = node.m_rt - RB.m_rt;

				n = lm[3]; if (n >= 0) R[n] += a.y*F.z-a.z*F.y; RB.m_Mr.x -= a.y*F.z-a.z*F.y;
				n = lm[4]; if (n >= 0) R[n] += a.z*F.x-a.x*F.z; RB.m_Mr.y -= a.z*F.x-a.x*F.z;
				n = lm[5]; if (n >= 0) R[n] += a.x*F.y-a.y*F.x; RB.m_Mr.z -= a.x*F.y-a.y*F.x;
/*
				// if the rotational degrees of freedom are constrained for a rigid node
				// then we need to add an additional component to the residual
				if (node.m_ID[DOF_RU] == lm[3])
				{
					d = node.m_Dt;
					n = lm[3]; if (n >= 0) R[n] += d.y*F.z-d.z*F.y; RB.m_Mr.x -= d.y*F.z-d.z*F.y;
					n = lm[4]; if (n >= 0) R[n] += d.z*F.x-d.x*F.z; RB.m_Mr.y -= d.z*F.x-d.x*F.z;
					n = lm[5]; if (n >= 0) R[n] += d.x*F.y-d.y*F.x; RB.m_Mr.z -= d.x*F.y-d.y*F.x;
				}
*/
				// add to global force vector
				n = lm[0]; if (n >= 0) R[n] += F.x; RB.m_Fr.x -= F.x;
				n = lm[1]; if (n >= 0) R[n] += F.y; RB.m_Fr.y -= F.y;
				n = lm[2]; if (n >= 0) R[n] += F.z; RB.m_Fr.z -= F.z;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the concentrated nodal forces

void FESolidSolver::NodalForces(vector<double>& F)
{
	int i, id, bc, lc, n;
	double s, f;
	vec3d a;
	int* lm;

	// zero nodal force vector
	zero(F);

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
		
			f = s*m_fem.GetLoadCurve(lc)->Value();
			
			// For pressure and concentration loads, multiply by dt
			// for consistency with evaluation of residual and stiffness matrix
			if ((bc == DOF_P) || (bc >= DOF_C))
				f *= m_fem.GetCurrentStep()->m_dt;

			if (n >= 0) F[n] = f;
			else if (node.m_rid >=0)
			{
				// this is a rigid body node
				FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.m_Obj[node.m_rid]);

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

void FESolidSolver::InertialForces(vector<double>& R)
{
	// get the mesh
	FEMesh& mesh = m_fem.m_mesh;

	// allocate F
	vector<double> F(3*mesh.Nodes());
	zero(F);

	// calculate F
	double dt = m_fem.GetCurrentStep()->m_dt;
	double a = 4.0 / dt;
	double b = a / dt;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vec3d& rt = node.m_rt;
		vec3d& rp = node.m_rp;
		vec3d& vp = node.m_vp;
		vec3d& ap = node.m_ap;

		F[3*i  ] = b*(rt.x - rp.x) - a*vp.x - ap.x;
		F[3*i+1] = b*(rt.y - rp.y) - a*vp.y - ap.y;
		F[3*i+2] = b*(rt.z - rp.z) - a*vp.z - ap.z;
	}

	// now multiply F with the mass matrix
	matrix ke;
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(nd));
		dom.InertialForces(this, R, F);
	}
}
