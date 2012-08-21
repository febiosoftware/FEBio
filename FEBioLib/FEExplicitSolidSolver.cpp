#include "stdafx.h"
#include "FEExplicitSolidSolver.h"
#include "FEElasticSolidDomain.h"
#include "FERigidBody.h"
#include "FERigid.h"
#include "FEPointBodyForce.h"
#include "log.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEExplicitSolidSolver, FESolver)
	ADD_PARAMETER(m_dyn_damping, FE_PARAM_DOUBLE, "dyn_damping");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEExplicitSolidSolver::FEExplicitSolidSolver(FEModel& fem) : FESolver(fem)
{
	m_dyn_damping = 0.99;

	m_niter = 0;
	m_nreq = 0;
}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver::Init()
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
	m_inv_mass.assign(neq, 1);

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
	matrix ke;
	int j, iel;
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

	return true;
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
bool FEExplicitSolidSolver::InitEquations()
{
	int i, j, n;

	// get the mesh
	FEMesh& mesh = m_fem.m_mesh;

	// initialize nr of equations
	int neq = 0;

	// give all free dofs an equation number
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (j=0; j<MAX_NDOFS; ++j)
			if (node.m_ID[j] >= 0) node.m_ID[j] = neq++;
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
//!  Assembles the element into the global residual. This function
//!  also checks for rigid dofs and assembles the residual using a condensing
//!  procedure in the case of rigid dofs.

void FEExplicitSolidSolver::AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R)
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
//! Updates the current state of the model
void FEExplicitSolidSolver::Update(vector<double>& ui)
{
	// update kinematics
	UpdateKinematics(ui);

	// update contact
	// Update all contact interfaces
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.ContactInterface(i)->Update(m_niter);

	// update element stresses
	UpdateStresses();

	// update other stuff that may depend on the deformation
	int NBF = m_fem.BodyForces();
	for (int i=0; i<NBF; ++i) m_fem.GetBodyForce(i)->Update();

	// dump all states to the plot file when requested
	if (m_fem.GetCurrentStep()->m_nplot == FE_PLOT_MINOR_ITRS) m_fem.Write();
}


//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEExplicitSolidSolver::UpdateKinematics(vector<double>& ui)
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
}

//-----------------------------------------------------------------------------
//! Updates the rigid body data
void FEExplicitSolidSolver::UpdateRigidBodies(vector<double>& ui)
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
void FEExplicitSolidSolver::UpdateStresses()
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
//! Save data to dump file

void FEExplicitSolidSolver::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nrhs;
		ar << m_niter;
		ar << m_nref << m_ntotref;
		ar << m_naug;
		ar << m_neq << m_nreq;

		ar << m_bfgs.m_maxups;
		ar << m_bfgs.m_maxref;
		ar << m_bfgs.m_cmax;
		ar << m_bfgs.m_nups;
	}
	else
	{
		ar >> m_nrhs;
		ar >> m_niter;
		ar >> m_nref >> m_ntotref;
		ar >> m_naug;
		ar >> m_neq >> m_nreq;

		ar >> m_bfgs.m_maxups;
		ar >> m_bfgs.m_maxref;
		ar >> m_bfgs.m_cmax;
		ar >> m_bfgs.m_nups;
	}
}

//-----------------------------------------------------------------------------
//!  This function mainly calls the Quasin routine 
//!  and deals with exceptions that require the immediate termination of
//!	quasi-Newton iterations.
bool FEExplicitSolidSolver::SolveStep(double time)
{
	bool bret;

	try
	{
		// let's try to solve the step
		bret = DoSolve(time);
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
void FEExplicitSolidSolver::PrepStep(double time)
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
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.ContactInterface(i)->Update(m_niter);

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
//! calculates the concentrated nodal forces

void FEExplicitSolidSolver::NodalForces(vector<double>& F)
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
bool FEExplicitSolidSolver::DoSolve(double time)
{
	int i, n;

	vector<double> u0(m_neq);
	vector<double> Rold(m_neq);

	// Get the current step
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// prepare for the first iteration
	PrepStep(time);

	// check for CTRL+C interruption before we do any work
	m_fem.CheckInterruption();

	// calculate initial residual
	if (Residual(m_bfgs.m_R0) == false) return false;

	m_bfgs.m_R0 += m_Fd;

	clog.printf("\n===== beginning time step %d : %lg =====\n", pstep->m_ntimesteps+1, m_fem.m_ftime);

	Logfile::MODE oldmode = clog.GetMode();
	if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
		(pstep->GetPrintLevel() != FE_PRINT_NEVER)) clog.SetMode(Logfile::FILE_ONLY);

	clog.printf(" %d\n", m_niter+1);
	clog.SetMode(oldmode);

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

	// need to update everything for the explicit solver
	// Update geometry
	Update(m_bfgs.m_ui);

	// calculate residual at this point
	Residual(m_bfgs.m_R1);

	// update total displacements
	int neq = m_Ui.size();
	for (i=0; i<neq; ++i) m_Ui[i] += m_bfgs.m_ui[i];

	// increase iteration number
	m_niter++;

	// let's flush the logfile to make sure the last output will not get lost
	clog.flush();

	// check for CTRL+C interruption
	m_fem.CheckInterruption();

	// when converged, 
	// print a convergence summary to the clog file
	Logfile::MODE mode = clog.SetMode(Logfile::FILE_ONLY);
	if (mode != Logfile::NEVER)
	{
		clog.printf("\nconvergence summary\n");
		clog.printf("    number of iterations   : %d\n", m_niter);
		clog.printf("    number of reformations : %d\n", m_nref);
	}
	clog.SetMode(mode);

	// if converged we update the total displacements
	m_Ut += m_Ui;

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEExplicitSolidSolver::Residual(vector<double>& R)
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
//! Calculates the contact forces
void FEExplicitSolidSolver::ContactForces(vector<double>& R)
{
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.ContactInterface(i)->ContactForces(R, this);
}


//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces 
void FEExplicitSolidSolver::NonLinearConstraintForces(vector<double> &R)
{
	int N = m_fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		plc->Residual(this, R);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the inertial forces for dynamic problems

void FEExplicitSolidSolver::InertialForces(vector<double>& R)
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
