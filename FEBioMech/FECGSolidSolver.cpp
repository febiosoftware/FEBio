#include "stdafx.h"
#include "FECGSolidSolver.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEMesh.h"
#include "FECore/log.h"
#include "FEContactInterface.h"
#include "FEUncoupledMaterial.h"
#include "FEResidualVector.h"
#include "FEElasticDomain.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FERigidSystem.h"
#include "FERigidBody.h"
#include "RigidBC.h"
#include <FECore/FEPrescribedBC.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FELinearConstraintManager.h>
#include "FEBodyForce.h"
#include "FECore/sys.h"
#include "FEMechModel.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FECGSolidSolver, FESolver)
	ADD_PARAMETER(m_Dtol  , "dtol");
	ADD_PARAMETER(m_Etol  , "etol");
	ADD_PARAMETER(m_Rtol  , "rtol");
	ADD_PARAMETER(m_Rmin  , "min_residual");
	ADD_PARAMETER(m_beta  , "beta");
	ADD_PARAMETER(m_gamma , "gamma");
	ADD_PARAMETER(m_LStol , "lstol");
	ADD_PARAMETER(m_LSmin , "lsmin");
	ADD_PARAMETER(m_LSiter, "lsiter");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FECGSolidSolver::FECGSolidSolver(FEModel* pfem) : FESolver(pfem)
{
	// default values
	m_Rtol = 0;	// deactivate residual convergence 
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_Rmin = 1.0e-20;

	m_LStol = 0;
	m_LSmin = 0;
	m_LSiter = 10;

	m_niter = 0;
	m_nreq = 0;

	// default Newmark parameters for unconditionally stable time integration
	m_beta = 0.25;
	m_gamma = 0.5;

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int varD = dofs.AddVariable("displacement", VAR_VEC3);
	dofs.SetDOFName(varD, 0, "x");
	dofs.SetDOFName(varD, 1, "y");
	dofs.SetDOFName(varD, 2, "z");
	int varQ = dofs.AddVariable("rotation", VAR_VEC3);
	dofs.SetDOFName(varQ, 0, "u");
	dofs.SetDOFName(varQ, 1, "v");
	dofs.SetDOFName(varQ, 2, "w");
	int varQR = dofs.AddVariable("rigid rotation", VAR_VEC3);
	dofs.SetDOFName(varQR, 0, "Ru");
	dofs.SetDOFName(varQR, 1, "Rv");
	dofs.SetDOFName(varQR, 2, "Rw");
	int varV = dofs.AddVariable("velocity", VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");

	// get the DOF indices
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
	m_dofVX = pfem->GetDOFIndex("vx");
	m_dofVY = pfem->GetDOFIndex("vy");
	m_dofVZ = pfem->GetDOFIndex("vz");
	m_dofU = pfem->GetDOFIndex("u");
	m_dofV = pfem->GetDOFIndex("v");
	m_dofW = pfem->GetDOFIndex("w");
	m_dofRU = pfem->GetDOFIndex("Ru");
	m_dofRV = pfem->GetDOFIndex("Rv");
	m_dofRW = pfem->GetDOFIndex("Rw");
}

//-----------------------------------------------------------------------------
void FECGSolidSolver::Clean()
{

}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FESolidSolver
//
bool FECGSolidSolver::Init()
{
	if (FESolver::Init() == false) return false;

	// check parameters
	if (m_Dtol <  0.0) { felog.printf("Error: dtol must be nonnegative.\n"   ); return false; }
	if (m_Etol <  0.0) { felog.printf("Error: etol must be nonnegative.\n"); return false; }
	if (m_Rtol <  0.0) { felog.printf("Error: rtol must be nonnegative.\n"); return false; }
	if (m_Rmin <  0.0) { felog.printf("Error: min_residual must be nonnegative.\n"  ); return false; }
	if (m_LStol  < 0.) { felog.printf("Error: lstol must be nonnegative.\n" ); return false; }
	if (m_LSmin  < 0.) { felog.printf("Error: lsmin must be nonnegative.\n" ); return false; }
	if (m_LSiter < 0 ) { felog.printf("Error: lsiter must be nonnegative.\n"  ); return false; }

	// get nr of equations
	int neq = m_neq;

	// allocate vectors
	m_Fn.assign(neq, 0);
	m_Fr.assign(neq, 0);
	m_Ui.assign(neq, 0);
	m_Ut.assign(neq, 0);

	// we need to fill the total displacement vector m_Ut
	FEMesh& mesh = m_fem.GetMesh();
	gather(m_Ut, mesh, m_dofX);
	gather(m_Ut, mesh, m_dofY);
	gather(m_Ut, mesh, m_dofZ);
	gather(m_Ut, mesh, m_dofU);
	gather(m_Ut, mesh, m_dofV);
	gather(m_Ut, mesh, m_dofW);

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
bool FECGSolidSolver::InitEquations()
{
	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// initialize nr of equations
	int neq = 0;

	// give all free dofs an equation number
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j<(int)node.m_ID.size(); ++j)
		{
			if (node.m_ID[j] == DOF_FIXED) { node.m_ID[j] = -1; }
			else if (node.m_ID[j] == DOF_OPEN) { node.m_ID[j] = neq++; }
			else if (node.m_ID[j] == DOF_PRESCRIBED) { node.m_ID[j] = -neq - 2; neq++; }
			else { assert(false); return false; }
		}
	}

	// Next, we assign equation numbers to the rigid body degrees of freedom
	m_nreq = neq;
	FEMechModel& fem = static_cast<FEMechModel&>(m_fem);
	FERigidSystem& rigid = *fem.GetRigidSystem();
	int nrb = rigid.Objects();
	for (int i = 0; i<nrb; ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		for (int j = 0; j<6; ++j)
		{
			int bcj = RB.m_BC[j];
			int lmj = RB.m_LM[j];
			if (bcj == DOF_OPEN) { RB.m_LM[j] = neq; neq++; }
			else if (bcj == DOF_PRESCRIBED) { RB.m_LM[j] = -neq - 2; neq++; }
			else if (bcj == DOF_FIXED) { RB.m_LM[j] = -1; }
			else { assert(false); return false; }
		}
	}

	// store the number of equations
	m_neq = neq;

	// we assign the rigid body equation number to
	// Also make sure that the nodes are NOT constrained!
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0)
		{
			FERigidBody& RB = *rigid.Object(node.m_rid);
			node.m_ID[m_dofX ] = (RB.m_LM[0] >= 0 ? -RB.m_LM[0] - 2 : RB.m_LM[0]);
			node.m_ID[m_dofY ] = (RB.m_LM[1] >= 0 ? -RB.m_LM[1] - 2 : RB.m_LM[1]);
			node.m_ID[m_dofZ ] = (RB.m_LM[2] >= 0 ? -RB.m_LM[2] - 2 : RB.m_LM[2]);
			node.m_ID[m_dofRU] = (RB.m_LM[3] >= 0 ? -RB.m_LM[3] - 2 : RB.m_LM[3]);
			node.m_ID[m_dofRV] = (RB.m_LM[4] >= 0 ? -RB.m_LM[4] - 2 : RB.m_LM[4]);
			node.m_ID[m_dofRW] = (RB.m_LM[5] >= 0 ? -RB.m_LM[5] - 2 : RB.m_LM[5]);
		}
	}

	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration. 
void FECGSolidSolver::PrepStep()
{
	TRACK_TIME("update");

	const FETimeInfo& tp = m_fem.GetTime();

	// initialize counters
	m_niter = 0;	// nr of iterations
	m_nrhs = 0;	// nr of RHS evaluations
	m_nref = 0;	// nr of stiffness reformations
	m_ntotref = 0;
	m_naug = 0;	// nr of augmentations

	// zero total displacements
	zero(m_Ui);

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEMesh& mesh = m_fem.GetMesh();
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		ni.m_rp = ni.m_rt;
		ni.m_vp = ni.get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
		ni.m_ap = ni.m_at;
	}

	// apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
	NodalForces(m_Fn, tp);

	// apply prescribed displacements
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_ui;
	zero(ui);
	int neq = m_neq;
	int nbc = m_fem.PrescribedBCs();
	for (int i = 0; i<nbc; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive()) dc.PrepStep(ui);
	}

	// initialize rigid bodies
	FEMechModel& fem = static_cast<FEMechModel&>(m_fem);
	FERigidSystem& rigid = *fem.GetRigidSystem();
	int NO = rigid.Objects();
	for (int i = 0; i<NO; ++i) rigid.Object(i)->Init();

	// calculate local rigid displacements
	for (int i = 0; i<(int)rigid.PrescribedBCs(); ++i)
	{
		FERigidBodyDisplacement& DC = *rigid.PrescribedBC(i);
		FERigidBody& RB = *rigid.Object(DC.id);
		if (DC.IsActive())
		{
			int I = DC.bc;
			int lc = DC.lc;
			if (lc >= 0)
			{
				RB.m_dul[I] = DC.Value() - RB.m_Ut[DC.bc];
			}
		}
	}

	// calculate global rigid displacements
	for (int i = 0; i<NO; ++i)
	{
		FERigidBody* prb = rigid.Object(i);
		if (prb)
		{
			FERigidBody& RB = *prb;
			if (RB.m_prb == 0)
			{
				for (int j = 0; j<6; ++j) RB.m_du[j] = RB.m_dul[j];
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
				quatd Q0 = RB.GetRotation();

				dr = Q0*dr;
				dq = Q0*dq*Q0.Inverse();

				while (pprb)
				{
					vec3d r1 = pprb->m_rt;
					dul = pprb->m_dul;

					quatd Q1 = pprb->GetRotation();

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
	for (int i = 0; i<NO; ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		for (int j = 0; j<6; ++j)
		{
			int I = -RB.m_LM[j] - 2;
			if (I >= 0) ui[I] = RB.m_du[j];
		}
	}

	FEAnalysis* pstep = m_fem.GetCurrentStep();
	if (pstep->m_nanalysis == FE_DYNAMIC)
	{
		FEMesh& mesh = m_fem.GetMesh();
		FEMechModel& fem = static_cast<FEMechModel&>(m_fem);
		FERigidSystem& rigid = *fem.GetRigidSystem();

		// set the initial velocities of all rigid nodes
		for (int i = 0; i<mesh.Nodes(); ++i)
		{
			FENode& n = mesh.Node(i);
			if (n.m_rid >= 0)
			{
				FERigidBody& rb = *rigid.Object(n.m_rid);
				vec3d V = rb.m_vt;
				vec3d W = rb.m_wt;
				vec3d r = n.m_rt - rb.m_rt;

				vec3d v = V + (W ^ r);
				n.m_vp = v;
				n.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, v);

				vec3d a = (W ^ V)*2.0 + (W ^ (W ^ r));
				n.m_ap = n.m_at = a;
			}
		}
	}

	// store the current rigid body reaction forces
	for (int i = 0; i<rigid.Objects(); ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		RB.m_Fp = RB.m_Fr;
		RB.m_Mp = RB.m_Mr;
	}

	// initialize contact
	if (m_fem.SurfacePairConstraints() > 0) UpdateContact();

	// initialize nonlinear constraints
	if (m_fem.NonlinearConstraints() > 0) UpdateConstraints();

	// intialize material point data
	for (int i = 0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);

	// update stresses
	UpdateModel();

	// see if we need to do contact augmentations
	m_baugment = false;
	for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface& ci = dynamic_cast<FEContactInterface&>(*m_fem.SurfacePairConstraint(i));
		if (ci.IsActive() && ci.m_blaugon) m_baugment = true;
	}

	// see if we need to do incompressible augmentations
	int nmat = m_fem.Materials();
	for (int i = 0; i<nmat; ++i)
	{
		FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_fem.GetMaterial(i));
		if (pmi && pmi->m_blaugon) m_baugment = true;
	}

	// see if we have to do nonlinear constraint augmentations
	for (int i = 0; i<m_fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint& ci = *m_fem.NonlinearConstraint(i);
		if (ci.IsActive()) m_baugment = true;
	}
}

//-----------------------------------------------------------------------------
bool FECGSolidSolver::SolveStep()
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
	bool breform = false;	// reformation flag
	bool sdflag = true;		// flag for steepest descent iterations in NLCG

	// Get the current step
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// prepare for the first iteration
	const FETimeInfo& tp = m_fem.GetTime();
	PrepStep();

	// update stresses
	UpdateModel();

	// calculate initial residual
	if (Residual(m_R0) == false) return false;

	m_R0 += m_Fd;

	// TODO: I can check here if the residual is zero.
	// If it is than there is probably no force acting on the system
	// if (m_R0*m_R0 < eps) bconv = true;

//	double r0 = m_R0*m_R0;

	// set the initial step length estimates to 1.0
	double s, olds, oldolds;  // line search step lengths from the current iteration and the two previous ones
	s=1; olds=1; oldolds=1;

	// loop until converged or when max nr of reformations reached
	bool bconv = false;		// convergence flag
	do
	{
		Logfile::MODE oldmode = felog.GetMode();
		if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::LOG_FILE);

		felog.printf(" %d\n", m_niter+1);
		felog.SetMode(oldmode);

		// assume we'll converge. 
		bconv = true;
		if ((m_niter>0)&&(breform==0))  // no need to restart CG
		{ 
			// calculate Hager- Zhang direction
        	double moddU=sqrt(u0*u0);  // needed later for the step length calculation

			// calculate yk
			vector<double> RR(m_neq);
			RR=m_R1-Rold;

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
					m_ui[i] = RR[i]-2.0*u0[i]*RR2/bdiv;	// yk-2*dk*yk^2/(dk.yk)
					}
				betapcg=m_ui*m_R1;	// m_ui*gk+1
				betapcg=-betapcg/bdiv;   
          		double modR=sqrt(m_R0*m_R0);
          		double etak=-1.0/(moddU*min(0.01,modR));
          		betapcg=max(etak,betapcg);
				// try Fletcher - Reeves instead
				// betapcg=(m_R0*m_R0)/(m_Rold*m_Rold);
				// betapcg=0.0;
				sdflag=false;
			}

			for (i=0; i<m_neq; ++i) 
			{
            	m_ui[i]=m_R1[i]+betapcg*u0[i];
			}
		}
		else 
		{
			// use steepest descent for first iteration or when a restart is needed
            m_ui=m_R0;
			breform=false;
			sdflag=true;
        	}
		Rold=m_R1;		// store residual for use next time
		u0=m_ui;		// store direction for use on the next iteration

		// check for nans
		double du = m_ui*m_ui;
		if (ISNAN(du)) throw NANDetected();

		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_R0*m_R0);
			normEi = fabs(m_ui*m_R0);
			normUi = fabs(m_ui*m_ui);
			normEm = normEi;
		}

		// perform a linesearch
		// the geometry is also updated in the line search
		// use the step length from two steps previously as the initial guess
		// note that it has its own linesearch, different from the BFGS one
		s = LineSearchCG(oldolds);
		// update the old step lengths for use as an initial guess in two iterations' time
		if (m_niter<1) oldolds=s;	// if this is the first iteration, use current step length
		else oldolds=olds;	// otherwise use the previous one
		olds=s;  // and store the current step to be used for the iteration after next

		// update total displacements
		int neq = (int)m_Ui.size();
		for (i=0; i<neq; ++i) m_Ui[i] += s*m_ui[i];

		// calculate norms
		normR1 = m_R1*m_R1;
		normu  = (m_ui*m_ui)*(s*s);
		normU  = m_Ui*m_Ui;
		normE1 = s*fabs(m_ui*m_R1);

		// check residual norm
		if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;	

		// check displacement norm
		if ((m_Dtol > 0) && (normu  > (m_Dtol*m_Dtol)*normU )) bconv = false;

		// check energy norm
		if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;

		// check linestep size
		if ((m_LStol > 0) && (s < m_LSmin)) bconv = false;

		// check energy divergence
		if (normE1 > normEm) bconv = false;

		// print convergence summary
		oldmode = felog.GetMode();
		if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::LOG_FILE);

		felog.printf(" Nonlinear solution status: time= %lg\n", tp.currentTime);
		felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
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
			if (s < m_LSmin)
			{
				// check for zero linestep size
				felog.printbox("WARNING", "Zero linestep size. Stiffness matrix will now be reformed");
				breform = true;
			}
			else if (normE1 > normEm)
			{
				// check for diverging
				normEm = normE1;
				normEi = normE1;
				normRi = normR1;
				breform = true;
			}

			// zero displacement increments
			// we must set this to zero before the reformation
			// because we assume that the prescribed displacements are stored 
			// in the m_ui vector.
			zero(m_ui);

			// copy last calculated residual
			m_R0 = m_R1;
		}
		else if (m_baugment)
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
				UpdateModel();
				Residual(m_R0);
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

	// if converged we update the total displacements
	if (bconv)
	{
		m_Ut += m_Ui;
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FECGSolidSolver::UpdateKinematics(vector<double>& ui)
{
	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// update rigid bodies
	UpdateRigidBodies(ui);

	// total displacements
	vector<double> U(m_Ut.size());
	for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];

	// update flexible nodes
	// translational dofs
	scatter(U, mesh, m_dofX);
	scatter(U, mesh, m_dofY);
	scatter(U, mesh, m_dofZ);
	// rotational dofs
	scatter(U, mesh, m_dofU);
	scatter(U, mesh, m_dofV);
	scatter(U, mesh, m_dofW);

	// make sure the prescribed displacements are fullfilled
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive()) dc.Update();
	}

	// enforce the linear constraints
	// TODO: do we really have to do this? Shouldn't the algorithm
	// already guarantee that the linear constraints are satisfied?
	FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
	if (LCM.LinearConstraints() > 0)
	{
		LCM.Update();
	}

	// Update the spatial nodal positions
	// Don't update rigid nodes since they are already updated
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid == -1)
			node.m_rt = node.m_r0 + node.get_vec3d(m_dofX, m_dofY, m_dofZ);
	}
}

//-----------------------------------------------------------------------------
double FECGSolidSolver::LineSearchCG(double s)
{
	//  s is now passed from the solver routine instead of defaulting to 1.0
	double smin = s;

	double FA, FB, FC, AA, AB, r, r0;
	bool failed = false;

	// max nr of line search iterations
	int nmax = m_LSiter;
	int n = 0;

	// initial energy
	FA = m_ui*m_R0;
	AA = 0.0;
	r0 = FA;

	double rmin = fabs(FA);

	vector<double> ul(m_ui.size());

	// so we can set AA = 0 and FA= r0
	// AB=s and we need to evaluate FB (called r1)
	// n is a count of the number of linesearch attempts

	// calculate residual at this point, reducing s if necessary
	do
	{
		// Update geometry using the initial guess s
		vcopys(ul, m_ui, s);
		failed = false;
		try
		{
			Update(ul);
			Residual(m_R1);
		}
		catch (...)
		{
			//					printf("reducing s at initial evaluation");
			failed = true;
			s = 0.1*s;
		}
	} while (failed == true);

	// calculate energies
	FB = m_ui*m_R1;
	AB = s;

	if (FB<rmin){
		rmin = FB;
		smin = s;
	}

	do
	{
		// make sure that r1 does not happen to be really close to zero,
		// since in that case we won't find any better solution.
		if (fabs(FB) < 1.e-20) r = 0;  // we've hit the converged solution and don't need to do any more
		else r = fabs(FB / r0);

		if (r > m_LStol)	// we need to search and find a better value of s
		{
			if (FB == FA) s = (AA + AB) * 1000;
			else {
				s = (AA*FB - AB*FA) / (FB - FA);
				s = min(s, 100 * max(AA, AB));
			}
			// calculate residual at this point, reducing s if necessary
			do
			{
				// Update geometry using the initial guess s
				vcopys(ul, m_ui, s);
				failed = false;
				try
				{
					Update(ul);
					Residual(m_R1);
				}
				catch (...)
				{
					//					printf("reducing s at FC");
					failed = true;
					s = 0.1*s;
				}
			} while ((failed == true) && (s>m_LSmin));

			// calculate energies
			FC = m_ui*m_R1;
			r = fabs(FC / r0);

			if (fabs(FC)>100 * min(fabs(FA), fabs(FB)))  //  it was a bad guess and we need to go back a bit
			{
				s = 0.1*s;

				// calculate residual at this point, reducing s if necessary
				do
				{
					// Update geometry using the initial guess s
					vcopys(ul, m_ui, s);
					failed = false;
					try
					{
						Update(ul);
						Residual(m_R1);
					}
					catch (...)
					{
						//					printf("reducing s after bad guess");
						failed = true;
						s = 0.1*s;
					}
				} while (failed == true);

				// calculate energies
				FC = m_ui*m_R1;
				r = fabs(FC / r0);
			}

			if (fabs(FA)<fabs(FB)) // use the new value and the closest of the previous ones
			{
				FB = FC;
				AB = s;
			}
			else
			{
				FA = FC;
				AA = s;
			}

			++n;
		}
	} while ((r > m_LStol) && (n < nmax));


	if (n >= nmax)
	{
		// max nr of iterations reached.
		// we choose the line step that reached the smallest energy
		s = smin;
		// calculate residual at this point, reducing s if necessary
		do
		{
			// Update geometry using the initial guess s
			vcopys(ul, m_ui, s);
			failed = false;
			try
			{
				Update(ul);
				Residual(m_R1);
			}
			catch (...)
			{
				//					printf("reducing s after failed line search");
				failed = true;
				s = 0.1*s;
			}
		} while (failed == true);
	}

	return s;
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FECGSolidSolver::Residual(vector<double>& R)
{
	TRACK_TIME("residual");

	// get the time information
	const FETimeInfo& tp = m_fem.GetTime();

	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// setup the global vector
	FEResidualVector RHS(GetFEModel(), R, m_Fr);

	// zero rigid body reaction forces
	FEMechModel& fem = static_cast<FEMechModel&>(m_fem);
	FERigidSystem& rigid = *fem.GetRigidSystem();
	int NRB = rigid.Objects();
	for (int i = 0; i<NRB; ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		RB.m_Fr = RB.m_Mr = vec3d(0, 0, 0);
	}

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// calculate the internal (stress) forces
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// calculate the body forces
	for (int j = 0; j<m_fem.BodyLoads(); ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
		{
			for (int i = 0; i<pbf->Domains(); ++i)
			{
				FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(*pbf->Domain(i));
				dom.BodyForce(RHS, *pbf);
			}
		}
	}

	// calculate inertial forces for dynamic problems
	if (m_fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC) InertialForces(RHS);

	// calculate forces due to surface loads
	int nsl = m_fem.SurfaceLoads();
	for (int i = 0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->Residual(tp, RHS);
	}

	// calculate contact forces
	if (m_fem.SurfacePairConstraints() > 0)
	{
		ContactForces(RHS);
	}

	// calculate nonlinear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	NonLinearConstraintForces(RHS, tp);

	// forces due to point constraints
	//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->Residual(this, R);

	// add model loads
	int NML = m_fem.ModelLoads();
	for (int i = 0; i<NML; ++i)
	{
		FEModelLoad& mli = *m_fem.ModelLoad(i);
		if (mli.IsActive())
		{
			mli.Residual(RHS, tp);
		}
	}

	// set the nodal reaction forces
	// TODO: Is this a good place to do this?
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		node.m_Fr = vec3d(0, 0, 0);

		int n;
		if ((n = -node.m_ID[m_dofX] - 2) >= 0) node.m_Fr.x = -m_Fr[n];
		if ((n = -node.m_ID[m_dofY] - 2) >= 0) node.m_Fr.y = -m_Fr[n];
		if ((n = -node.m_ID[m_dofZ] - 2) >= 0) node.m_Fr.z = -m_Fr[n];
	}

	// increase RHS counter
	m_nrhs++;

	return true;
}

//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FECGSolidSolver::UpdateModel()
{
	FEMesh& mesh = m_fem.GetMesh();
	const FETimeInfo& tp = m_fem.GetTime();

	// update the stresses on all domains
	for (int i = 0; i<mesh.Domains(); ++i) mesh.Domain(i).Update(tp);
}

//-----------------------------------------------------------------------------
//! Update contact interfaces.
void FECGSolidSolver::UpdateContact()
{
	// Update all contact interfaces
	const FETimeInfo& tp = m_fem.GetTime();
	for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->Update(m_niter, tp);
	}
}

//-----------------------------------------------------------------------------
//! Update nonlinear constraints
void FECGSolidSolver::UpdateConstraints()
{
	const FETimeInfo& tp = m_fem.GetTime();

	// Update all nonlinear constraints
	for (int i = 0; i<m_fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* pci = m_fem.NonlinearConstraint(i);
		if (pci->IsActive()) pci->Update(m_niter, tp);
	}
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
bool FECGSolidSolver::Augment()
{
	const FETimeInfo& tp = m_fem.GetTime();

	// Assume we will pass (can't hurt to be optimistic)
	bool bconv = true;

	// Do contact augmentations
	if (m_fem.SurfacePairConstraints() > 0)
	{
		// loop over all contact interfaces
		for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
		{
			FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
			if (pci->IsActive()) bconv = (pci->Augment(m_naug, tp) && bconv);
		}
	}

	// do nonlinear constraint augmentations
	int n = m_fem.NonlinearConstraints();
	for (int i = 0; i<n; ++i)
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		if (plc->IsActive()) bconv = plc->Augment(m_naug, tp) && bconv;
	}

	// do incompressibility multipliers for 3Field domains
	FEMesh& mesh = m_fem.GetMesh();
	int ND = mesh.Domains();
	for (int i = 0; i<ND; ++i)
	{
		FE3FieldElasticSolidDomain* pd = dynamic_cast<FE3FieldElasticSolidDomain*>(&mesh.Domain(i));
		if (pd) bconv = (pd->Augment(m_naug) && bconv);
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FECGSolidSolver::ContactForces(FEGlobalVector& R)
{
	const FETimeInfo& tp = m_fem.GetTime();
	for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->Residual(R, tp);
	}
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces 
void FECGSolidSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
	int N = m_fem.NonlinearConstraints();
	for (int i = 0; i<N; ++i)
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->Residual(R, tp);
	}
}

//-----------------------------------------------------------------------------
//! calculates the concentrated nodal forces

void FECGSolidSolver::NodalForces(vector<double>& F, const FETimeInfo& tp)
{
	// zero nodal force vector
	zero(F);

	// loop over nodal loads
	int NNL = m_fem.NodalLoads();
	for (int i = 0; i<NNL; ++i)
	{
		const FENodalLoad& fc = *m_fem.NodalLoad(i);
		if (fc.IsActive())
		{
			int dof = fc.GetDOF();

			int N = fc.Nodes();
			for (int j=0; j<N; ++j)
			{
				int nid = fc.NodeID(j);

				// get the nodal load value
				double f = fc.NodeValue(j);

				// assemble into residual
				AssembleResidual(nid, dof, f, F);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the inertial forces for dynamic problems

void FECGSolidSolver::InertialForces(FEGlobalVector& R)
{
	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// allocate F
	vector<double> F(3 * mesh.Nodes());
	zero(F);

	// calculate F
    double dt = m_fem.GetTime().timeIncrement;
	double a = 1.0 / (m_beta*dt);
	double b = a / dt;
	double c = 1.0 - 0.5 / m_beta;
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vec3d& rt = node.m_rt;
		vec3d& rp = node.m_rp;
		vec3d& vp = node.m_vp;
		vec3d& ap = node.m_ap;

		F[3 * i] = b*(rt.x - rp.x) - a*vp.x + c * ap.x;
		F[3 * i + 1] = b*(rt.y - rp.y) - a*vp.y + c * ap.y;
		F[3 * i + 2] = b*(rt.z - rp.z) - a*vp.z + c * ap.z;
	}

	// now multiply F with the mass matrix
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(nd));
		dom.InertialForces(R, F);
	}
}

//-----------------------------------------------------------------------------
// \todo I'd like to do something different with this. Right now, if a nodal load
//       it applied to a rigid body, the load has to be translated to a force and 
//       torque applied to the rigid body. Perhaps we should really define two types
//       of nodal loads, one for the deformable body and for the rigid body. This can
//       be done in a pre-processor phase. That way, standard assembly routines can be
//       used to assemble to loads into the global vector.
void FECGSolidSolver::AssembleResidual(int node_id, int dof, double f, vector<double>& R)
{
	// get the mesh
	FEMechModel& fem = static_cast<FEMechModel&>(m_fem);
	FEMesh& mesh = fem.GetMesh();
	FERigidSystem& rigid = *fem.GetRigidSystem();

	// get the equation number
	FENode& node = mesh.Node(node_id);
	int n = node.m_ID[dof];

	// assemble into global vector
	if (n >= 0) R[n] += f;
	else if (node.m_rid >= 0)
	{
		// this is a rigid body node
		FERigidBody& RB = *rigid.Object(node.m_rid);

		// get the relative position
		vec3d a = node.m_rt - RB.m_rt;

		int* lm = RB.m_LM;
		if (dof == m_dofX)
		{
			if (lm[0] >= 0) R[lm[0]] += f;
			if (lm[4] >= 0) R[lm[4]] += a.z*f;
			if (lm[5] >= 0) R[lm[5]] += -a.y*f;
		}
		else if (dof == m_dofY)
		{
			if (lm[1] >= 0) R[lm[1]] += f;
			if (lm[3] >= 0) R[lm[3]] += -a.z*f;
			if (lm[5] >= 0) R[lm[5]] += a.x*f;
		}
		else if (dof == m_dofZ)
		{
			if (lm[2] >= 0) R[lm[2]] += f;
			if (lm[3] >= 0) R[lm[3]] += a.y*f;
			if (lm[4] >= 0) R[lm[4]] += -a.x*f;
		}
	}
}

//-----------------------------------------------------------------------------
//! Updates the rigid body data
void FECGSolidSolver::UpdateRigidBodies(vector<double>& ui)
{
	// get the number of rigid bodies
	FEMechModel& fem = static_cast<FEMechModel&>(m_fem);
	FERigidSystem& rigid = *fem.GetRigidSystem();
	const int NRB = rigid.Objects();

	// first calculate the rigid body displacement increments
	for (int i = 0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *rigid.Object(i);
		int *lm = RB.m_LM;
		double* du = RB.m_du;

		if (RB.m_prb == 0)
		{
			for (int j = 0; j<6; ++j)
			{
				du[j] = (lm[j] >= 0 ? m_Ui[lm[j]] + ui[lm[j]] : 0);
			}
		}
	}

	// for prescribed displacements, the displacement increments are evaluated differently
	// TODO: Is this really necessary? Why can't the ui vector contain the correct values?
	const int NRD = rigid.PrescribedBCs();
	for (int i = 0; i<NRD; ++i)
	{
		FERigidBodyDisplacement& dc = *rigid.PrescribedBC(i);
		if (dc.IsActive())
		{
			FERigidBody& RB = *rigid.Object(dc.id);
			if (RB.m_prb == 0)
			{
				RB.m_du[dc.bc] = (dc.lc < 0 ? 0 : dc.Value() - RB.m_Up[dc.bc]);
			}
		}
	}

	// update the rigid bodies
	for (int i = 0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *rigid.Object(i);
		double* du = RB.m_du;

		// This is the "new" update algorithm which addressesses a couple issues
		// with the old method, namely that prescribed rotational dofs aren't update correctly.
		// Unfortunately, it seems to produce worse convergence in some cases, especially with line search
		// and it doesn't work when rigid bodies are used in a hierarchy
		if (RB.m_prb) du = RB.m_dul;
		RB.m_Ut[0] = RB.m_Up[0] + du[0];
		RB.m_Ut[1] = RB.m_Up[1] + du[1];
		RB.m_Ut[2] = RB.m_Up[2] + du[2];
		RB.m_Ut[3] = RB.m_Up[3] + du[3];
		RB.m_Ut[4] = RB.m_Up[4] + du[4];
		RB.m_Ut[5] = RB.m_Up[5] + du[5];

		RB.m_rt = RB.m_r0 + vec3d(RB.m_Ut[0], RB.m_Ut[1], RB.m_Ut[2]);

		vec3d Rt(RB.m_Ut[3], RB.m_Ut[4], RB.m_Ut[5]);
		RB.SetRotation(quatd(Rt));
	}

	// we need to update the position of rigid nodes
	rigid.UpdateMesh();

	// Since the rigid nodes are repositioned we need to update the displacement DOFS
	FEMesh& mesh = m_fem.GetMesh();
	int N = mesh.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0)
		{
			vec3d ut = node.m_rt - node.m_r0;
			node.set_vec3d(m_dofX, m_dofY, m_dofZ, ut);
		}
	}
}
