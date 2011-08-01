#include "stdafx.h"
#include <math.h>
#include "FESolver.h"
#include "FEBioLib/FERigid.h"
#include "FESolidSolver.h"
#include "log.h"
#include "Interrupt.h"

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

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FESolidSolver::SolveStep
//  Solves a single time step
//  This function mainly calls the Quasin routine 
//  and deals with exceptions that require the immediate termination of
//	quasi-Newton iterations.
//

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
		if (m_fem.GetDebugFlag()) m_fem.m_plot->Write(m_fem);
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

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolidSolver::PrepStep
// Prepares the data for the first QN iteration. 
//

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
		m_fem.m_mesh.Node(i).m_cp = m_fem.m_mesh.Node(i).m_ct;
	}

	// apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
	NodalForces(m_Fn);

	// apply prescribed displacements
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_bfgs.m_ui;
	zero(ui);
	int neq = m_fem.m_neq;
	for (i=0; i<(int) m_fem.m_DC.size(); ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;

			double dq = s*m_fem.GetLoadCurve(lc)->Value();

			int I;

			FENode& node = m_fem.m_mesh.Node(n);

			switch (bc)
			{
			case 0: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<neq) 
					ui[I] = dq - (node.m_rt.x - node.m_r0.x);
				break;
			case 1: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<neq) 
					ui[I] = dq - (node.m_rt.y - node.m_r0.y); 
				break;
			case 2: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<neq) 
					ui[I] = dq - (node.m_rt.z - node.m_r0.z); 
				break;
			// ---> TODO: move to the FEPoroSolidSolver
			case 6: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<neq) 
					ui[I] = dq - node.m_pt; 
				break;
			case 11: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<neq) 
					ui[I] = dq - node.m_ct; 
				break;
			// --->
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
			}
		}
	}

	// initialize rigid bodies
	int nrb = m_fem.m_RB.size();
	for (i=0; i<nrb; ++i)
	{
		FERigidBody& RB = m_fem.m_RB[i];

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

	// calculate local rigid displacements
	for (i=0; i<(int) m_fem.m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *m_fem.m_RDC[i];
		FERigidBody& RB = m_fem.m_RB[DC.id];
		if (RB.m_bActive && DC.IsActive())
		{
			int I = DC.bc;
			int lc = DC.lc;
			if (lc > 0)
			{
				RB.m_dul[I] = DC.sf*m_fem.GetLoadCurve(lc-1)->Value() - RB.m_Ut[DC.bc];
			}
		}
	}

	// calculate global rigid displacements
	for (i=0; i<(int) m_fem.m_RB.size(); ++i)
	{
		FERigidBody& RB = m_fem.m_RB[i];
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

	// store rigid displacements in Ui vector
	for (i=0; i<(int) m_fem.m_RB.size(); ++i)
	{
		FERigidBody& RB = m_fem.m_RB[i];
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
		FERigidBody& RB = m_fem.m_RB[FC.id];
		if (RB.m_bActive && FC.IsActive())
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
	if (m_fem.ContactInterfaces() > 0) m_fem.UpdateContact();

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	FEMaterialPoint::dt = m_fem.m_pStep->m_dt;
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
//
bool FESolidSolver::Quasin(double time)
{
	int i;
	double s;

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

	// prepare for the first iteration
	PrepStep(time);

	// check for CTRL+C interruption before we do any work
	if (m_fem.m_bInterruptable)
	{
		Interruption itr;
		if (itr.m_bsig)
		{
			itr.m_bsig = false;
			itr.interrupt();
		}
	}

	// calculate initial stiffness matrix
	if (ReformStiffness() == false) return false;

	// calculate initial residual
	if (Residual(m_bfgs.m_R0) == false) return false;

	m_bfgs.m_R0 += m_Fd;

	// TODO: I can check here if the residual is zero.
	// If it is than there is probably no force acting on the system
	// if (m_R0*m_R0 < eps) bconv = true;

//	double r0 = m_R0*m_R0;

	Logfile::MODE oldmode;

	clog.printf("\n===== beginning time step %d : %lg =====\n", m_fem.m_pStep->m_ntimesteps+1, m_fem.m_ftime);

	// loop until converged or when max nr of reformations reached
	do
	{
		oldmode = clog.GetMode();
		if ((m_fem.m_pStep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.m_pStep->GetPrintLevel() != FE_PRINT_NEVER)) clog.SetMode(Logfile::FILE_ONLY);

		clog.printf(" %d\n", m_niter+1);
		clog.SetMode(oldmode);

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
		if (m_bfgs.m_LStol > 0) s = m_bfgs.LineSearch();
		else
		{
			s = 1;

			// Update geometry
			Update(m_bfgs.m_ui);

			// calculate residual at this point
			Residual(m_bfgs.m_R1);
		}

		// update total displacements
		int neq = m_Ui.size();
		for (i=0; i<neq; ++i) m_Ui[i] += s*m_bfgs.m_ui[i];

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
		oldmode = clog.GetMode();
		if ((m_fem.m_pStep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.m_pStep->GetPrintLevel() != FE_PRINT_NEVER)) clog.SetMode(Logfile::FILE_ONLY);

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
				clog.printbox("WARNING", "Problem is diverging. Stiffness matrix will now be reformed");
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
			if (breform)
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
		else if (m_fem.m_pStep->m_baugment)
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
				if (m_fem.m_pStep->m_psolver->m_bfgs.m_maxups == 0)
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
		if (m_fem.m_bInterruptable)
		{
			Interruption itr;
			if (itr.m_bsig)
			{
				itr.m_bsig = false;
				itr.interrupt();
			}
		}
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

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolidSolver::ReformStiffness
// Reforms a stiffness matrix and factorizes it
//

bool FESolidSolver::ReformStiffness()
{
	// first, let's make sure we have not reached the max nr of reformations allowed
	if (m_nref >= m_bfgs.m_maxref) throw MaxStiffnessReformations();

	// recalculate the shape of the stiffness matrix if necessary
	if (m_breshape)
	{
		// TODO: I don't think I need to update here
//		if (m_fem.m_bcontact) m_fem.UpdateContact();

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

