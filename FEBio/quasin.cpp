#include "stdafx.h"
#include "FESolver.h"
#include <math.h>
#include "stack.h"

#ifdef WIN32
	#include <float.h>
	#define ISNAN(x) _isnan(x)
#endif

#ifdef LINUX
	#include <math.h>
	#define ISNAN(x) isnan(x)
#endif


///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FESolver::SolveStep
//  Solves a single time step
//  This function mainly calls the Quasin routine 
//  and deals with exceptions that require the immediate termination of
//	quasi-Newton iterations.
//

bool FESolver::SolveStep(double time)
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
		m_log.printbox("ERROR","Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng, e.m_vol);
		if (m_fem.m_debug) m_fem.m_plot.Write(m_fem);
		return false;
	}
	catch (MaxStiffnessReformations)
	{
		// max nr of reformations is reached
		m_log.printbox("ERROR", "Max nr of reformations reached.");
		return false;
	}
	catch (ForceConversion)
	{
		// user forced conversion of problem
		m_log.printbox("WARNING", "User forced conversion.\nSolution might not be stable.");
		return true;
	}
	catch (IterationFailure)
	{
		// user caused a forced iteration failure
		m_log.printbox("WARNING", "User forced iteration failure.");
		return false;
	}
	catch (ZeroLinestepSize)
	{
		// a zero line step size was detected
		m_log.printbox("ERROR", "Zero line step size.");
		return false;
	}
	catch (EnergyDiverging)
	{
		// problem was diverging after stiffness reformation
		m_log.printbox("ERROR", "Problem diverging uncontrollably.");
		return false;
	}

	return bret;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolver::PrepStep
// Prepares the data for the first QN iteration. 
//

void FESolver::PrepStep(double time)
{
	int i, j;

	// initialize counters
	m_niter = 0;	// nr of iterations
	m_nrhs  = 0;	// nr of RHS evaluations
	m_nref  = 0;	// nr of stiffness reformations
	m_nups	= 0;	// nr of stiffness updates between reformations
	m_naug  = 0;	// nr of augmentations

	// evaluate load curve values at current time
	for (i=0; i<m_fem.LoadCurves(); ++i) m_fem.GetLoadCurve(i)->Evaluate(time);

	// zero total displacements/pressures
	m_Ui.zero();

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	for (i=0; i<m_fem.m_mesh.Nodes(); ++i)
	{
		m_fem.m_mesh.Node(i).m_rp = m_fem.m_mesh.Node(i).m_rt;
		m_fem.m_mesh.Node(i).m_vp = m_fem.m_mesh.Node(i).m_vt;
		m_fem.m_mesh.Node(i).m_ap = m_fem.m_mesh.Node(i).m_at;
	}

	// apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
	NodalForces(m_Fn);

	// update the body acceleration values
	// the "-" sign is to be consistent with NIKE3D's convention
	if (m_fem.m_BF[0].lc >= 0) m_fem.m_acc.x = -m_fem.GetLoadCurve(m_fem.m_BF[0].lc)->Value()*m_fem.m_BF[0].s;
	if (m_fem.m_BF[1].lc >= 0) m_fem.m_acc.y = -m_fem.GetLoadCurve(m_fem.m_BF[1].lc)->Value()*m_fem.m_BF[1].s;
	if (m_fem.m_BF[2].lc >= 0) m_fem.m_acc.z = -m_fem.GetLoadCurve(m_fem.m_BF[2].lc)->Value()*m_fem.m_BF[2].s;

	// apply prescribed displacements
	// we save the prescribed displacements increments in the ui vector
	m_ui.zero();
	for (i=0; i<m_fem.m_DC.size(); ++i)
	{
		if (m_fem.m_DC[i].IsActive())
		{
			int n    = m_fem.m_DC[i].node;
			int lc   = m_fem.m_DC[i].lc;
			int bc   = m_fem.m_DC[i].bc;
			double s = m_fem.m_DC[i].s;

			double dq = s*m_fem.GetLoadCurve(lc)->Value();

			int I;

			FENode& node = m_fem.m_mesh.Node(n);

			switch (bc)
			{
			case 0: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<m_fem.m_neq) 
					m_ui[I] = dq - (node.m_rt.x - node.m_r0.x);
				break;
			case 1: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<m_fem.m_neq) 
					m_ui[I] = dq - (node.m_rt.y - node.m_r0.y); 
				break;
			case 2: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<m_fem.m_neq) 
					m_ui[I] = dq - (node.m_rt.z - node.m_r0.z); 
				break;
			case 6: 
				I = -node.m_ID[bc]-2;
				if (I>=0 && I<m_fem.m_neq) 
					m_ui[I] = dq - node.m_pt; 
				break;
			case 20:
				{
					vec3d dr = node.m_r0;
					dr.x = 0; dr.unit(); dr *= dq;

					I = -node.m_ID[1]-2;
					if (I>=0 && I<m_fem.m_neq) 
						m_ui[I] = dr.y - (node.m_rt.y - node.m_r0.y); 
					I = -node.m_ID[2]-2;
					if (I>=0 && I<m_fem.m_neq) 
						m_ui[I] = dr.z - (node.m_rt.z - node.m_r0.z); 
				}
				break;
			}
		}
	}

	// apply prescribed rigid body constraints
	for (i=0; i<m_fem.m_nrb; ++i)
	{
		FERigidBody& RB = m_fem.m_RB[i];

		// clear reaction forces
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);

		for (j=0; j<6; ++j)
		{
			int I = -RB.m_LM[j]-2;
			int lc = RB.m_bc[j];
			if ((I >= 0) && (lc > 0))
			{
				m_ui[I] = m_fem.GetLoadCurve(lc-1)->Value() - RB.m_Ut[j];
			}
		}

		RB.m_rp = RB.m_rt;
		RB.m_qp = RB.m_qt;
		RB.m_Up[0] = RB.m_Ut[0];
		RB.m_Up[1] = RB.m_Ut[1];
		RB.m_Up[2] = RB.m_Ut[2];
		RB.m_Up[3] = RB.m_Ut[3];
		RB.m_Up[4] = RB.m_Ut[4];
		RB.m_Up[5] = RB.m_Ut[5];

		// apply prescribed rigid body forces
		FERigid* pm = dynamic_cast<FERigid*>(m_fem.GetMaterial(RB.m_mat));
		for (j=0; j<6; ++j)
		{
			int lc = pm->m_fc[j];
			int I  = RB.m_LM[j];
			if ((I>=0) && (lc >= 0))
			{
				double f = m_fem.GetLoadCurve(lc)->Value()*pm->m_fs[j];
				m_Fn[I] += f;

				switch (j)
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
	if (m_fem.m_bcontact) m_fem.UpdateContact();

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	FEMaterialPoint::dt = m_fem.m_pStep->m_dt;
	FEMaterialPoint::time = m_fem.m_ftime;

	FEMesh& mesh = m_fem.m_mesh;
	for (i=0; i<mesh.SolidElements(); ++i)
	{
		FESolidElement& el = mesh.SolidElement(i);
		int n = el.GaussPoints();
		for (j=0; j<n; ++j) el.m_State[j]->Init(false);
	}
	for (i=0; i<mesh.ShellElements(); ++i)
	{
		FEShellElement& el = mesh.ShellElement(i);
		int n = el.GaussPoints();
		for (j=0; j<n; ++j) el.m_State[j]->Init(false);
	}


	// intialize the stresses
	// TODO: is this a good place to update the stresses?
	// Perhaps I should place this back in the residual routine?
	m_fem.UpdateStresses();
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FESolver::Quasin
// Implements the BFGS algorithm to solve the nonlinear FE equations.
// The details of this implementation of the BFGS method can be found in:
//   "Finite Element Procedures", K.J. Bathe, p759 and following
//

bool FESolver::Quasin(double time)
{
	int i;
	double s;

	// initialize flags
	bool bconv = false;		// convergence flag
	bool breform = false;	// reformation flag

	// prepare for the first iteration
	PrepStep(time);

	// calculate initial stiffness matrix
	if (ReformStiffness() == false) return false;

	// calculate initial residual
	if (Residual(m_R0) == false) return false;

	m_R0 += m_Fd;

	// TODO: I can check here if the residual is zero.
	// If it is than there is probably no force acting on the system
	// if (m_R0*m_R0 < eps) bconv = true;

//	double r0 = m_R0*m_R0;

	Logfile::MODE oldmode;

	// loop until converged or when max nr of reformations reached
	do
	{
		oldmode = m_log.GetMode();
		if (m_fem.m_pStep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) m_log.SetMode(Logfile::FILE_ONLY);
		m_log.printf(" %d\n", m_niter+1);
		m_log.SetMode(oldmode);

		// assume we'll converge. 
		bconv = true;

		// solve the equations
		SolveEquations(m_ui, m_R0);

		// check for nans
		if (m_fem.GetDebugFlag())
		{
			double du = m_ui*m_ui;
			if (ISNAN(du)) throw NANDetected();
		}

		// set initial convergence norms
		if (m_niter == 0)
		{
			m_normRi = fabs(m_R0*m_R0);
			m_normEi = fabs(m_ui*m_R0);
			m_normUi = fabs(m_ui*m_ui);
			m_normEm = m_normEi;
		}

		// perform a linesearch
		// the geometry is also updated in the line search
		s = LineSearch();

		// update total displacements
		for (i=0; i<m_fem.m_neq; ++i) m_Ui[i] += s*m_ui[i];

		// calculate norms
		m_normR1 = m_R1*m_R1;
		m_normu  = (m_ui*m_ui)*(s*s);
		m_normU  = m_Ui*m_Ui;
		m_normE1 = s*fabs(m_ui*m_R1);
	
		// check residual norm
		if (m_normR1 > m_Rtol*m_normRi) bconv = false;	

		// check displacement norm
		if (m_normu  > (m_Dtol*m_Dtol)*m_normU ) bconv = false;

		// check energy norm
		if (m_normE1 > m_Etol*m_normEi) bconv = false;

		// check linestep size
		if (s < m_LSmin) bconv = false;

		// check energy divergence
		if (m_normE1 > m_normEm) bconv = false;

		// print convergence summary
		oldmode = m_log.GetMode();
		if (m_fem.m_pStep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) m_log.SetMode(Logfile::FILE_ONLY);
		m_log.printf(" Nonlinear solution status: time= %lg\n", time); 
		m_log.printf("\tstiffness updates             = %d\n", m_nups);
		m_log.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		m_log.printf("\tstiffness matrix reformations = %d\n", m_nref);
		m_log.printf("\tstep from line search         = %lf\n", s);
		m_log.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		m_log.printf("\t   residual         %15le %15le %15le \n", m_normRi, m_normR1, m_Rtol*m_normRi);
		m_log.printf("\t   energy           %15le %15le %15le \n", m_normEi, m_normE1, m_Etol*m_normEi);
		m_log.printf("\t   displacement     %15le %15le %15le \n", m_normUi, m_normu ,(m_Dtol*m_Dtol)*m_normU );
		m_log.SetMode(oldmode);

		// check if we have converged. 
		// If not, calculate the BFGS update vectors
		if (bconv == false)
		{
			if ((m_niter == 0) && (m_normR1 < 1.0e-20))
			{
				// check for almost zero-residual on the first iteration
				// this might be an indication that there is no force on the system
				m_log.printbox("WARNING", "No force acting on the system.");
				bconv = true;
			}
			else if (s < m_LSmin)
			{
				// check for zero linestep size
				m_log.printbox("WARNING", "Zero linestep size. Stiffness matrix will now be reformed");
				breform = true;
			}
			else if (m_normE1 > m_normEm)
			{
				// check for diverging
				m_log.printbox("WARNING", "Problem is diverging. Stiffness matrix will now be reformed");
				m_normEm = m_normE1;
				m_normEi = m_normE1;
				m_normRi = m_normR1;
				breform = true;
			}
			else
			{
				// If we havn't reached max nr of BFGS updates
				// do an update
				if (!breform)
				{
					if (m_nups < m_maxups-1)
					{
						if (BFGSUpdate(s) == false)
						{
							// Stiffness update has failed.
							// this might be due a too large condition number
							// or the update was no longer positive definite.
							m_log.printbox("WARNING", "The BFGS update has failed.\nStiffness matrix will now be reformed.");
							breform = true;
						}
					}
					else
					{
						// we've reached the max nr of BFGS updates, so
						// we need to do a stiffness reformation
						breform = true;

						// print a warning only if the user did not intent full-Newton
						if (m_maxups > 0)
							m_log.printbox("WARNING", "Max nr of iterations reached.\nStiffness matrix will now be reformed.");

					}
				}
			}	

			// zero displacement increments
			// we must set this to zero before the reformation
			// because we assume that the prescribed displacements are stored 
			// in the m_ui vector.
			m_ui.zero();

			// reform stiffness matrices if necessary
			if (breform)
			{
				m_log.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);

				// reform the matrix
				if (ReformStiffness() == false) break;
	
				// reset reformation flag
				breform = false;
			}

			// copy last calculated residual
			m_R0 = m_R1;
		}
		else if (m_fem.m_pStep->m_baugment)
		{
			// we have converged, so let's see if the augmentations have converged as well

			m_log.printf("\n........................ augmentation # %d\n", m_naug+1);

			// do the augmentations
			bconv = Augment();

			// update counter
			++m_naug;
	
			// If we havn't converged we prepare for the next iteration
			if (!bconv) 
			{
				// Since the Lagrange multipliers have changed, we can't just copy 
				// the last residual but have to recalculate the residual
				// we also recalculate the stresses in case we are doing augmentations
				// for incompressible materials
				m_fem.UpdateStresses();
				Residual(m_R0);
			}
		}
	
		// increase iteration number
		m_niter++;

		// let's flush the logfile to make sure the last output will not get lost
		m_log.flush();

		// check for CTRL+C interruption
		if (m_bsig)
		{
			m_bsig = false;
			interrupt();
		}
	}
	while (bconv == false);

	// when converged, 
	// print a convergence summary to the log file
	if (bconv)
	{
		Logfile::MODE mode = m_log.SetMode(Logfile::FILE_ONLY);
		m_log.printf("\nconvergence summary\n");
		m_log.printf("    number of iterations   : %d\n", m_niter);
		m_log.printf("    number of reformations : %d\n", m_nref);

		m_log.SetMode(mode);
	}

	// if converged we update the total displacements
	if (bconv)
	{
		m_Ut += m_Ui;
	}

	return bconv;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolver::ReformStiffness
// Reforms a stiffness matrix and factorizes it
//

bool FESolver::ReformStiffness()
{
	// first, let's make sure we have not reached the max nr of reformations allowed
	if (m_nref >= m_maxref) throw MaxStiffnessReformations();

	// recalculate the shape of the stiffness matrix if necessary
	if (m_breshape)
	{
		if (m_fem.m_bcontact) m_fem.UpdateContact();

		// reshape the stiffness matrix
		if (!CreateStiffness(m_niter == 0)) return false;

		// reset reshape flag, except for contact
		m_breshape = (m_fem.m_bcontact ? true : false);
	}

	// calculate the stiffness matrices
	bool bret = StiffnessMatrix();

	if (bret)
	{
		m_SolverTime.start();
		{
			// factorize the stiffness matrix
			m_psolver->Factor(*m_pK);
		}
		m_SolverTime.stop();

		// increase total nr of reformations
		m_nref++;

		// reset bfgs update counter
		m_nups = 0;
	}

	return bret;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FESolver::LineSearch
// Performs a linesearch on a NR iteration
// The description of this method can be found in:
//    "Nonlinear Continuum Mechanics for Finite Element Analysis", 	Bonet & Wood.
//
// TODO: Find a different way to update the deformation based on the ls.
// For instance, define a di so that ui = s*di. Also, define the 
// position of the nodes at the previous iteration.

double FESolver::LineSearch()
{
	double s = 1.0;
	double smin = s;

	double a, A, B, D;
	double r0, r1, r;

	// max nr of line search iterations
	int nmax = m_LSiter;
	int n = 0;

	// initial energy
	r0 = m_ui*m_R0;

	double rmin = fabs(r0);

	if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
	{
		m_fem.m_log.printf("\nEntering line search\n");
		m_fem.m_log.printf("       STEPSIZE     INITIAL        CURRENT         REQUIRED\n");
	}

	do
	{
		// Update geometry
		Update(m_ui, s);

		// calculate residual at this point
		Residual(m_R1);

		// make sure we are still in a valid range
		if (s < m_LSmin) 
		{
			// it appears that we are not converging
			// I found in the NIKE3D code that when this happens,
			// the line search step is simply set to 0.5.
			// so let's try it here too
			s = 0.5;

			// reupdate  
			Update(m_ui, s);

			// recalculate residual at this point
			Residual(m_R1);

			// return and hope for the best
			break;
		}

		// calculate energies
		r1 = m_ui*m_R1;

		if ((n==0) || (fabs(r1) < rmin))
		{
			smin = s;
			rmin = fabs(r1);
		}

		// make sure that r1 does not happen to be really close to zero,
		// since in that case we won't find any better solution.
		if (fabs(r1) < 1.e-20) r = 0;
		else r = fabs(r1/r0);

		if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
		{
			m_fem.m_log.printf("%15lf%15lg%15lg%15lg\n", s, fabs(r0), fabs(r1), fabs(r0*m_LStol));
		}

		if (r > m_LStol)
		{
			// calculate the line search step
			a = r0/r1;

			A = 1 + a*(s-1);
			B = a*(s*s);
			D = B*B - 4*A*B;

			// update the line step
			if (D >= 0) 
			{
				s = (B + sqrt(D))/(2*A);
				if (s < 0) s = (B - sqrt(D))/(2*A);
				if (s < 0) s = 0;
			}
			else 
			{
				s = 0.5*B/A;
			}

			++n;
		}
	}
	while ((r > m_LStol) && (n < nmax));

	if (n >= nmax)
	{
		// max nr of iterations reached.
		// we choose the line step that reached the smallest energy
		s = smin;
		Update(m_ui, s);
		Residual(m_R1);
	}

	if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
	{
		if ((r < m_LStol) && (n < nmax))
		{
			m_fem.m_log.printf("------------------------------------------------------------\n");
			m_fem.m_log.printf("%15lf%15lg%15lg%15lg\n\n", s, fabs(r0), fabs(r1), fabs(r0*m_LStol));
		}
		else if (n >= nmax)
		{
			m_fem.m_log.printf("Line search failed: max iterations reached.\n\n");
		}
	}

	return s;
}


///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolver::SolveEquations
// This function solves a system of equations using the BFGS update vectors
// The variable m_nups keeps track of how many updates have been made so far.
// 

void FESolver::SolveEquations(vector<double>& x, vector<double>& b)
{
	int i, j;
	double *vi, *wi, vr, wr;

	// get the nr of equations
	int neq = m_fem.m_neq;

	// create temporary storage
	static vector<double> tmp;
	tmp = b;

	// loop over all update vectors
	for (i=m_nups-1; i>=0; --i)
	{
		vi = m_V[i];
		wi = m_W[i];

		wr = 0;
		for (j=0; j<neq; j++) wr += wi[j]*tmp[j];
		for (j=0; j<neq; j++) tmp[j] += vi[j]*wr;
	}

	// perform a backsubstitution
	m_SolverTime.start();
	{
		m_psolver->Solve(*m_pK, x, tmp);
	}
	m_SolverTime.stop();

	// loop again over all update vectors
	for (i=0; i<m_nups; ++i)
	{
		vi = m_V[i];
		wi = m_W[i];

		vr = 0;
		for (j=0; j<neq; ++j) vr += vi[j]*x[j];
		for (j=0; j<neq; ++j) x[j] += wi[j]*vr;
	}
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolver::SolveEquations
// This function solves a system of equations using the BFGS update vectors
// The variable m_nups keeps track of how many updates have been made so far.
// 

void FESolver::SolveEquations(matrix& x, matrix& b)
{
	int i, j, k;
	double *vi, *wi, vr, wr;

	// get the nr of equations
	int neq = m_fem.m_neq;

	// create temporary storage
	matrix tmp;
	tmp = b;

	// get the nr of rows
	// note that this is actually the nr of columns since we 
	// store the transpose of matrices in order to have
	// column major ordering
	int nr = x.rows();

	// loop over all update vectors
	for (i=m_nups-1; i>=0; --i)
	{
		vi = m_V[i];
		wi = m_W[i];

		for (k=0; k<nr; ++k)
		{
			wr = 0;
			for (j=0; j<neq; j++) wr += wi[j]*tmp[k][j];
			for (j=0; j<neq; j++) tmp[k][j] += vi[j]*wr;
		}
	}

	// perform a backsubstitution
	m_SolverTime.start();
	{
		m_psolver->Solve(*m_pK, x, tmp);
	}
	m_SolverTime.stop();

	// loop again over all update vectors
	for (i=0; i<m_nups; ++i)
	{
		vi = m_V[i];
		wi = m_W[i];

		for (k=0; k<nr; ++k)
		{
			vr = 0;
			for (j=0; j<neq; ++j) vr += vi[j]*x[k][j];
			for (j=0; j<neq; ++j) x[k][j] += wi[j]*vr;
		}
	}
}

//-----------------------------------------------------------------------------
//! This function performs a BFGS stiffness update.
//! The last line search step is input to this function.
//! This function performs the update assuming the stiffness matrix
//! is positive definite. In the case that this is not the case
//! the function returns false. The function also returns false if 
//! the condition number is too big. A too big condition number might
//! be an indication of an ill-conditioned matrix and the update should
//! not be performed.

bool FESolver::BFGSUpdate(double s)
{
	int i;
	double dg, dh,dgi, c, r;
	double *vn, *wn;

	int neq = m_fem.m_neq;

	// calculate the BFGS update vectors
	for (i=0; i<neq; ++i)	
	{
		m_D[i] = s*m_ui[i];
		m_G[i] = m_R0[i] - m_R1[i];
		m_H[i] = m_R0[i]*s;
	}

	dg = m_D*m_G;
	dh = m_D*m_H;
	dgi = 1.0 / dg;
	r = dg/dh;

	// check to see if this is still a pos definite update
//	if (r <= 0) 
//	{
//		return false;
//	}

	// calculate the condition number
//	c = sqrt(r);
	c = sqrt(fabs(r));

	// make sure c is less than the the maximum.
	if (c > m_cmax) return false;

	vn = m_V[m_nups];
	wn = m_W[m_nups];

	// TODO: There might be a bug here. Check signs!
	for (i=0; i<neq; ++i)	
	{
		vn[i] = -m_H[i]*c - m_G[i];
		wn[i] = m_D[i]*dgi;
	}

	// increment update counter
	++m_nups;

	return true;
}
