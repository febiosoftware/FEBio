#include "FEMultiphasicSolver.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FEBiphasicSolidDomain.h"
#include "FEBiphasicSoluteDomain.h"
#include "FEMultiphasicDomain.h"
#include "FETriphasicDomain.h"
#include "FEBiphasicSolidDomain.h"
#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FESlidingInterfaceMP.h"
#include "FEBioMech/FEPressureLoad.h"
#include "FEBioMech/FEResidualVector.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "FECore/FEGlobalMatrix.h"

#ifdef WIN32
	#include <float.h>
	#define ISNAN(x) _isnan(x)
#endif

#ifdef LINUX
	#define ISNAN(x) std::isnan(x)
#endif

#ifdef __APPLE__
#include <math.h>
#define ISNAN(x) isnan(x)
#endif

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEMultiphasicSolver, FESolidSolver2)
	ADD_PARAMETER(m_Ptol         , FE_PARAM_DOUBLE, "ptol"        );
	ADD_PARAMETER(m_Ctol         , FE_PARAM_DOUBLE, "ctol"        );
	ADD_PARAMETER(m_bsymm        , FE_PARAM_BOOL  , "symmetric_biphasic");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEMultiphasicSolver::FEMultiphasicSolver(FEModel* pfem) : FESolidSolver2(pfem)
{
	m_Ctol = 0.01;
    
	m_bsymm = false; // assume non-symmetric stiffness matrix by default
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures.
//
bool FEMultiphasicSolver::Init()
{
	// initialize base class
	if (FESolidSolver2::Init() == false) return false;

	// allocate poro-vectors
	assert(m_ndeq > 0);
	m_di.assign(m_ndeq, 0);
	m_Di.assign(m_ndeq, 0);

//	assert(m_npeq > 0);
	if (m_npeq > 0) {
		m_pi.assign(m_npeq, 0);
		m_Pi.assign(m_npeq, 0);

		// we need to fill the total displacement vector m_Ut
		// TODO: I need to find an easier way to do this
		FEMesh& mesh = m_fem.GetMesh();
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);

			// pressure dofs
			int n = node.m_ID[DOF_P]; if (n >= 0) m_Ut[n] = node.m_pt;
		}
	}

    // get number of DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
	// allocate concentration-vectors
	m_ci.assign(MAX_CDOFS,vector<double>(0,0));
	m_Ci.assign(MAX_CDOFS,vector<double>(0,0));
	for (int i=0; i<MAX_CDOFS; ++i) {
		m_ci[i].assign(m_nceq[i], 0);
		m_Ci[i].assign(m_nceq[i], 0);
	}
	
	// we need to fill the total displacement vector m_Ut
	// TODO: I need to find an easier way to do this
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		
		// concentration dofs
		for (int j=0; j<(int)m_nceq.size(); ++j) {
			if (m_nceq[j]) {
				int n = node.m_ID[DOF_C+j];
				if (n >= 0) m_Ut[n] = node.m_ct[j];
			}
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize equations
bool FEMultiphasicSolver::InitEquations()
{
	// base class does most of the work
	FESolidSolver2::InitEquations();

	// determined the nr of pressure and concentration equations
	FEMesh& mesh = m_fem.GetMesh();
	m_ndeq = m_npeq = 0;
	
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& n = mesh.Node(i);
		if (n.m_ID[DOF_X] != -1) m_ndeq++;
		if (n.m_ID[DOF_Y] != -1) m_ndeq++;
		if (n.m_ID[DOF_Z] != -1) m_ndeq++;
		if (n.m_ID[DOF_P] != -1) m_npeq++;
	}
	
	// determine the nr of concentration equations
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetCDOFS();
    m_nceq.assign(MAX_CDOFS, 0);
	
    // get number of DOFS
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& n = mesh.Node(i);
		for (int j=0; j<(int)m_nceq.size(); ++j)
			if (n.m_ID[DOF_C+j] != -1) m_nceq[j]++;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first QN iteration. 
//!
void FEMultiphasicSolver::PrepStep(double time)
{
	for (int j=0; j<(int)m_nceq.size(); ++j) if (m_nceq[j]) zero(m_Ci[j]);

	zero(m_Pi);
	zero(m_Di);

	FESolidSolver2::PrepStep(time);
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
//! The details of this implementation of the BFGS method can be found in:
//!   "Finite Element Procedures", K.J. Bathe, p759 and following
//!
bool FEMultiphasicSolver::Quasin(double time)
{
	int i, j;
	double s;

	// convergence norms
	double	normR1;		// residual norm
	double	normE1;		// energy norm
	double	normD;		// displacement norm
	double	normd;		// displacement increment norm
	double	normRi;		// initial residual norm
	double	normEi;		// initial energy norm
	double	normEm;		// max energy norm
	double	normDi;		// initial displacement norm

	// poro convergence norms data
	double	normPi;		// initial pressure norm
	double	normP;		// current pressure norm
	double	normp;		// incremement pressure norm

    // get number of DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
	// solute convergence data
	vector<double>	normCi(MAX_CDOFS);	// initial concentration norm
	vector<double>	normC(MAX_CDOFS);	// current concentration norm
	vector<double>	normc(MAX_CDOFS);	// incremement concentration norm

	// initialize flags
	bool bconv = false;		// convergence flag
	bool breform = false;	// reformation flag

	// get the current step
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// make sure this is poro-solute problem
	assert(pstep->GetType() == FE_MULTIPHASIC);

	// prepare for the first iteration
	PrepStep(time);

	// do minor iterations callbacks
	m_fem.DoCallback(CB_MINOR_ITERS);

	// calculate initial stiffness matrix
	FETimePoint tp = m_fem.GetTime();
	if (ReformStiffness(tp) == false) return false;

	// calculate initial residual
	if (Residual(m_R0) == false) return false;

	// Add the "reaction forces" from prescribed dofs.
	// This vector is created by bringing the stiffness contributions
	// from the prescribed dofs to the RHS. 
	m_R0 += m_Fd;

	// TODO: I can check here if the residual is zero.
	// If it is than there is probably no force acting on the system
	// if (m_R0*m_R0 < eps) bconv = true;

//	double r0 = m_R0*m_R0;

	Logfile::MODE oldmode;

	felog.printf("\n===== beginning time step %d : %lg =====\n", pstep->m_ntimesteps+1, m_fem.m_ftime);

	// loop until converged or when max nr of reformations reached
	do
	{
		oldmode = felog.GetMode();
		if ((m_fem.GetCurrentStep()->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.GetCurrentStep()->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::FILE_ONLY);

		felog.printf(" %d\n", m_niter+1);
		felog.SetMode(oldmode);

		// assume we'll converge. 
		bconv = true;

		// solve the equations
		m_SolverTime.start();
		{
			m_bfgs.SolveEquations(m_ui, m_R0);
		}
		m_SolverTime.stop();

		// check for nans
		double du = m_ui*m_ui;
		if (ISNAN(du)) throw NANDetected();

		// extract the pressure increments
		GetDisplacementData(m_di, m_ui);

		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_R0*m_R0);
			normEi = fabs(m_ui*m_R0);
			normDi = fabs(m_di*m_di);
			normEm = normEi;
		}

		// perform a linesearch
		// the geometry is also updated in the line search
		if (m_LStol > 0) s = LineSearch(1.0);
		else
		{
			s = 1;

			// Update geometry
			Update(m_ui);

			// calculate residual at this point
			Residual(m_R1);
		}

		// update all degrees of freedom
		for (i=0; i<m_neq; ++i) m_Ui[i] += s*m_ui[i];

		// update displacements
		for (i=0; i<m_ndeq; ++i) m_Di[i] += s*m_di[i];

		// calculate norms
		normR1 = m_R1*m_R1;
		normd  = (m_di*m_di)*(s*s);
		normD  = m_Di*m_Di;
		normE1 = s*fabs(m_ui*m_R1);

		// check residual norm
		if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;	

		// check displacement norm
		if ((m_Dtol > 0) && (normd  > (m_Dtol*m_Dtol)*normD )) bconv = false;

		// check energy norm
		if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;

		// check linestep size
		if ((m_LStol > 0) && (s < m_LSmin)) bconv = false;

		// check energy divergence
		if (normE1 > normEm) bconv = false;

		// check poroelastic convergence
		{
			// extract the pressure increments
			GetPressureData(m_pi, m_ui);

			// set initial norm
			if (m_niter == 0) normPi = fabs(m_pi*m_pi);

			// update total pressure
			for (i=0; i<m_npeq; ++i) m_Pi[i] += s*m_pi[i];

			// calculate norms
			normP = m_Pi*m_Pi;
			normp = (m_pi*m_pi)*(s*s);

			// check convergence
			if ((m_Ptol > 0) && (normp > (m_Ptol*m_Ptol)*normP)) bconv = false;
		}

		// check solute convergence
		{
			// extract the concentration increments
			for (j=0; j<(int)m_nceq.size(); ++j) {
				if (m_nceq[j]) {
					GetConcentrationData(m_ci[j], m_ui,j);
					
					// set initial norm
					if (m_niter == 0)
						normCi[j] = fabs(m_ci[j]*m_ci[j]);
					
					// update total concentration
					for (i=0; i<m_nceq[j]; ++i) m_Ci[j][i] += s*m_ci[j][i];
					
					// calculate norms
					normC[j] = m_Ci[j]*m_Ci[j];
					normc[j] = (m_ci[j]*m_ci[j])*(s*s);
					
				}
			}
			
			// check convergence
			if (m_Ctol > 0) {
				for (j=0; j<(int)m_nceq.size(); ++j)
					if (m_nceq[j]) bconv = bconv && (normc[j] <= (m_Ctol*m_Ctol)*normC[j]);
			}
		}

		// print convergence summary
		oldmode = felog.GetMode();
		if ((m_fem.GetCurrentStep()->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.GetCurrentStep()->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::FILE_ONLY);

		felog.printf(" Nonlinear solution status: time= %lg\n", time); 
		felog.printf("\tstiffness updates             = %d\n", m_bfgs.m_nups);
		felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
		felog.printf("\tconvergence norms :        INITIAL         CURRENT         REQUIRED\n");
		felog.printf("\t residual               %15le %15le %15le\n", normRi, normR1, m_Rtol*normRi);
		felog.printf("\t energy                 %15le %15le %15le\n", normEi, normE1, m_Etol*normEi);
		felog.printf("\t displacement           %15le %15le %15le\n", normDi, normd ,(m_Dtol*m_Dtol)*normD );
		felog.printf("\t fluid pressure         %15le %15le %15le\n", normPi, normp ,(m_Ptol*m_Ptol)*normP );
		for (j=0; j<(int)m_nceq.size(); ++j) {
			if (m_nceq[j])
				felog.printf("\t solute %d concentration %15le %15le %15le\n", j+1, normCi[j], normc[j] ,(m_Ctol*m_Ctol)*normC[j] );
		}

		felog.SetMode(oldmode);

		// check if we have converged. 
		// If not, calculate the BFGS update vectors
		if (bconv == false)
		{
			if ((normR1 < m_Rmin))
			{
				// check for almost zero-residual on the first iteration
				// this might be an indication that there is no force on the system
				felog.printbox("WARNING", "No force acting on the system.");
				bconv = true;
			}
			else if (s < m_LSmin)
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
				normDi = normd;
				normPi = normp;
				for (j=0; j<(int)m_nceq.size(); ++j)
					if (m_nceq[j]) normCi[j] = normc[j];
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
						if (m_bfgs.Update(s, m_ui, m_R0, m_R1) == false)
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
			zero(m_ui);

			// reform stiffness matrices if necessary
			if (breform)
			{
				felog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);

				// reform the matrix
				if (ReformStiffness(tp) == false) break;
	
				// reset reformation flag
				breform = false;
			}

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
				UpdateStresses();
				Residual(m_R0);

				// reform the matrix if we are using full-Newton
				if (m_bfgs.m_maxups == 0)
				{
					felog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);
					if (ReformStiffness(tp) == false) break;
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
		m_Ut += m_Ui;
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEMultiphasicSolver::Residual(vector<double>& R)
{
	int i;
	double dt = m_fem.GetCurrentStep()->m_dt;

	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// setup global RHS vector
	FEResidualVector RHS(GetFEModel(), R, m_Fr);

	// zero rigid body reaction forces
	FERigidSystem& rs = *m_fem.GetRigidSystem();
	int NRB = rs.Objects();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& RB = *rs.Object(i);
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);
	}

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// internal stress work
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	if (m_fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
	{
		for (i=0; i<mesh.Domains(); ++i) 
		{
			FEDomain& dom = mesh.Domain(i);
			FEBiphasicSolidDomain*  pbd = dynamic_cast<FEBiphasicSolidDomain* >(&dom);
			FEBiphasicSoluteDomain* pbs = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
			FETriphasicDomain*      ptd = dynamic_cast<FETriphasicDomain*     >(&dom);
			FEMultiphasicDomain*    pmd = dynamic_cast<FEMultiphasicDomain*   >(&dom);

			if (pbd)
			{
				pbd->InternalFluidWorkSS(R, dt);
			} 
			else if (pbs)
			{
				pbs->InternalFluidWorkSS(R, dt);
				pbs->InternalSoluteWorkSS(R, dt);
			} 
			else if (ptd)
			{
				ptd->InternalFluidWorkSS(R, dt);
				ptd->InternalSoluteWorkSS(R, dt);
			}
			else if (pmd)
			{
				pmd->InternalFluidWorkSS (R, dt);
				pmd->InternalSoluteWorkSS(R, dt);
			}
		}
	}
	else
	{
		for (i=0; i<mesh.Domains(); ++i) 
		{
			FEDomain& dom = mesh.Domain(i);
			FEBiphasicSolidDomain*  pbd = dynamic_cast<FEBiphasicSolidDomain* >(&dom);
			FEBiphasicSoluteDomain* pbs = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
			FETriphasicDomain*      ptd = dynamic_cast<FETriphasicDomain*     >(&dom);
			FEMultiphasicDomain*    pmd = dynamic_cast<FEMultiphasicDomain*   >(&dom);

			if (pbd)
			{
				pbd->InternalFluidWork(R, dt);
			}
			else if (pbs)
			{
				pbs->InternalFluidWork(R, dt);
				pbs->InternalSoluteWork(R, dt);
			}
			else if (ptd)
			{
				ptd->InternalFluidWork(R, dt);
				ptd->InternalSoluteWork(R, dt);
			}
			else if (pmd)
			{
 				pmd->InternalFluidWork (R, dt);
				pmd->InternalSoluteWork(R, dt);
			}
		}
	}

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

	// get the time information
	FETimePoint tp = m_fem.GetTime();

	// calculate linear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	NonLinearConstraintForces(RHS, tp);

	// add model loads
	int NML = m_fem.ModelLoads();
	for (i=0; i<NML; ++i)
	{
		FEModelLoad& mli = *m_fem.ModelLoad(i);
		if (mli.IsActive())
		{
			mli.Residual(RHS, tp);
		}
	}

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
//! Calculates global stiffness matrix.

bool FEMultiphasicSolver::StiffnessMatrix(const FETimePoint& tp)
{
	// get the stiffness matrix
	SparseMatrix& K = *m_pK;

	// zero stiffness matrix
	K.zero();

	// zero the residual adjustment vector
	zero(m_Fd);

	// element stiffness matrix
	matrix ke;

	// nodal degrees of freedom
	int i, j, I;

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// calculate the stiffness matrix for each domain
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	bool bsymm = m_bsymm;
	double dt = pstep->m_dt;
	if (pstep->m_nanalysis == FE_STEADY_STATE)
	{
		for (i=0; i<mesh.Domains(); ++i) 
		{
			FEDomain& dom = mesh.Domain(i);
			FEElasticSolidDomain*   pde = dynamic_cast<FEElasticSolidDomain*  >(&dom);
			FEBiphasicSolidDomain*  pbd = dynamic_cast<FEBiphasicSolidDomain* >(&dom);
			FEBiphasicSoluteDomain* pbs = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
			FETriphasicDomain*      ptd = dynamic_cast<FETriphasicDomain*     >(&dom);
			FEMultiphasicDomain*    pmd = dynamic_cast<FEMultiphasicDomain*   >(&dom);

			if      (pbd) pbd->StiffnessMatrixSS(this, bsymm, tp.dt);
			else if (pbs) pbs->StiffnessMatrixSS(this, bsymm, tp);
			else if (ptd) ptd->StiffnessMatrixSS(this, bsymm, tp);
			else if (pmd) pmd->StiffnessMatrixSS(this, bsymm, tp);
            else if (pde) pde->StiffnessMatrix(this);
		}
	}
	else
	{
		for (i=0; i<mesh.Domains(); ++i) 
		{
			FEDomain& dom = mesh.Domain(i);
			FEElasticSolidDomain*   pde = dynamic_cast<FEElasticSolidDomain*  >(&dom);
			FEBiphasicSolidDomain*  pbd = dynamic_cast<FEBiphasicSolidDomain* >(&dom);
			FEBiphasicSoluteDomain* pbs = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
			FETriphasicDomain*      ptd = dynamic_cast<FETriphasicDomain*     >(&dom);
			FEMultiphasicDomain*    pmd = dynamic_cast<FEMultiphasicDomain*   >(&dom);

			if      (pbd) pbd->StiffnessMatrix(this, bsymm, tp.dt);
			else if (pbs) pbs->StiffnessMatrix(this, bsymm, tp);
			else if (ptd) ptd->StiffnessMatrix(this, bsymm, tp);
			else if (pmd) pmd->StiffnessMatrix(this, bsymm, tp);
            else if (pde) pde->StiffnessMatrix(this);
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
	NonLinearConstraintStiffness(tp);

	// we still need to set the diagonal elements to 1
	// for the prescribed rigid body dofs.
	FERigidSystem& rs = *m_fem.GetRigidSystem();
	int NRB = rs.Objects();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& rb = *rs.Object(i);
		for (j=0; j<6; ++j)
			if (rb.m_LM[j] < -1)
			{
				I = -rb.m_LM[j]-2;
				K.set(I,I, 1);
			}
	}

	// let's check the stiffness matrix for zero diagonal elements
	int neq = K.Size();
	for (i=0; i<neq; ++i)
	{
		if (K.diag(i) == 0) throw ZeroDiagonal(-1, -1);
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolver::GetDisplacementData(vector<double> &di, vector<double> &ui)
{
	int N = m_fem.GetMesh().Nodes(), nid, m = 0;
	zero(di);
	for (int i=0; i<N; ++i)
	{
		FENode& n = m_fem.GetMesh().Node(i);
		nid = n.m_ID[DOF_X];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
		nid = n.m_ID[DOF_Y];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
		nid = n.m_ID[DOF_Z];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
	}
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolver::GetPressureData(vector<double> &pi, vector<double> &ui)
{
	int N = m_fem.GetMesh().Nodes(), nid, m = 0;
	zero(pi);
	for (int i=0; i<N; ++i)
	{
		FENode& n = m_fem.GetMesh().Node(i);
		nid = n.m_ID[DOF_P];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			pi[m++] = ui[nid];
			assert(m <= (int) pi.size());
		}
	}
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolver::GetConcentrationData(vector<double> &ci, vector<double> &ui, const int sol)
{
	int N = m_fem.GetMesh().Nodes(), nid, m = 0;
	zero(ci);
	for (int i=0; i<N; ++i)
	{
		FENode& n = m_fem.GetMesh().Node(i);
		nid = n.m_ID[DOF_C+sol];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			ci[m++] = ui[nid];
			assert(m <= (int) ci.size());
		}
	}
}


//-----------------------------------------------------------------------------
//! Update the model's kinematic data. This is overriden from FEBiphasicSolver so
//! that solute data is updated
void FEMultiphasicSolver::UpdateKinematics(vector<double>& ui)
{
	// first update all solid-mechanics kinematics
	FESolidSolver2::UpdateKinematics(ui);

	// update poroelastic data
	UpdatePoro(ui);

	// update solute-poroelastic data
	UpdateSolute(ui);
}

//-----------------------------------------------------------------------------
//! Updates the poroelastic data
void FEMultiphasicSolver::UpdatePoro(vector<double>& ui)
{
	int i, n;

	FEMesh& mesh = m_fem.GetMesh();
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update nodal pressures
		n = node.m_ID[DOF_P];
		if (n >= 0) node.m_pt = 0 + m_Ut[n] + m_Ui[n] + ui[n];
	}

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update velocities
		node.m_vt  = (node.m_rt - node.m_rp) / pstep->m_dt;
	}
}

//-----------------------------------------------------------------------------
//! Updates the solute data
void FEMultiphasicSolver::UpdateSolute(vector<double>& ui)
{
	int i, j, n;
	
	FEMesh& mesh = m_fem.GetMesh();
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	
    // get number of DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
	// update solute data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		
		// update nodal concentration
		for (j=0; j<MAX_CDOFS; ++j) {
			n = node.m_ID[DOF_C+j];
//			if (n >= 0) node.m_ct[j] = 0 + m_Ut[n] + m_Ui[n] + ui[n];
			// Force the concentrations to remain positive
			if (n >= 0) {
				node.m_ct[j] = 0 + m_Ut[n] + m_Ui[n] + ui[n];
				if (node.m_ct[j] < 0) {
					node.m_ct[j] = 0;
				}
			}
		}
	}
	
	// update solute data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		
		// update velocities
		node.m_vt  = (node.m_rt - node.m_rp) / pstep->m_dt;
	}
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolver::UpdateContact()
{
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// mark all free-draining surfaces
	for (int i=0; i<m_fem.SurfacePairInteractions(); ++i) 
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairInteraction(i));

		FESlidingInterface2* psi2 = dynamic_cast<FESlidingInterface2*>(pci);
		if (psi2) psi2->MarkFreeDraining();
		FESlidingInterface3* psi3 = dynamic_cast<FESlidingInterface3*>(pci);
		if (psi3) psi3->MarkAmbient();
		FESlidingInterfaceMP* psiMP = dynamic_cast<FESlidingInterfaceMP*>(pci);
		if (psiMP) psiMP->MarkAmbient();
	}

	// Update all contact interfaces
	FESolidSolver2::UpdateContact();

	// set free-draining boundary conditions
	for (int i=0; i<m_fem.SurfacePairInteractions(); ++i) 
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairInteraction(i));

		FESlidingInterface2* psi2 = dynamic_cast<FESlidingInterface2*>(pci);
		if (psi2) psi2->SetFreeDraining();
		FESlidingInterface3* psi3 = dynamic_cast<FESlidingInterface3*>(pci);
		if (psi3) psi3->SetAmbient();
		FESlidingInterfaceMP* psiMP = dynamic_cast<FESlidingInterfaceMP*>(pci);
		if (psiMP) psiMP->SetAmbient();
	}
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEMultiphasicSolver::Serialize(DumpFile& ar)
{
	FESolidSolver2::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_Ptol;
		ar << m_ndeq << m_npeq;
	}
	else
	{
		ar >> m_Ptol;
		ar >> m_ndeq >> m_npeq;
	}
    
	if (ar.IsSaving())
	{
		ar << m_Ctol;
		ar << m_nceq;
	}
	else
	{
		ar >> m_Ctol;
		ar >> m_nceq;
	}
}
