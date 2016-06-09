#include "FEBiphasicSolver.h"
#include "FEBiphasicSolidDomain.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FESlidingInterfaceBiphasic.h"
#include "FEBioMech/FEElasticDomain.h"
#include "FEBioMech/FEPressureLoad.h"
#include "FEBioMech/FEResidualVector.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/sys.h"
#include <FECore/FEModel.h>
#include <FECore/FEModelLoad.h>
#include <FECore/BC.h>
#include <FECore/FERigidSystem.h>
#include <FECore/FEAnalysis.h>

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEBiphasicSolver, FESolidSolver2)
	ADD_PARAMETER(m_Ptol         , FE_PARAM_DOUBLE, "ptol"        );
	ADD_PARAMETER(m_bsymm        , FE_PARAM_BOOL  , "symmetric_biphasic");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEBiphasicSolver::FEBiphasicSolver(FEModel* pfem) : FESolidSolver2(pfem)
{
	m_Ptol = 0.01;
	m_ndeq = 0;
	m_npeq = 0;

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int varP = dofs.AddVariable("fluid pressure");
	dofs.SetDOFName(varP, 0, "p");

	// get pressure dof
	m_dofP = pfem->GetDOFIndex("p");
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures.
//
bool FEBiphasicSolver::Init()
{
	// initialize base class
	if (FESolidSolver2::Init() == false) return false;

	// allocate poro-vectors
//	assert(m_ndeq > 0);
	m_di.assign(m_ndeq, 0);
	m_Di.assign(m_ndeq, 0);

//	assert(m_npeq > 0);
	if (m_npeq > 0) {
		m_pi.assign(m_npeq, 0);
		m_Pi.assign(m_npeq, 0);

		// we need to fill the total displacement vector m_Ut
		// TODO: I need to find an easier way to do this
		FEMesh& mesh = m_fem.GetMesh();
		int i, n;
		for (i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);

			// pressure dofs
			n = node.m_ID[m_dofP]; if (n >= 0) m_Ut[n] = node.get(m_dofP);
		}
	}

	return true;
}


//-----------------------------------------------------------------------------
//! Initialize equations
bool FEBiphasicSolver::InitEquations()
{
	// base class does most of the work
	FESolidSolver2::InitEquations();

	int i;

	// determined the nr of pressure and concentration equations
	FEMesh& mesh = m_fem.GetMesh();
	m_ndeq = m_npeq = 0;
	
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& n = mesh.Node(i);
		if (n.m_ID[m_dofX] != -1) m_ndeq++;
		if (n.m_ID[m_dofY] != -1) m_ndeq++;
		if (n.m_ID[m_dofZ] != -1) m_ndeq++;
		if (n.m_ID[m_dofP] != -1) m_npeq++;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first QN iteration. 
//!
//! \todo There is some more stuff in the base method that 
//!       I need to move to this method, but since it will
//!       change the order of some operations I need to make
//!       sure it won't break anything
void FEBiphasicSolver::PrepStep(const FETimeInfo& timeInfo)
{
	zero(m_Pi);
	zero(m_Di);

	FESolidSolver2::PrepStep(timeInfo);
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
//! The details of this implementation of the BFGS method can be found in:
//!   "Finite Element Procedures", K.J. Bathe, p759 and following
//!
bool FEBiphasicSolver::Quasin(double time)
{
	int i;
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

	// initialize flags
	bool bconv = false;		// convergence flag
	bool breform = false;	// reformation flag

	// get the current step
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// prepare for the first iteration
	FETimeInfo tp = m_fem.GetTime();
	PrepStep(tp);

	// calculate initial stiffness matrix
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
			m_pbfgs->SolveEquations(m_ui, m_R0);
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

		// print convergence summary
		oldmode = felog.GetMode();
		if ((m_fem.GetCurrentStep()->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.GetCurrentStep()->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::FILE_ONLY);

		felog.printf(" Nonlinear solution status: time= %lg\n", time); 
		felog.printf("\tstiffness updates             = %d\n", m_pbfgs->m_nups);
		felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
		felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		felog.printf("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
		felog.printf("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
		felog.printf("\t   displacement     %15le %15le %15le \n", normDi, normd ,(m_Dtol*m_Dtol)*normD );
		felog.printf("\t   fluid pressure   %15le %15le %15le \n", normPi, normp ,(m_Ptol*m_Ptol)*normP );

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
			else if ((normE1 > normEm) && m_bdivreform)
			{
				// check for diverging
				felog.printbox("WARNING", "Problem is diverging. Stiffness matrix will now be reformed");
				normEm = normE1;
				normEi = normE1;
				normRi = normR1;
				normDi = normd;
				normPi = normp;
				breform = true;
			}
			else
			{
				// If we havn't reached max nr of BFGS updates
				// do an update
				if (!breform)
				{
					if (m_pbfgs->m_nups < m_pbfgs->m_maxups-1)
					{
						if (m_pbfgs->Update(s, m_ui, m_R0, m_R1) == false)
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
						if (m_pbfgs->m_maxups > 0)
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
				if (m_pbfgs->m_maxups == 0)
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
//! calculates the concentrated nodal forces
void FEBiphasicSolver::NodalForces(vector<double>& F, const FETimeInfo& tp)
{
	// zero nodal force vector
	zero(F);

	// loop over nodal loads
	int NNL = m_fem.NodalLoads();
	for (int i=0; i<NNL; ++i)
	{
		const FENodalLoad& fc = *m_fem.NodalLoad(i);
		if (fc.IsActive())
		{
			int dof = fc.GetDOF();

			int N = fc.Nodes();
			for (int j=0; j<N; ++j)
			{
				int nid	= fc.NodeID(j);	// node ID

				// get the nodal load value
				double f = fc.NodeValue(j);
			
				// For pressure and concentration loads, multiply by dt
				// for consistency with evaluation of residual and stiffness matrix
				if (dof == m_dofP) f *= tp.timeIncrement;

				// assemble into residual
				AssembleResidual(nid, dof, f, F);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEBiphasicSolver::Residual(vector<double>& R)
{
	TimerTracker t(m_RHSTime);

	int i;

	// get the time information
	FETimeInfo tp = m_fem.GetTime();
	double dt = tp.timeIncrement;

	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// setup global RHS vector
	FEResidualVector RHS(GetFEModel(), R, m_Fr);

	// zero rigid body reaction forces
	FERigidSystem& rigid = *m_fem.GetRigidSystem();
	int NRB = rigid.Objects();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);
	}

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// calculate internal stress force
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// calculate internal fluid work
	if (m_fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
	{
		for (i=0; i<mesh.Domains(); ++i)
		{
			FEBiphasicSolidDomain* pdom = dynamic_cast<FEBiphasicSolidDomain*>(&mesh.Domain(i));
			if (pdom) pdom->InternalFluidWorkSS(R, dt);
		}
	}
	else
	{
		for (i=0; i<mesh.Domains(); ++i)
		{
			FEBiphasicSolidDomain* pdom = dynamic_cast<FEBiphasicSolidDomain*>(&mesh.Domain(i));
			if (pdom) pdom->InternalFluidWork(R, dt);
		}
	}

    // calculate the body forces
    for (i=0; i<mesh.Domains(); ++i)
    {
        FEBiphasicSolidDomain* pbdom = dynamic_cast<FEBiphasicSolidDomain*>(&mesh.Domain(i));
        FEElasticSolidDomain* pedom = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(i));
        for (int j=0; j<m_fem.BodyLoads(); ++j)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
            if (pbf) {
                if (pbdom) pbdom->BodyForce(RHS, *pbf);
                else if (pedom) pedom->BodyForce(RHS, *pbf);
            }
        }
    }
    
	// calculate forces due to surface loads
	int nsl = m_fem.SurfaceLoads();
	for (i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->Residual(tp, RHS);
	}

	// calculate contact forces
	if (m_fem.SurfacePairInteractions() > 0)
	{
		ContactForces(RHS);
	}

	// calculate nonlinear constraint forces
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
		if ((n = -node.m_ID[m_dofX]-2) >= 0) node.m_Fr.x = -m_Fr[n];
		if ((n = -node.m_ID[m_dofY]-2) >= 0) node.m_Fr.y = -m_Fr[n];
		if ((n = -node.m_ID[m_dofZ]-2) >= 0) node.m_Fr.z = -m_Fr[n];
	}

	// increase RHS counter
	m_nrhs++;

	return true;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FEBiphasicSolver::StiffnessMatrix(const FETimeInfo& tp)
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
	double dt = tp.timeIncrement;
	if (pstep->m_nanalysis == FE_STEADY_STATE)
	{
		for (i=0; i<mesh.Domains(); ++i) 
		{
            // Biphasic analyses may include biphasic and elastic domains
			FEBiphasicSolidDomain* pbdom = dynamic_cast<FEBiphasicSolidDomain*>(&mesh.Domain(i));
			if (pbdom) pbdom->StiffnessMatrixSS(this, bsymm, dt);
            else
			{
				FEElasticDomain* pedom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
				if (pedom) pedom->StiffnessMatrix(this);
			}
		}
	}
	else
	{
		for (i=0; i<mesh.Domains(); ++i) 
		{
            // Biphasic analyses may include biphasic and elastic domains
			FEBiphasicSolidDomain* pbdom = dynamic_cast<FEBiphasicSolidDomain*>(&mesh.Domain(i));
			if (pbdom) pbdom->StiffnessMatrix(this, bsymm, dt);
            else 
			{
				FEElasticDomain* pedom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
				if (pedom) pedom->StiffnessMatrix(this);
			}
		}
	}

    // calculate the body force stiffness matrix for each domain
	// TODO: This is not  going to work with FEDiscreteSpringDomain
    for (i=0; i<mesh.Domains(); ++i)
    {
        FEBiphasicSolidDomain* pbdom = dynamic_cast<FEBiphasicSolidDomain*>(&mesh.Domain(i));
        FEElasticSolidDomain* pedom = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(i));
        int NBL = m_fem.BodyLoads();
        for (int j=0; j<NBL; ++j)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
            if (pbf) {
                if (pbdom) pbdom->BodyForceStiffness(this, *pbf);
                else if (pedom) pedom->BodyForceStiffness(this, *pbf);
            }
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
		if ((dynamic_cast<FEPressureLoad*>(psl) == 0) || (m_fem.GetCurrentStep()->m_istiffpr != 0)) psl->StiffnessMatrix(tp, this); 
	}

	// calculate nonlinear constraint stiffness
	// note that this is the contribution of the 
	// constrainst enforced with augmented lagrangian
	NonLinearConstraintStiffness(tp);

	// we still need to set the diagonal elements to 1
	// for the prescribed rigid body dofs.
	FERigidSystem& rigid = *m_fem.GetRigidSystem();
	int NRB = rigid.Objects();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& rb = *rigid.Object(i);
		for (j=0; j<6; ++j)
			if (rb.m_LM[j] < -1)
			{
				I = -rb.m_LM[j]-2;
				K.set(I,I, 1);
			}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Update the model's kinematic data. This is overriden from FESolidSolver2 so
//! that biphasic data is updated
void FEBiphasicSolver::UpdateKinematics(vector<double>& ui)
{
	// first update all solid-mechanics kinematics
	FESolidSolver2::UpdateKinematics(ui);

	// update poroelastic data
	UpdatePoro(ui);
}

//-----------------------------------------------------------------------------
//! Updates the poroelastic data
void FEBiphasicSolver::UpdatePoro(vector<double>& ui)
{
	int i, n;

	FEMesh& mesh = m_fem.GetMesh();
	double dt = m_fem.GetTime().timeIncrement;

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update nodal pressures
		n = node.m_ID[m_dofP];
		if (n >= 0) node.set(m_dofP, 0 + m_Ut[n] + m_Ui[n] + ui[n]);
	}

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update velocities
		vec3d vt = (node.m_rt - node.m_rp) / dt;
		node.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, vt); 
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolver::UpdateContact()
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
        FESlidingInterfaceBiphasic* psib = dynamic_cast<FESlidingInterfaceBiphasic*>(pci);
        if (psib) psib->MarkFreeDraining();
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
        FESlidingInterfaceBiphasic* psib = dynamic_cast<FESlidingInterfaceBiphasic*>(pci);
        if (psib) psib->SetFreeDraining();
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolver::GetDisplacementData(vector<double> &di, vector<double> &ui)
{
	int N = m_fem.GetMesh().Nodes(), nid, m = 0;
	zero(di);
	for (int i=0; i<N; ++i)
	{
		FENode& n = m_fem.GetMesh().Node(i);
		nid = n.m_ID[m_dofX];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
		nid = n.m_ID[m_dofY];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
		nid = n.m_ID[m_dofZ];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolver::GetPressureData(vector<double> &pi, vector<double> &ui)
{
	int N = m_fem.GetMesh().Nodes(), nid, m = 0;
	zero(pi);
	for (int i=0; i<N; ++i)
	{
		FENode& n = m_fem.GetMesh().Node(i);
		nid = n.m_ID[m_dofP];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			pi[m++] = ui[nid];
			assert(m <= (int) pi.size());
		}
	}
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEBiphasicSolver::Serialize(DumpStream& ar)
{
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

	FESolidSolver2::Serialize(ar);
}
