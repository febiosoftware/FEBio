#include "FEThermoElasticSolver.h"
#include "FEBioMech/FEResidualVector.h"
#include <FEBioMech/FEElasticDomain.h>
#include "FEThermoElasticSolidDomain.h"
#include <FECore/FERigidBody.h>
#include <FECore/log.h>

#ifdef WIN32
	#include <float.h>
	#define ISNAN(x) _isnan(x)
#endif

#ifdef LINUX
	#ifdef CENTOS
		#define ISNAN(x) isnan(x)
	#else
		#define ISNAN(x) std::isnan(x)
	#endif
#endif

#ifdef __APPLE__
#include <math.h>
#define ISNAN(x) isnan(x)
#endif

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEThermoElasticSolver, FESolidSolver)
	ADD_PARAMETER(m_Ttol         , FE_PARAM_DOUBLE, "Ttol"        );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEThermoElasticSolver::FEThermoElasticSolver(FEModel* pfem) : FESolidSolver(pfem)
{
	m_Ttol = 0.01;

	m_ndeq = 0;
	m_nteq = 0;

	// Allocate degrees of freedom
	// (X,Y,Z) dofs are allocated in base class
	DOFS& dofs = pfem->GetDOFS();
	int nT = dofs.AddVariable("temperature");
	dofs.SetDOFName(nT, 0, "T");

	// get the temperature degree of freedom index
	m_dofT = m_fem.GetDOFS().GetDOF("T");
}

//-----------------------------------------------------------------------------
// destructor
FEThermoElasticSolver::~FEThermoElasticSolver()
{
	// nothing to do here
}

//-----------------------------------------------------------------------------
// This function initializes the solver data. It is called before the Solve routine.
bool FEThermoElasticSolver::Init()
{
	// we need a non-symmetric stiffness matrix
	m_bsymm = false;

	// call the base class version first
	if (FESolidSolver::Init() == false) return false;

	// allocate vectors
	if (m_ndeq > 0)
	{
		m_di.assign(m_ndeq, 0);
		m_Di.assign(m_ndeq, 0);
	}
	if (m_nteq > 0) {
		m_ti.assign(m_nteq, 0);
		m_Ti.assign(m_nteq, 0);

		// we need to fill the total displacement vector m_Ut
		// TODO: I need to find an easier way to do this
		FEMesh& mesh = m_fem.GetMesh();
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);

			// temperature dofs
			int n = node.m_ID[m_dofT]; if (n >= 0) m_Ut[n] = node.get(m_dofT);
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize equations.
//! This function is called by the FEAnalysis::Init function. Base class does
//! most of the work. Just need to count the number of displacement and temperature
//! equations.
bool FEThermoElasticSolver::InitEquations()
{
	// base class does most of the work
	FESolidSolver::InitEquations();

	// determine the nr of displacement and temperature equations
	FEMesh& mesh = m_fem.GetMesh();
	m_ndeq = m_nteq = 0;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& n = mesh.Node(i);
		if (n.m_ID[m_dofX] != -1) m_ndeq++;
		if (n.m_ID[m_dofY] != -1) m_ndeq++;
		if (n.m_ID[m_dofZ] != -1) m_ndeq++;
		if (n.m_ID[m_dofT] != -1) m_nteq++;
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
void FEThermoElasticSolver::PrepStep(double time)
{
	zero(m_Ti);
	zero(m_Di);

	FESolidSolver::PrepStep(time);
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
bool FEThermoElasticSolver::Quasin(double time)
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

	// temperature convergence norms data
	double	normTi;		// initial temperature norm
	double	normT;		// current temperature norm
	double	normt;		// incremement temperature norm

	// initialize flags
	bool bconv = false;		// convergence flag
	bool breform = false;	// reformation flag

	// get the current step
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// prepare for the first iteration
	PrepStep(time);

	// do minor iterations callbacks
	m_fem.DoCallback(CB_MINOR_ITERS);

	// calculate initial stiffness matrix
	FETimePoint tp = m_fem.GetTime();
	if (ReformStiffness(tp) == false) return false;

	// calculate initial residual
	if (Residual(m_R0) == false) return false;

	// add the stiffness-contributions to the RHS (from displacement BC's)
	m_R0 += m_Fd;

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

		// extract the displacement increments
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
			GetTemperatureData(m_ti, m_ui);

			// set initial norm
			if (m_niter == 0) normTi = fabs(m_ti*m_ti);

			// update total pressure
			for (i=0; i<m_nteq; ++i) m_Ti[i] += s*m_ti[i];

			// calculate norms
			normT = m_Ti*m_Ti;
			normt = (m_ti*m_ti)*(s*s);

			// check convergence
			if ((m_Ttol > 0) && (normt > (m_Ttol*m_Ttol)*normT)) bconv = false;
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
		felog.printf("\t   temperature      %15le %15le %15le \n", normTi, normt ,(m_Ttol*m_Ttol)*normT );

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
				normTi = normt;
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
// Extract the displacements from the incremental update vector.
void FEThermoElasticSolver::GetDisplacementData(vector<double> &di, const vector<double> &ui)
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
// Extract the temperature from the incremental update vector.
void FEThermoElasticSolver::GetTemperatureData(vector<double> &ti, const vector<double> &ui)
{
	int N = m_fem.GetMesh().Nodes(), nid, m = 0;
	zero(ti);
	for (int i=0; i<N; ++i)
	{
		FENode& n = m_fem.GetMesh().Node(i);
		nid = n.m_ID[m_dofT];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			ti[m++] = ui[nid];
			assert(m <= (int) ti.size());
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEThermoElasticSolver::Residual(vector<double>& R)
{
	TimerTracker t(m_RHSTime);

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

	// calculate internal stress force
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// calculate internal thermal work
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEThermoElasticSolidDomain* pdom = dynamic_cast<FEThermoElasticSolidDomain*>(&mesh.Domain(i));
		if (pdom) pdom->InternalThermalWork(R);
	}

    // calculate the body forces
    for (i=0; i<mesh.Domains(); ++i)
    {
        FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
        for (int j=0; j<m_fem.BodyLoads(); ++j)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
            if (pbf) dom.BodyForce(RHS, *pbf);
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
void FEThermoElasticSolver::Update(vector<double>& ui)
{
	// update all elastic data
	FESolidSolver::Update(ui);

	// update nodal temperatures
	FEMesh& mesh = m_fem.GetMesh();
	const int NN = mesh.Nodes();
	for (int i=0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);

		// current temperature = initial + total at prev conv step + total increment so far + current increment
		int n = node.m_ID[m_dofT];
		if (n >= 0) node.set(m_dofT, m_Ut[n] + m_Ui[n] + ui[n]);
	}

	// make sure the prescribed displacements are fullfilled
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive()) dc.Update();
	}
}
