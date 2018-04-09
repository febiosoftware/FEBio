#include "FEThermoElasticSolver.h"
#include "FEBioMech/FEResidualVector.h"
#include <FEBioMech/FEElasticDomain.h>
#include "FEThermoElasticSolidDomain.h"
#include <FECore/log.h>
#include <FECore/sys.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FEBioMech/FEBodyForce.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelLoad.h>

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

		// we need to fill the total solution vector
		FEMesh& mesh = m_fem.GetMesh();
		gather(m_Ut, mesh, m_dofT);
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
void FEThermoElasticSolver::PrepStep()
{
	zero(m_Ti);
	zero(m_Di);

	FESolidSolver::PrepStep();
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
bool FEThermoElasticSolver::Quasin()
{
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

	// get the current step
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// prepare for the first iteration
	const FETimeInfo& tp = m_fem.GetTime();
	PrepStep();

	// init QN method
	if (QNInit() == false) return false;

	// loop until converged or when max nr of reformations reached
	bool bconv = false;		// convergence flag
	do
	{
		Logfile::MODE oldmode = felog.GetMode();
		if ((m_fem.GetCurrentStep()->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.GetCurrentStep()->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::LOG_FILE);

		felog.printf(" %d\n", m_niter+1);
		felog.SetMode(oldmode);

		// assume we'll converge. 
		bconv = true;

		// solve the equations (returns line search; solution stored in m_ui)
		double s = QNSolve();

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

		// update all degrees of freedom
		for (int i=0; i<m_neq; ++i) m_Ui[i] += s*m_ui[i];

		// update displacements
		for (int i = 0; i<m_ndeq; ++i) m_Di[i] += s*m_di[i];

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
		if ((m_lineSearch->m_LStol > 0) && (s < m_lineSearch->m_LSmin)) bconv = false;

		// check energy divergence
		if (normE1 > normEm) bconv = false;

		// check poroelastic convergence
		{
			// extract the pressure increments
			GetTemperatureData(m_ti, m_ui);

			// set initial norm
			if (m_niter == 0) normTi = fabs(m_ti*m_ti);

			// update total pressure
			for (int i = 0; i<m_nteq; ++i) m_Ti[i] += s*m_ti[i];

			// calculate norms
			normT = m_Ti*m_Ti;
			normt = (m_ti*m_ti)*(s*s);

			// check convergence
			if ((m_Ttol > 0) && (normt > (m_Ttol*m_Ttol)*normT)) bconv = false;
		}

		// print convergence summary
		oldmode = felog.GetMode();
		if ((m_fem.GetCurrentStep()->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.GetCurrentStep()->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::LOG_FILE);

		felog.printf(" Nonlinear solution status: time= %lg\n", tp.currentTime); 
		felog.printf("\tstiffness updates             = %d\n", m_strategy->m_nups);
		felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_lineSearch->m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
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
				QNForceReform(true);
			}
			else if (s < m_lineSearch->m_LSmin)
			{
				// check for zero linestep size
				felog.printbox("WARNING", "Zero linestep size. Stiffness matrix will now be reformed");
				QNForceReform(true);
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
				QNForceReform(true);
			}

			// Do the QN update (This may also do a stiffness reformation if necessary)
			bool bret = QNUpdate();

			// something went wrong with the update, so we'll need to break
			if (bret == false) break;
		}
		else if (m_baugment)
		{
			// Do augmentations
			bconv = DoAugmentations();
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
	TRACK_TIME("residual");

	// get the time information
	const FETimeInfo& tp = m_fem.GetTime();

	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// setup global RHS vector
	FEResidualVector RHS(GetFEModel(), R, m_Fr);

	// zero rigid body reaction forces
	m_rigidSolver.Residual();

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// calculate internal stress force
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// calculate internal thermal work
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEThermoElasticSolidDomain* pdom = dynamic_cast<FEThermoElasticSolidDomain*>(&mesh.Domain(i));
		if (pdom) pdom->InternalThermalWork(R);
	}

    // calculate the body forces
    for (int i=0; i<mesh.Domains(); ++i)
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
	for (int i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->Residual(tp, RHS);
	}

	// calculate contact forces
	ContactForces(RHS);

	// calculate nonlinear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	NonLinearConstraintForces(RHS, tp);

	// add model loads
	int NML = m_fem.ModelLoads();
	for (int i=0; i<NML; ++i)
	{
		FEModelLoad& mli = *m_fem.ModelLoad(i);
		if (mli.IsActive())
		{
			mli.Residual(RHS, tp);
		}
	}

	// set the nodal reaction forces
	// TODO: Is this a good place to do this?
	for (int i=0; i<mesh.Nodes(); ++i)
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

	// total temperature
	vector<double> u = m_Ut + m_Ui + ui;

	// update nodal temperatures
	FEMesh& mesh = m_fem.GetMesh();
	scatter(u, mesh, m_dofT);
}
