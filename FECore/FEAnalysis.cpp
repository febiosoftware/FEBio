#include "stdafx.h"
#include "FEAnalysis.h"
#include "FEModel.h"
#include "FECoreKernel.h"
#include "log.h"
#include "DOFS.h"
#include "MatrixProfile.h"
#include "FEBoundaryCondition.h"
#include "DumpMemStream.h"
#include "FELinearConstraintManager.h"
#include "FEShellDomain.h"
#include "fecore_error.h"

REGISTER_SUPER_CLASS(FEAnalysis, FEANALYSIS_ID);

BEGIN_FECORE_CLASS(FEAnalysis, FECoreBase)
	ADD_PARAMETER(m_ntime       , FE_RANGE_GREATER_OR_EQUAL(-1) , "time_steps");
	ADD_PARAMETER(m_dt0         , FE_RANGE_GREATER_OR_EQUAL(0.0), "step_size");
	ADD_PARAMETER(m_final_time  , FE_RANGE_GREATER_OR_EQUAL(0.0), "final_time");
	ADD_PARAMETER(m_bplotZero   , "plot_zero_state");
	ADD_PARAMETER(m_nplotRange  , 2, "plot_range");
	ADD_PARAMETER(m_nplot       , "plot_level", 0, "PLOT_NEVER\0PLOT_MAJOR_ITRS\0PLOT_MINOR_ITRS\0PLOT_MUST_POINTS\0PLOT_FINAL\0PLOT_AUGMENTATIONS\0PLOT_STEP_FINAL\0");
	ADD_PARAMETER(m_noutput     , "output_level", 0, "OUTPUT_NEVER\0OUTPUT_MAJOR_ITRS\0OUTPUT_MINOR_ITRS\0OUTPUT_MUST_POINTS\0OUTPUT_FINAL\0");
	ADD_PARAMETER(m_nplot_stride, "plot_stride");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
FEAnalysis::FEAnalysis(FEModel* fem) : FECoreBase(fem), m_timeController(this)
{
	m_psolver = nullptr;
	m_tend = 0.0;

	// --- Analysis data ---
	m_nanalysis = FE_STATIC;	// do quasi-static analysis

	// --- Time Step Data ---
	m_ntime = -1;
	m_final_time = 0.0;
	m_dt0 = 0;
	m_dt = 0;
	m_bautostep = false;

	// initialize counters
	m_ntotref    = 0;		// total nr of stiffness reformations
	m_ntotiter   = 0;		// total nr of non-linear iterations
	m_ntimesteps = 0;		// time steps completed
	m_ntotrhs    = 0;		// total nr of right hand side evaluations

	// --- I/O Data ---
	m_nplot   = FE_PLOT_MAJOR_ITRS;
	m_noutput = FE_OUTPUT_MAJOR_ITRS;
	m_nplot_stride = 1;
	m_nplotRange[0] = 0; // by default, will store step zero.
	m_nplotRange[1] = -1; // by default, store last time step
	m_bplotZero = false; // don't force plotting step zero.

	m_bactive = false;
}

//-----------------------------------------------------------------------------
FEAnalysis::~FEAnalysis()
{
	if (m_psolver) delete m_psolver;
}

//-----------------------------------------------------------------------------
//! copy data from another analysis
void FEAnalysis::CopyFrom(FEAnalysis* step)
{
	m_nanalysis = step->m_nanalysis;

	m_ntime      = step->m_ntime;
	m_final_time = step->m_final_time;
	m_dt         = step->m_dt;
	m_dt0        = step->m_dt0;
	m_tstart     = step->m_tstart;
	m_tend       = step->m_tend;
	m_bautostep  = step->m_bautostep;

	m_timeController.CopyFrom(&step->m_timeController);
}

//-----------------------------------------------------------------------------
//! Return a domain
FEDomain* FEAnalysis::Domain(int i)
{
	return &(GetFEModel()->GetMesh().Domain(m_Dom[i])); 
}

//-----------------------------------------------------------------------------
void FEAnalysis::AddModelComponent(FEModelComponent* pmc)
{
	if (pmc) m_MC.push_back(pmc);
}

//-----------------------------------------------------------------------------
int FEAnalysis::ModelComponents() const
{
	return (int) m_MC.size();
}

//-----------------------------------------------------------------------------
//! sets the plot level
void FEAnalysis::SetPlotLevel(int n) { m_nplot = n; }

//-----------------------------------------------------------------------------
//! sets the plot stride
void FEAnalysis::SetPlotStride(int n) { m_nplot_stride = n; }

//-----------------------------------------------------------------------------
//! sets the plot range
void FEAnalysis::SetPlotRange(int n0, int n1)
{
	m_nplotRange[0] = n0;
	m_nplotRange[1] = n1; 
}

//-----------------------------------------------------------------------------
//! sets the zero-state plot flag
void FEAnalysis::SetPlotZeroState(bool b)
{
	m_bplotZero = b;	
}

//-----------------------------------------------------------------------------
//! get the plot level
int FEAnalysis::GetPlotLevel() { return m_nplot; }

//! Set the output level
void FEAnalysis::SetOutputLevel(int n) { m_noutput = n; }

//! Get the output level
int FEAnalysis::GetOutputLevel() { return m_noutput; }

//-----------------------------------------------------------------------------
void FEAnalysis::Reset()
{
	m_ntotref    = 0;		// total nr of stiffness reformations
	m_ntotiter   = 0;		// total nr of non-linear iterations
	m_ntimesteps = 0;		// time steps completed
	m_ntotrhs    = 0;		// total nr of right hand side evaluations

	m_dt = m_dt0;

	m_timeController.Reset();

	// Deactivate the step
	Deactivate();
}

//-----------------------------------------------------------------------------
FESolver* FEAnalysis::GetFESolver() 
{ 
	return m_psolver; 
}

//-----------------------------------------------------------------------------
void FEAnalysis::SetFESolver(FESolver* psolver)
{
	if (m_psolver) delete m_psolver;
	m_psolver = psolver;
}

//-----------------------------------------------------------------------------
//! Data initialization and data checking.
bool FEAnalysis::Init()
{
	m_dt = m_dt0;

	if (m_timeController.Init() == false) return false;
	if (m_nplot_stride <= 0) return false;
	return Validate();
}

//-----------------------------------------------------------------------------
//! See if this step is active
bool FEAnalysis::IsActive()
{
	return m_bactive;
}

//-----------------------------------------------------------------------------
//! This function gets called right before the step needs to be solved.
bool FEAnalysis::Activate()
{
	FEModel& fem = *GetFEModel();

	// Make sure we are not activated yet
	// This can happen after a restart during FEModel::Solve
	if (m_bactive) return true;

	// activate the time step
	m_bactive = true;

	// set first time step
	// We can't do this since it will mess up the value from a restart
//	m_dt = m_dt0;

	// determine the end time
	double Dt;
	if (m_ntime == -1) Dt = m_final_time; else Dt = m_dt0*m_ntime;
	m_tstart = fem.GetStartTime();
	m_tend = m_tstart + Dt;

	// For now, add all domains to the analysis step
	FEMesh& mesh = fem.GetMesh();
	int ndom = mesh.Domains();
	ClearDomains();
	for (int i=0; i<ndom; ++i) AddDomain(i);

	// activate the model components assigned to this step
	// NOTE: This currently does not ensure that initial conditions are
	// applied first. This is important since relative prescribed displacements must 
	// be applied after initial conditions.
	for (int i=0; i<(int) m_MC.size(); ++i) m_MC[i]->Activate();

	// Next, we need to determine which degrees of freedom are active. 
	// We start by resetting all nodal degrees of freedom.
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<(int)node.dofs(); ++j) node.set_inactive(j);
	}

    // Then, we activate the domains.
    // This will activate the relevant degrees of freedom
    // NOTE: this must be done after the model components are activated.
    // This is to make sure that all initial and prescribed values are applied.
    // Activate all domains
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.Class() != FE_DOMAIN_SHELL)
            dom.Activate();
    }
    // but activate shell domains last (to deal with sandwiched shells)
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.Class() == FE_DOMAIN_SHELL)
            dom.Activate();
    }

	// initialize equations
	FESolver* psolver = GetFESolver();
	if (psolver->InitEquations() == false) return false;

	// do one time initialization of solver data
	if (psolver->Init() == false) return false;

	// initialize linear constraints
	// Must be done after equations are initialized
	if (fem.GetLinearConstraintManager().Initialize() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! This function deactivates all boundary conditions and contact interfaces.
//! It also gives the linear solver to clean its data.
//! This is called at the completion of an analysis step.
void FEAnalysis::Deactivate()
{
	// deactivate the model components
	for (size_t i=0; i<(int) m_MC.size(); ++i) m_MC[i]->Deactivate();

	// clean up solver data (i.e. destroy linear solver)
	GetFESolver()->Clean();

	// deactivate the time step
	m_bactive = false;
}

//-----------------------------------------------------------------------------
bool FEAnalysis::Solve()
{
	FEModel& fem = *GetFEModel();

	// convergence flag
	// we initialize it to true so that when a restart is performed after 
	// the last time step we terminate normally.
	bool bconv = true;

	// calculate end time value
	double starttime = fem.GetStartTime();
//	double endtime = fem.m_ftime0 + m_ntime*m_dt0;
	double endtime = m_tend;
	const double eps = endtime*1e-7;

	// if we restarted we need to update the timestep
	// before continuing
	if (m_ntimesteps != 0)
	{
		// update time step
		if (m_bautostep && (fem.GetCurrentTime() + eps < endtime)) m_timeController.AutoTimeStep(GetFESolver()->m_niter);
	}
	else
	{
		// make sure that the timestep is at least the min time step size
		if (m_bautostep) m_timeController.AutoTimeStep(0);
	}

	// dump stream for running restarts
	DumpMemStream dmp(fem);

	// repeat for all timesteps
	m_timeController.m_nretries = 0;
	while (endtime - fem.GetCurrentTime() > eps)
	{
		// keep a copy of the current state, in case
		// we need to retry this time step
		if (m_bautostep && (m_timeController.m_maxretries > 0))
		{ 
			dmp.clear();
			fem.Serialize(dmp); 
		}

		// Inform that the time is about to change. (Plugins can use 
		// this callback to modify time step)
		fem.DoCallback(CB_UPDATE_TIME);

		// update time
		FETimeInfo& tp = fem.GetTime();
		double newTime = tp.currentTime + m_dt;
		tp.currentTime = newTime;
		tp.timeIncrement = m_dt;
		feLog("\n===== beginning time step %d : %lg =====\n", m_ntimesteps + 1, newTime);

		// initialize the solver step
		// (This basically evaluates all the parameter lists, but let's the solver
		//  customize this process to the specific needs of the solver)
		if (GetFESolver()->InitStep(newTime) == false)
		{
			bconv = false;
			break;
		}

		// call the FE Solver
		int ierr = CallFESolver();

		// see if we want to abort
		if (ierr == 2) 
		{
			bconv = false;
			break;
		}

		// update counters
		FESolver* psolver = GetFESolver();
		m_ntotref  += psolver->m_ntotref;
		m_ntotiter += psolver->m_niter;
		m_ntotrhs  += psolver->m_nrhs;

		// update model's data
		fem.UpdateModelData();

		// see if we have converged
		if (ierr == 0)
		{
			bconv = true;

			// Yes! We have converged!
			feLog("\n------- converged at time : %lg\n\n", fem.GetCurrentTime());

			// update nr of completed timesteps
			m_ntimesteps++;

			// call callback function
			if (fem.DoCallback(CB_MAJOR_ITERS) == false)
			{
				bconv = false;
				feLogWarning("Early termination on user's request");
				break;
			}

			// reset retry counter
			m_timeController.m_nretries = 0;

			// update time step
			if (m_bautostep && (fem.GetCurrentTime() + eps < endtime)) m_timeController.AutoTimeStep(psolver->m_niter);
		}
		else 
		{
			// We failed to converge. 
			bconv = false;

			// This will allow states that have negative Jacobians to be stored
			fem.DoCallback(CB_MINOR_ITERS);

			// Report the sad news to the user.
			feLog("\n\n------- failed to converge at time : %lg\n\n", fem.GetCurrentTime());

			// If we have auto time stepping, decrease time step and let's retry
			if (m_bautostep && (m_timeController.m_nretries < m_timeController.m_maxretries))
			{
				// restore the previous state
				dmp.Open(false, true);
				fem.Serialize(dmp);
				
				// let's try again
				m_timeController.Retry();

				// rewind the solver
				GetFESolver()->Rewind();
			}
			else 
			{
				// can't retry, so abort
				if (m_timeController.m_nretries >= m_timeController.m_maxretries)
					feLog("Max. nr of retries reached.\n\n");

				break;
			}
		}
	}

	// TODO: Why is this here?
	fem.SetStartTime(fem.GetCurrentTime());

	return bconv;
}

//-----------------------------------------------------------------------------
// This function calls the FE Solver for solving this analysis and also handles
// all the exceptions. 
int FEAnalysis::CallFESolver()
{
	int nerr = 0;
	try
	{

		// solve this timestep,
		bool bconv = GetFESolver()->SolveStep();
		nerr = (bconv ? 0 : 1);
	}
	catch (LinearSolverFailed)
	{
		feLogError("FATAL ERROR", "Linear solver failed to find solution. Aborting run.");
		nerr = 2;
	}
	catch (ZeroDiagonal e)
	{
		// TODO: Fix this feature
		feLogError("FATAL ERROR", "Zero diagonal detected. Aborting run.");
		nerr = 2;
	}
	catch (NANDetected)
	{
		feLogError("ERROR", "NAN Detected.");
		nerr = 1;	// don't abort, instead let's retry the step
	}
	catch (FEMultiScaleException)
	{
		feLogError("FATAL ERROR", "The RVE problem has failed. Aborting macro run.");
		nerr = 2;
	}
	catch (std::bad_alloc e)
	{
		feLogError("FATAL ERROR", "A memory allocation failure has occured.\nThe program will now be terminated.");
		nerr = 2;
	}

	return nerr;
}

//-----------------------------------------------------------------------------
void FEAnalysis::Serialize(DumpStream& ar)
{
	FEModel& fem = *GetFEModel();

	// don't serialize for shallow copies
	if (ar.IsShallow()) return;

	// --- analysis data ---
	ar & m_nanalysis;
	ar & m_bactive;

	// --- Time Step Data ---
	ar & m_ntime;
	ar & m_final_time;
	ar & m_dt0 & m_dt;
	ar & m_tstart & m_tend;
	ar & m_bautostep;
	ar & m_ntotrhs;
	ar & m_ntotref;
	ar & m_ntotiter;
	ar & m_ntimesteps;

	// --- I/O Data ---
	ar & m_nplot;
	ar & m_noutput;
	ar & m_nplot_stride;
	ar & m_nplotRange;
	ar & m_bplotZero;

	// Serialize solver data
	ar & m_psolver;

	if (ar.IsSaving())
	{
		// store the class IDs for the active model components
		ar << (int) m_MC.size();
		for (int i=0; i< (int) m_MC.size(); ++i) ar << m_MC[i]->GetClassID();
	}
	else
	{
		// read the active model components
		int n, nid;
		ar >> n;
		m_MC.clear();
		for (int i=0; i<n; ++i)
		{
			ar >> nid;
			FEModelComponent* pmc = fem.FindModelComponent(nid);
			assert(pmc);
			AddModelComponent(pmc);
		}

		// For now, add all domains to the analysis step
		FEMesh& mesh = fem.GetMesh();
		int ndom = mesh.Domains();
		ClearDomains();
		for (int i = 0; i<ndom; ++i) AddDomain(i);
	}

	// serialize time controller
	m_timeController.Serialize(ar);
}
