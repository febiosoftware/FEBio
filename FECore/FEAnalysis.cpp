#include "stdafx.h"
#include "FEAnalysis.h"
#include "FEModel.h"
#include "FECoreKernel.h"
#include "log.h"
#include "DOFS.h"
#include "MatrixProfile.h"
#include "FEBoundaryCondition.h"
#include "DumpMemStream.h"
#include "FEDataLoadCurve.h"
#include "FELinearConstraintManager.h"
#include "FEShellDomain.h"

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

BEGIN_PARAMETER_LIST(FEAnalysis, FECoreBase)
	ADD_PARAMETER2(m_ntime     , FE_PARAM_INT   , FE_RANGE_GREATER_OR_EQUAL(-1) , "time_steps");
	ADD_PARAMETER2(m_dt0       , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "step_size");
	ADD_PARAMETER2(m_final_time, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "final_time");
END_PARAMETER_LIST()

//-----------------------------------------------------------------------------
FEAnalysis::FEAnalysis(FEModel* pfem) : m_fem(*pfem), FECoreBase(FEANALYSIS_ID)
{
	m_psolver = 0;
	m_tend = 0.0;
	m_nmust = -1;
	m_next_must = -1;

	// --- Analysis data ---
	m_nanalysis = FE_STATIC;	// do quasi-static analysis

	// --- Time Step Data ---
	m_ntime = -1;
	m_final_time = 0.0;
	m_dt = 0;
	m_dt0 = 0;
	m_dtp = 0;
	m_bautostep = false;
	m_iteopt = 11;
	m_dtmax = m_dtmin = 0;
	m_ddt = 0;
	m_nmplc = -1;
	m_naggr = 0;

	// --- Quasi-Newton Solver Variables ---
	m_nretries = 0;
	m_maxretries = 5;

	// initialize counters
	m_ntotref    = 0;		// total nr of stiffness reformations
	m_ntotiter   = 0;		// total nr of non-linear iterations
	m_ntimesteps = 0;		// time steps completed
	m_ntotrhs    = 0;		// total nr of right hand side evaluations

	// --- I/O Data ---
	m_ndump   = FE_DUMP_NEVER;
	m_nplot   = FE_PLOT_MAJOR_ITRS;
	m_nprint  = FE_PRINT_MINOR_ITRS;
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
//! Return a domain
FEDomain* FEAnalysis::Domain(int i)
{
	return &(m_fem.GetMesh().Domain(m_Dom[i])); 
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
void FEAnalysis::Reset()
{
	m_dt = m_dt0;
	m_dtp = m_dt0;
	m_ntotref    = 0;		// total nr of stiffness reformations
	m_ntotiter   = 0;		// total nr of non-linear iterations
	m_ntimesteps = 0;		// time steps completed
	m_ntotrhs    = 0;		// total nr of right hand side evaluations

	// Deactivate the step
	Deactivate();
}

//-----------------------------------------------------------------------------
void FEAnalysis::SetFESolver(FESolver* psolver)
{
	if (m_psolver) delete m_psolver;
	m_psolver = psolver;
}

//-----------------------------------------------------------------------------
//! Data initialization and data chekcing.
bool FEAnalysis::Init()
{
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
	m_tstart = m_fem.GetStartTime();
	m_tend = m_tstart + Dt;

	// For now, add all domains to the analysis step
	FEMesh& mesh = m_fem.GetMesh();
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
		for (int j=0; j<(int)node.m_ID.size(); ++j)	node.m_ID[j] = DOF_INACTIVE;
	}

    // Then, we activate the domains.
    // This will activate the relevant degrees of freedom
    // NOTE: this must be done after the model components are activated.
    // This is to make sure that all initial and prescribed values are applied.
    // Activate all domains
    // Activate shell domains
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dynamic_cast<FEShellDomain*>(&dom))
            dom.Activate();
    }
    // but activate solid domains last (to deal with sandwiched shells)
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (!dynamic_cast<FEShellDomain*>(&dom))
            dom.Activate();
    }
    
	// Now we apply the BC's to the active dofs
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int nbc = (int)node.m_ID.size();
		for (int j=0; j<nbc; ++j)
		{
			if (node.m_ID[j] == DOF_ACTIVE) node.m_ID[j] = node.m_BC[j];
			else node.m_ID[j] = DOF_FIXED;
		}
	}

	// initialize equations
	FESolver* psolver = GetFESolver();
	if (psolver->InitEquations() == false) return false;

	// do one time initialization of solver data
	if (psolver->Init() == false)
	{
		const char* szerr = fecore_get_error_string();
		felog.printf("ERROR: %s\n", szerr);
		felog.printbox("FATAL ERROR","Failed to initialize solver.\nAborting run.\n");
		return false;
	}

	// initialize linear constraints
	// Must be done after equations are initialized
	if (m_fem.GetLinearConstraintManager().Initialize() == false) return false;

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
	// convergence flag
	// we initialize it to true so that when a restart is performed after 
	// the last time step we terminate normally.
	bool bconv = true;

	// calculate end time value
	double starttime = m_fem.GetStartTime();
//	double endtime = m_fem.m_ftime0 + m_ntime*m_dt0;
	double endtime = m_tend;
	const double eps = endtime*1e-7;

	// print initial progress bar
	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		printf("\nProgress:\n");
		for (int i=0; i<50; ++i) printf("\xB0"); printf("\r");
		felog.SetMode(Logfile::LOG_FILE);
	}

	// if we restarted we need to update the timestep
	// before continuing
	if (m_ntimesteps != 0)
	{
		// update time step
		if (m_bautostep && (m_fem.GetCurrentTime() + eps < endtime)) AutoTimeStep(GetFESolver()->m_niter);
	}
	else
	{
		// make sure that the timestep is at least the min time step size
		if (m_bautostep) AutoTimeStep(0);
	}

	// dump stream for running restarts
	DumpMemStream dmp(m_fem);

	// repeat for all timesteps
	m_nretries = 0;
	while (endtime - m_fem.GetCurrentTime() > eps)
	{
		// keep a copy of the current state, in case
		// we need to retry this time step
		if (m_bautostep && (m_maxretries > 0)) 
		{ 
			FESolver* solver = GetFESolver();
			solver->m_UpdateTime.start();
			dmp.clear(); m_fem.Serialize(dmp); 
			solver->m_UpdateTime.stop();
		}

		// update time
		double newTime = m_fem.GetCurrentTime() + m_dt;
		m_fem.SetCurrentTime(newTime);
		felog.printf("\n===== beginning time step %d : %lg =====\n", m_ntimesteps + 1, newTime);

		// initialize the solver step
		// (This basically evaluates all the parameter lists, but let's the solver
		//  customize this process to the specific needs of the solver)
		if (GetFESolver()->InitStep(newTime) == false)
		{
			bconv = false;
			break;
		}

		try
		{
			// do the callback
			m_fem.DoCallback(CB_UPDATE_TIME);

			// solve this timestep,
			bconv = GetFESolver()->SolveStep(newTime);
		}
		catch (LinearSolverFailed)
		{
			bconv = false;
			felog.printbox("FATAL ERROR", "Linear solver failed to find solution. Aborting run.");
			break;
		}
		catch (ExitRequest)
		{
			bconv = false;
			felog.printbox("WARNING", "Early termination on user's request");
			break;
		}
		catch (ZeroDiagonal e)
		{
			bconv = false;
			// TODO: Fix this feature
			felog.printbox("FATAL ERROR", "Zero diagonal detected. Aborting run.");
			break;
		}
		catch (NANDetected)
		{
			bconv = false;
			felog.printbox("ERROR", "NAN Detected.");
		}
		catch (MemException e)
		{
			bconv = false;
			if (e.m_falloc < 1024*1024)
				felog.printbox("FATAL ERROR", "Failed allocating %lg bytes", e.m_falloc);
			else
			{
				double falloc = e.m_falloc / (1024.0*1024.0);
				felog.printbox("FATAL ERROR", "Failed allocating %lg MB", falloc);
			}
			break;
		}
		catch (FEMultiScaleException)
		{
			bconv = false;
			felog.printbox("FATAL ERROR", "The RVE problem has failed. Aborting macro run.");
			break;
		}
		catch (std::bad_alloc e)
		{
			bconv = false;
			felog.printbox("FATAL ERROR", "A memory allocation failure has occured.\nThe program will now be terminated.");
			break;
		}

// We only catch all exceptions for debug versions
#ifndef _DEBUG
		catch (...)
		{
			bconv = false;
			felog.printbox("FATAL ERROR", "An unknown exception has occured.\nThe program will now be terminated.");
			break;
		}
#endif
		// update counters
		FESolver* psolver = GetFESolver();
		m_ntotref  += psolver->m_ntotref;
		m_ntotiter += psolver->m_niter;
		m_ntotrhs  += psolver->m_nrhs;

		// update model's data
		m_fem.UpdateModelData();

		// see if we have converged
		if (bconv)
		{
			// Yes! We have converged!
			if (GetPrintLevel() != FE_PRINT_NEVER)
				felog.printf("\n\n------- converged at time : %lg\n\n", m_fem.GetCurrentTime());

			// update nr of completed timesteps
			m_ntimesteps++;

			// call callback function
			if (m_fem.DoCallback(CB_MAJOR_ITERS) == false)
			{
				bconv = false;
				felog.printbox("WARNING", "Early termination on user's request");
				break;
			}

			// reset retry counter
			m_nretries = 0;

			// update time step
			if (m_bautostep && (m_fem.GetCurrentTime() + eps < endtime)) AutoTimeStep(psolver->m_niter);
		}
		else 
		{
			// This will allow states that have negative Jacobians to be stored
			m_fem.DoCallback(CB_MINOR_ITERS);

			// Report the sad news to the user.
			felog.printf("\n\n------- failed to converge at time : %lg\n\n", m_fem.GetCurrentTime());

			// If we have auto time stepping, decrease time step and let's retry
			if (m_bautostep && (m_nretries < m_maxretries))
			{
				// restore the previous state
				dmp.Open(false, true);
				m_fem.Serialize(dmp);
				
				// let's try again
				Retry();
			}
			else 
			{
				// can't retry, so abort
				if (m_nretries >= m_maxretries)	felog.printf("Max. nr of retries reached.\n\n");

				break;
			}
		}

		// print a progress bar
		if (GetPrintLevel() == FE_PRINT_PROGRESS)
		{
			int l = (int)(50 * m_fem.GetCurrentTime() / endtime);
			for (int i=0; i<l; ++i) printf("\xB2"); printf("\r");
			fflush(stdout);
		}

		// flush the m_log file, so we don't loose anything if 
		// the next timestep goes wrong
		felog.flush();
	}

	// TODO: Why is this here?
	m_fem.SetStartTime(m_fem.GetCurrentTime());

	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		felog.SetMode(Logfile::LOG_FILE_AND_SCREEN);
	}

	if (GetPrintLevel() != FE_PRINT_NEVER)
	{
		// output report
		felog.printf("\n\nN O N L I N E A R   I T E R A T I O N   I N F O R M A T I O N\n\n");
		felog.printf("\tNumber of time steps completed .................... : %d\n\n", m_ntimesteps);
		felog.printf("\tTotal number of equilibrium iterations ............ : %d\n\n", m_ntotiter);
		felog.printf("\tAverage number of equilibrium iterations .......... : %lg\n\n", (double) m_ntotiter / (double) m_ntimesteps);
		felog.printf("\tTotal number of right hand evaluations ............ : %d\n\n", m_ntotrhs);
		felog.printf("\tTotal number of stiffness reformations ............ : %d\n\n", m_ntotref);

		// get and print elapsed time
		char sztime[64];

		GetFESolver()->m_SolverTime.time_str(sztime);
		felog.printf("\tTime in linear solver: %s\n\n", sztime);
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Restores data for a running restart

void FEAnalysis::Retry()
{
	felog.printf("Retrying time step. Retry attempt %d of max %d\n\n", m_nretries+1, m_maxretries);

	// adjust time step
	double dtn;

	if (m_nretries == 0) m_ddt = (m_dt) / (m_maxretries+1);

	if (m_naggr == 0) dtn = m_dt - m_ddt;
	else dtn = m_dt*0.5;

	felog.printf("\nAUTO STEPPER: retry step, dt = %lg\n\n", dtn);

	// increase retry counter
	m_nretries++;

	// the new time step cannot be a must-point
	if (m_nmust != -1)
	{
		// if we were at a must-point, make sure we can hit this must-point again
		m_next_must--;
		m_nmust = -1;
	}

	m_dtp = dtn;
	m_dt = dtn;
}

//-----------------------------------------------------------------------------
//! Adjusts the time step size based on the convergence information.
//!	If the previous time step was able to converge in less than
//! m_fem.m_iteopt iterations the step size is increased, else it
//! is decreased.

void FEAnalysis::AutoTimeStep(int niter)
{
	double dtn = m_dtp;
	double told = m_fem.GetCurrentTime();

	// make sure the timestep size is at least the minimum
	if (dtn < m_dtmin) dtn = m_dtmin;

	// get the max time step
	double dtmax = m_dtmax;
	
	// If we have a must-point load curve
	// we take the max step size from the lc
	if (m_nmplc >= 0)
	{
		FELoadCurve& lc = *m_fem.GetLoadCurve(m_nmplc);
		dtmax = lc.Value(told);
	}

	// adjust time step size
	if (niter > 0)
	{
		double scale = sqrt((double) m_iteopt / (double) niter);

		// Adjust time step size
		if (scale >= 1)
		{	
			dtn = dtn + (dtmax - dtn)*MIN(.20, scale - 1);
			dtn = MIN(dtn, 5.0*m_dtp);
			dtn = MIN(dtn, dtmax);
		}
		else	
		{
			dtn = dtn - (dtn - m_dtmin)*(1 - scale);
			dtn = MAX(dtn, m_dtmin);
			dtn = MIN(dtn, dtmax);
		}

		// Report new time step size
		if (dtn > m_dt)
			felog.printf("\nAUTO STEPPER: increasing time step, dt = %lg\n\n", dtn);
		else if (dtn < m_dt)
			felog.printf("\nAUTO STEPPER: decreasing time step, dt = %lg\n\n", dtn);
	}

	// Store this time step value. This is the value that will be used to evaluate
	// the next time step increment. This will not include adjustments due to the must-point
	// controller since this could create really small time steps that may be difficult to
	// recover from. 
	m_dtp = dtn;

	// check for mustpoints
	if (m_nmplc >= 0) dtn = CheckMustPoints(told, dtn);

	// make sure we are not exceeding the final time
	if (told + dtn > m_tend)
	{
		dtn = m_tend - told;
		felog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
	}

	// store time step size
	m_dt = dtn;
}

//-----------------------------------------------------------------------------
//! This function makes sure that no must points are passed. It returns an
//! updated value (less than dt) if t + dt would pass a must point. Otherwise
//! it returns dt.
//! \param t current time
//! \param dt current time step
//! \return updated time step.
double FEAnalysis::CheckMustPoints(double t, double dt)
{
	double tnew = t + dt;
	double dtnew = dt;
	const double eps = m_tend*1e-12;
	double tmust = tnew + eps;
	FEDataLoadCurve& lc = dynamic_cast<FEDataLoadCurve&>(*m_fem.GetLoadCurve(m_nmplc));
	m_nmust = -1;
	if (m_next_must < lc.Points())
	{
		FEDataLoadCurve::LOADPOINT lp;
		if (m_next_must < 0)
		{
			// find the first must-point that is on or past this time
			m_next_must = 0;
			bool bfound = false;
			do
			{
				lp = lc.LoadPoint(m_next_must);
				if ((tmust > lp.time) && (fabs(tnew - lp.time) > 1e-12)) ++m_next_must;
				else bfound = true;
			}
			while ((bfound == false)&&(m_next_must < lc.Points()));

			// make sure we did not pass all must points
			if (m_next_must >= lc.Points()) return dt;
		}
		else lp = lc.LoadPoint(m_next_must);

		// TODO: what happens when dtnew < dtmin and the next time step fails??
		if (tmust > lp.time)
		{
			dtnew = lp.time - t;
			felog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
			m_nmust = m_next_must++;
		}
		else if (fabs(tnew - lp.time) < 1e-12)
		{
			m_nmust = m_next_must++;
			tnew = lp.time;
			dtnew = tnew - t;
			felog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
		}
		else if (tnew > m_tend)
		{
			dtnew = m_tend - t;
			felog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
			m_nmust = m_next_must++;
		}
	}
	return dtnew;
}

//-----------------------------------------------------------------------------
void FEAnalysis::Serialize(DumpStream& ar)
{
	// don't serialize for shallow copies
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		// --- analysis data ---
		ar << m_nanalysis;
		ar << m_bactive;

		// --- Time Step Data ---
		ar << m_ntime;
		ar << m_final_time;
		ar << m_dt;
		ar << m_dt0;
		ar << m_dtp;
		ar << m_tstart << m_tend;
		ar << m_bautostep;
		ar << m_iteopt;
		ar << m_dtmin;
		ar << m_dtmax;
		ar << m_nmplc;
		ar << m_naggr;

		// --- Quasi-Newton Solver variables ---
		ar << m_nretries;
		ar << m_maxretries;
		ar << m_ntotrhs;
		ar << m_ntotref;
		ar << m_ntotiter;
		ar << m_ntimesteps;

		// --- I/O Data ---
		ar << m_nplot;
		ar << m_nprint;
		ar << m_noutput;
		ar << m_ndump;
		ar << m_nplot_stride;
		ar << m_nplotRange[0] << m_nplotRange[1];
		ar << m_bplotZero;

		// store the class IDs for the active model components
		ar << (int) m_MC.size();
		for (int i=0; i< (int) m_MC.size(); ++i) ar << m_MC[i]->GetClassID();

		// Seriaize solver data
		FESolver* psolver = GetFESolver();
		const char* szsolver = psolver->GetTypeStr();
		ar << psolver->GetTypeStr();
		psolver->Serialize(ar);
	}
	else
	{
		// --- analysis data ---
		ar >> m_nanalysis;
		ar >> m_bactive;

		// --- Time Step Data ---
		ar >> m_ntime;
		ar >> m_final_time;
		ar >> m_dt;
		ar >> m_dt0;
		ar >> m_dtp;
		ar >> m_tstart >> m_tend;
		ar >> m_bautostep;
		ar >> m_iteopt;
		ar >> m_dtmin;
		ar >> m_dtmax;
		ar >> m_nmplc;
		ar >> m_naggr;

		// --- Quasi-Newton Solver variables ---
		ar >> m_nretries;
		ar >> m_maxretries;
		ar >> m_ntotrhs;
		ar >> m_ntotref;
		ar >> m_ntotiter;
		ar >> m_ntimesteps;

		// --- I/O Data ---
		ar >> m_nplot;
		ar >> m_nprint;
		ar >> m_noutput;
		ar >> m_ndump;
		ar >> m_nplot_stride;
		ar >> m_nplotRange[0] >> m_nplotRange[1];
		ar >> m_bplotZero;

#ifdef _DEBUG
		m_ndump = FE_DUMP_NEVER;
#endif

		// read the active model components
		int n, nid;
		ar >> n;
		m_MC.clear();
		for (int i=0; i<n; ++i)
		{
			ar >> nid;
			FEModelComponent* pmc = m_fem.FindModelComponent(nid);
			assert(pmc);
			AddModelComponent(pmc);
		}

		// Serialize solver data
		char szsolver[256] = {0};
		ar >> szsolver;
		assert(m_psolver == 0);
		m_psolver = fecore_new<FESolver>(FESOLVER_ID, szsolver, &m_fem); assert(m_psolver);
		m_psolver->Serialize(ar);

		// For now, add all domains to the analysis step
		FEMesh& mesh = m_fem.GetMesh();
		int ndom = mesh.Domains();
		ClearDomains();
		for (int i = 0; i<ndom; ++i) AddDomain(i);
	}
}
