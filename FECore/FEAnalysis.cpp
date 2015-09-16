#include "stdafx.h"
#include "FEAnalysis.h"
#include "FEModel.h"
#include "FECoreKernel.h"
#include "log.h"
#include "DOFS.h"

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

//-----------------------------------------------------------------------------
FEAnalysis::FEAnalysis(FEModel* pfem, int ntype) : FECoreBase(FEANALYSIS_ID), m_fem(*pfem), m_ntype(ntype) 
{
	m_psolver = 0;
	m_tend = 0.0;
	m_nmust = -1;

	// --- Analysis data ---
	m_nanalysis = FE_STATIC;	// do quasi-static analysis
	m_istiffpr = 1;				// use pressure stiffness

	// --- Time Step Data ---
	m_ntime = -1;
	m_final_time = 0.0;
	m_dt = 0;
	m_dt0 = 0;
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
	m_bDump = false;
	m_nplot   = FE_PLOT_MAJOR_ITRS;
	m_nprint  = FE_PRINT_MINOR_ITRS;
	m_noutput = FE_OUTPUT_MAJOR_ITRS;
	m_nplot_stride = 1;
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
void FEAnalysis::Reset()
{
	m_dt = m_dt0;
	m_ntotref    = 0;		// total nr of stiffness reformations
	m_ntotiter   = 0;		// total nr of non-linear iterations
	m_ntimesteps = 0;		// time steps completed
	m_ntotrhs    = 0;		// total nr of right hand side evaluations
}

//-----------------------------------------------------------------------------
//! Data initialization and data chekcing.
bool FEAnalysis::Init()
{
	if ((m_ntime <= 0) && (m_final_time <= 0.0)) { felog.printf("Invalid number of time steps for analysis step.\n"); return false; }
	if ((m_ntime >  0) && (m_final_time >  0.0)) { felog.printf("You must either set the number of time steps or the final time but not both.\n"); return false; }
	if (m_dt0   <= 0) { felog.printf("Invalid time step size for analysis step\n"); return false; }
	if (m_bautostep)
	{
//		if (m_pStep->m_dtmin <= 0) return err("Invalid minimum time step size");
//		if (m_pStep->m_dtmax <= 0) return err("Invalid maximum time step size");
	}
	return true;
}

//-----------------------------------------------------------------------------
//! This function gets called right before the step needs to be solved.
bool FEAnalysis::Activate()
{
	// set first time step
	// We can't do this since it will mess up the value from a restart
//	m_dt = m_dt0;

	// determine the end time
	double Dt;
	if (m_ntime == -1) Dt = m_final_time; else Dt = m_dt0*m_ntime;
	m_tstart = m_fem.m_ftime0;
	m_tend = m_fem.m_ftime0 + Dt;

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
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.Activate();
	}

	// Now we apply the BC's to the active dofs
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int nbc = node.m_ID.size();
		for (int j=0; j<nbc; ++j)
		{
			if (node.m_ID[j] == DOF_ACTIVE) node.m_ID[j] = node.m_BC[j];
			else node.m_ID[j] = DOF_FIXED;
		}
	}

	// initialize equations
	if (m_psolver->InitEquations() == false) return false;

	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		felog.printbox("FATAL ERROR","Failed to initialize solver.\nAborting run.\n");
		return false;
	}

	// initialize linear constraints
	// Must be done after equations are initialized
	if (InitLinearConstraints() == false) return false;

	// activate the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il) il->Activate();
	}

	// write the initial state
	m_fem.Write(FE_STEP_INITIALIZED);
	
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
	m_psolver->Clean();
}

//-----------------------------------------------------------------------------
//! This function initializes the linear constraint table (LCT). This table
//! contains for each dof the linear constraint it belongs to. (or -1 if it is
//! not constraint)

bool FEAnalysis::InitLinearConstraints()
{
	int nlin = m_fem.m_LinC.size();
	if (nlin == 0) return true;
	int i;

	FEMesh& mesh = m_fem.GetMesh();

	// set the equation numbers for the linear constraints
	list<FELinearConstraint>::iterator it = m_fem.m_LinC.begin();
	for (i=0; i<nlin; ++i, ++it)
	{
		FELinearConstraint& lc = *it;
		lc.master.neq = mesh.Node(lc.master.node).m_ID[lc.master.bc];

		// make sure the master did not get assigned an equation
		assert(lc.master.neq == -1);

		// set the slave equation numbers
		list<FELinearConstraint::SlaveDOF>::iterator is = lc.slave.begin();
		int nn = lc.slave.size();
		for (int n=0; n<nn; ++n, ++is)
		{
			FELinearConstraint::SlaveDOF& sn = *is;
			sn.neq = mesh.Node(sn.node).m_ID[sn.bc];
		}		
	}

	// create the linear constraint table
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_NDOFS = fedofs.GetNDOFS();
	m_fem.m_LCT.assign(mesh.Nodes()*MAX_NDOFS, -1);

	list<FELinearConstraint>::iterator ic = m_fem.m_LinC.begin();
	for (i=0; i<nlin; ++i, ++ic)
	{
		FELinearConstraint& lc = *ic;
		int n = lc.master.node;
		int m = lc.master.bc;
		
		m_fem.m_LCT[n*MAX_NDOFS+m] = i;
	}

	// to simplify accessing the linear constraint data
	// we store all pointers in an array
	// TODO: perhaps I should store the linear constraints that way
	// anyways and get rid of the list
	m_fem.m_LCA.resize(nlin);
	ic = m_fem.m_LinC.begin();
	for (i=0; i<nlin; ++i, ++ic) m_fem.m_LCA[i] = &(*ic);

	return true;
}

//-----------------------------------------------------------------------------
bool FEAnalysis::Solve()
{
	// convergence flag
	// we initialize it to true so that when a restart is performed after 
	// the last time step we terminate normally.
	bool bconv = true;

	// calculate end time value
	double starttime = m_fem.m_ftime0;
//	double endtime = m_fem.m_ftime0 + m_ntime*m_dt0;
	double endtime = m_tend;
	const double eps = endtime*1e-7;

	// print initial progress bar
	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		printf("\nProgress:\n");
		for (int i=0; i<50; ++i) printf("\xB0"); printf("\r");
		felog.SetMode(Logfile::FILE_ONLY);
	}

	// if we restarted we need to update the timestep
	// before continuing
	if (m_ntimesteps != 0)
	{
		// update time step
		if (m_bautostep && (m_fem.m_ftime + eps < endtime)) AutoTimeStep(m_psolver->m_niter);
	}
	else
	{
		// make sure that the timestep is at least the min time step size
		if (m_bautostep) AutoTimeStep(0);
	}

	// dump stream for running restarts
	DumpStream dmp;

	// Flags used to note whether the last converged time step was stored.
	// This is used by to make sure that the final time step is stored to the plot file and data files.
	bool boutput = false;

	// repeat for all timesteps
	m_nretries = 0;
	while (endtime - m_fem.m_ftime > eps)
	{
		// keep a copy of the current state, in case
		// we need to retry this time step
		if (m_bautostep) { dmp.clear(); m_fem.ShallowCopy(dmp, true); }

		// update time
		m_fem.m_ftime += m_dt;

		// initialize the solver step
		// (This basically evaluates all the parameter lists, but let's the solver
		//  customize this process to the specific needs of the solver)
		if (m_psolver->InitStep(m_fem.m_ftime) == false)
		{
			bconv = false;
			break;
		}

		// solve this timestep,
		try
		{
			bconv = m_psolver->SolveStep(m_fem.m_ftime);
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
		m_ntotref  += m_psolver->m_ntotref;
		m_ntotiter += m_psolver->m_niter;
		m_ntotrhs  += m_psolver->m_nrhs;

		// see if we have converged
		if (bconv)
		{
			// Yes! We have converged!
			if (GetPrintLevel() != FE_PRINT_NEVER)
				felog.printf("\n\n------- converged at time : %lg\n\n", m_fem.m_ftime);

			// update nr of completed timesteps
			m_ntimesteps++;

			// output results to plot database
			m_fem.Write(FE_CONVERGED);

			// Dump converged state to the archive
			if (m_bDump) m_fem.DumpData();

			// store additional data to the logfile
			boutput = false;
			if ((m_noutput != FE_OUTPUT_NEVER) && (m_noutput != FE_OUTPUT_FINAL))
			{
				if (m_noutput == FE_OUTPUT_MUST_POINTS)
				{
					if (m_nmust >= 0) { m_fem.WriteData(); boutput = true; }
				}
				else
				{
					m_fem.WriteData(); boutput = true;
				}
			}

			// update time step
			if (m_bautostep && (m_fem.m_ftime + eps < endtime)) AutoTimeStep(m_psolver->m_niter);

			// reset retry counter
			m_nretries = 0;

			// call callback function
			m_fem.DoCallback(CB_MAJOR_ITERS);
		}
		else 
		{
			// Report the sad news to the user.
			felog.printf("\n\n------- failed to converge at time : %lg\n\n", m_fem.m_ftime);

			// plot the state when the debug flag is on
			m_fem.Write(FE_UNCONVERGED);

			// If we have auto time stepping, decrease time step and let's retry
			if (m_bautostep && (m_nretries < m_maxretries))
			{
				// restore the previous state
				dmp.set_position(0);
				m_fem.ShallowCopy(dmp, false);
				
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
			int l = (int)(50*m_fem.m_ftime / endtime);
			for (int i=0; i<l; ++i) printf("\xB2"); printf("\r");
			fflush(stdout);
		}

		// flush the m_log file, so we don't loose anything if 
		// the next timestep goes wrong
		felog.flush();
	}

	// write the final time step, if it hasn't been written yet.
	if (bconv) m_fem.Write(FE_STEP_SOLVED);

	// output the final time step, if it hasn't been written yet.
	if (m_noutput != FE_OUTPUT_NEVER)
	{
		if (bconv && (boutput == false)) m_fem.WriteData();
	}

	m_fem.m_ftime0 = m_fem.m_ftime;

	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		felog.SetMode(Logfile::FILE_AND_SCREEN);
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

		m_psolver->m_SolverTime.time_str(sztime);
		felog.printf("\tTime in solver: %s\n\n", sztime);
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
	m_nmust = -1;

	m_dt = dtn;
}

//-----------------------------------------------------------------------------
//! Adjusts the time step size based on the convergence information.
//!	If the previous time step was able to converge in less than
//! m_fem.m_iteopt iterations the step size is increased, else it
//! is decreased.

void FEAnalysis::AutoTimeStep(int niter)
{
	double dtn = m_dt;
	double told = m_fem.m_ftime;

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
			dtn = MIN(dtn, 5.0*m_dt);
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
	const double eps = m_tend*1e-07;
	FELoadCurve& lc = *m_fem.GetLoadCurve(m_nmplc);
	double tmust = tnew;
	m_nmust = -1;
	int n = lc.FindPoint(t+eps, tmust);
	if (n >= 0)
	{
		// TODO: what happens when dtnew < dtmin and the next time step fails??
		if (tnew > tmust)
		{
			dtnew = tmust - t;
			felog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
			m_nmust = n;
		}
		else if (tnew == tmust) { m_nmust = n; }
		else if (tnew > m_tend)
		{
			dtnew = m_tend - t;
			felog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
			m_nmust = n;
		}
	}
	return dtnew;
}

//-----------------------------------------------------------------------------
void FEAnalysis::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		// --- analysis data ---
		ar << m_ntype;
		ar << m_nanalysis;
		ar << m_istiffpr;

		// --- Time Step Data ---
		ar << m_ntime;
		ar << m_final_time;
		ar << m_dt;
		ar << m_dt0;
		ar << m_tend;
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
		ar << m_bDump;
		ar << m_nplot_stride;

		// store the class IDs for the active model components
		ar << (int) m_MC.size();
		for (int i=0; i< (int) m_MC.size(); ++i) ar << m_MC[i]->GetClassID();

		// Seriaize solver data
		ar << m_psolver->GetTypeStr();
		m_psolver->Serialize(ar);
	}
	else
	{
		// --- analysis data ---
		ar >> m_ntype;
		ar >> m_nanalysis;
		ar >> m_istiffpr;

		// --- Time Step Data ---
		ar >> m_ntime;
		ar >> m_final_time;
		ar >> m_dt;
		ar >> m_dt0;
		ar >> m_tend;
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
		ar >> m_bDump;
		ar >> m_nplot_stride;

#ifdef _DEBUG
		m_bDump = false;
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
	}
}
