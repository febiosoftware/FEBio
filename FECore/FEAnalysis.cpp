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

	// --- Analysis data ---
	m_nanalysis = FE_STATIC;	// do quasi-static analysis
	m_istiffpr = 1;				// use pressure stiffness
	m_baugment = false;			// no augmentations

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
	bool bdebug = m_fem.GetDebugFlag();
	m_nplot  = (bdebug?FE_PLOT_MINOR_ITRS     : FE_PLOT_MAJOR_ITRS );
	m_nprint = (bdebug?FE_PRINT_MINOR_ITRS_EXP: FE_PRINT_MINOR_ITRS);
}

//-----------------------------------------------------------------------------
//! Return a domain
FEDomain* FEAnalysis::Domain(int i)
{
	return &(m_fem.GetMesh().Domain(m_Dom[i])); 
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
bool FEAnalysis::Init()
{
	// set first time step
	// We can't do this since it will mess up the value from a restart
//	m_dt = m_dt0;

	// determine the end time
	double Dt;
	if (m_ntime == -1) Dt = m_final_time; else Dt = m_dt0*m_ntime;
	m_tend = m_fem.m_ftime0 + Dt;

	// For now, add all domains to the analysis step
	FEMesh& mesh = m_fem.GetMesh();
	int ndom = mesh.Domains();
	ClearDomains();
	for (int i=0; i<ndom; ++i) AddDomain(i);

	// activate the boundary conditions
	for (int i=0; i<(int) m_BC.size(); ++i) m_BC[i]->Activate();

	// activate contact interface
	for (int i=0; i<(int) m_CI.size(); ++i) m_CI[i]->Activate();

	// activate non-linear constraints
	for (int i=0; i<(int) m_NLC.size(); ++i) m_NLC[i]->Activate();
	
	return true;
}

//-----------------------------------------------------------------------------
//! This function deactivates all boundary conditions and contact interfaces.
//! It also gives the linear solver to clean its data.
//! This is called at the completion of an analysis step.
void FEAnalysis::Finish()
{
	// deactivate the boundary conditions
	for (size_t i=0; i<m_BC.size(); ++i) m_BC[i]->Deactivate();

	// deactivate contact interfaces
	for (size_t i=0; i<m_CI.size(); ++i) m_CI[i]->Deactivate();

	// deactive non-linear constraints
	for (size_t i=0; i<m_NLC.size(); ++i) m_NLC[i]->Deactivate();

	// clean up solver data (i.e. destroy linear solver)
	m_psolver->Clean();
}

//-----------------------------------------------------------------------------
//! This function initializes the linear constraint table (LCT). This table
//! contains for each dof the linear constraint it belongs to. (or -1 if it is
//! not constraint)

bool FEAnalysis::InitConstraints()
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

	// let's do the aug lag linear constraints
	// TODO: This is also done in FEM::Init and FEAnalysis::Init. Where do I really need to do this?
	int N = m_fem.NonlinearConstraints();
	for (i=0; i<N; ++i) m_fem.NonlinearConstraint(i)->Init();

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

	// repeat for all timesteps
	m_nretries = 0;
	while (endtime - m_fem.m_ftime > eps)
	{
		// keep a copy of the current state, in case
		// we need to retry this time step
		if (m_bautostep) { dmp.clear(); m_fem.ShallowCopy(dmp, true); }

		// update time
		m_fem.m_ftime += m_dt;

		int i;

		// evaluate load curve values at current time
		for (i=0; i<m_fem.LoadCurves(); ++i) m_fem.GetLoadCurve(i)->Evaluate(m_fem.m_ftime);

		// evaluate material parameter lists
		for (i=0; i<m_fem.Materials(); ++i)
		{
			// get the material
			FEMaterial* pm = m_fem.GetMaterial(i);

			// evaluate its parameter list
			m_fem.EvaluateMaterialParameters(pm);
		}

		// initialize materials
		// TODO: I need to do this since the material parameters can have changed and thus a new initialization
		//       needs to be done to see if the material parameters are still valid. I would like to add value checking
		//       directly in the parameter evaluation above so this can be removed.
		if (m_fem.InitMaterials() == false)
		{
			bconv = false;
			break;
		}

		// evaluate body-force parameter lists
		for (i=0; i<m_fem.BodyLoads(); ++i)
		{
			FEParameterList& pl = m_fem.GetBodyLoad(i)->GetParameterList();
			m_fem.EvaluateParameterList(pl);
		}

		// evaluate contact interface parameter lists
		for (i=0; i<m_fem.SurfacePairInteractions(); ++i)
		{
			FEParameterList& pl = m_fem.SurfacePairInteraction(i)->GetParameterList();
			m_fem.EvaluateParameterList(pl);
		}

		// evaluate constraint parameter lists
		for (i=0; i<m_fem.NonlinearConstraints(); ++i)
		{
			FEParameterList& pl = m_fem.NonlinearConstraint(i)->GetParameterList();
			m_fem.EvaluateParameterList(pl);
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
			felog.printbox("FATAL ERROR", "NAN Detected. Run aborted.");
			break;
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
			felog.printf("\n\n------- converged at time : %lg\n\n", m_fem.m_ftime);

			// update nr of completed timesteps
			m_ntimesteps++;

			// output results to plot database
			if ((m_nplot != FE_PLOT_NEVER) && (m_nplot != FE_PLOT_FINAL))
			{
				if ((m_nplot == FE_PLOT_MUST_POINTS) && (m_nmplc >= 0))
				{
					FELoadCurve& lc = *m_fem.GetLoadCurve(m_nmplc);
					if (lc.HasPoint(m_fem.m_ftime)) m_fem.Write();
				}
				else m_fem.Write();
			}

			// Dump converged state to the archive
			if (m_bDump) m_fem.DumpData();

			// store additional data to the logfile
			m_fem.WriteData();

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

	if ((m_nplot == FE_PLOT_FINAL) && bconv) m_fem.Write();

	m_fem.m_ftime0 = m_fem.m_ftime;

	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		felog.SetMode(Logfile::FILE_AND_SCREEN);
	}

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
	if (m_nmplc >= 0)
	{
		double tnew = told + dtn;
		const double eps = m_tend*1e-07;
		FELoadCurve& lc = *m_fem.GetLoadCurve(m_nmplc);
		int n = lc.FindPoint(told+eps);
		if (n >= 0)
		{
			// TODO: what happens when dtn < dtmin and the next time step fails??
			double tmust = lc.LoadPoint(n).time;
			if (tnew > tmust)
			{
				dtn = tmust - told;
				felog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
			}
			else if (tnew > m_tend)
			{
				dtn = m_tend - told;
				felog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
			}
		}
	}

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
void FEAnalysis::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		// --- analysis data ---
		ar << m_ntype;
		ar << m_nanalysis;
		ar << m_istiffpr;
		ar << m_baugment;

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
		ar << m_bDump;

		// boundary conditions
		ar << (int) m_BC.size();
		for (int i=0; i< (int) m_BC.size(); ++i) ar << m_BC[i]->GetID();

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
		ar >> m_baugment;

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
		ar >> m_bDump;
#ifdef _DEBUG
		m_bDump = false;
#endif

		// boundary conditions
		int n, nbc;
		ar >> n;
		m_BC.clear();
		for (int i=0; i<n; ++i)
		{
			ar >> nbc;
			FEBoundaryCondition* pbc = m_fem.FindBC(nbc);
			assert(pbc);
			m_BC.push_back(pbc);
		}

		// Serialize solver data
		char szsolver[256] = {0};
		ar >> szsolver;
		assert(m_psolver == 0);
		m_psolver = fecore_new<FESolver>(FESOLVER_ID, szsolver, &m_fem); assert(m_psolver);
		m_psolver->Serialize(ar);
	}
}
