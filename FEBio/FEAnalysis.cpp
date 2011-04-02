#include "stdafx.h"
#include "FEAnalysis.h"
#include "fem.h"
#include "console.h"
#include "FERigid.h"
#include "FEUncoupledMaterial.h"
#include "log.h"
#include "FESolidSolver.h"
#include "FEHeatSolver.h"
#include "FEPoroElastic.h"
#include "FEPoroSolidSolver.h"
#include "FEPoroSoluteSolver.h"

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

//-----------------------------------------------------------------------------
//! constructor
FEAnalysis::FEAnalysis(FEM& fem) : m_fem(fem)
{
	// --- Analysis data ---
	m_nModule = FE_SOLID;		// solid-mechanics problem
	m_nanalysis = FE_STATIC;	// do quasi-static analysis
	m_istiffpr = 1;				// use pressure stiffness
	m_baugment = false;			// no augmentations
	m_hg = 1.0;

	// --- Time Step Data ---
	m_ntime = 0;
	m_dt = 0;
	m_dt0 = 0;
	m_bautostep = false;
	m_iteopt = 11;
	m_dtmax = m_dtmin = 0;
	m_nmplc = -1;
	m_naggr = 0;

	// --- Quasi-Newton Solver Variables ---
	m_psolver = 0;
	m_nretries = 0;
	m_maxretries = 5;

	// initialize counters
	m_ntotref    = 0;		// total nr of stiffness reformations
	m_ntotiter   = 0;		// total nr of non-linear iterations
	m_ntimesteps = 0;		// time steps completed
	m_ntotrhs    = 0;		// total nr of right hand side evaluations

	// --- I/O Data ---
	m_bDump = false;
	bool bdebug = fem.GetDebugFlag();
	m_nplot  = (bdebug?FE_PLOT_MINOR_ITRS     : FE_PLOT_MAJOR_ITRS );
	m_nprint = (bdebug?FE_PRINT_MINOR_ITRS_EXP: FE_PRINT_MINOR_ITRS);
}

//-----------------------------------------------------------------------------
FEAnalysis::~FEAnalysis(void)
{
}

//-----------------------------------------------------------------------------

void FEAnalysis::Finish()
{
	// deactivate the boundary conditions
	for (size_t i=0; i<m_BC.size(); ++i) m_BC[i]->Deactivate();

	// clean up solver data (i.e. destroy linear solver)
	m_psolver->Clean();
}

//-----------------------------------------------------------------------------
bool FEAnalysis::Init()
{
	int i, j, n;

	// set first time step
	// We can't do this since it will mess up the value from a restart
//	m_dt = m_dt0;

	m_tend = m_fem.m_ftime0 + m_dt0*m_ntime;

	// init must point curve
	if (m_nmplc < 0)
	{
		FELoadCurve* plc = new FELoadCurve();
		plc->Create(2);
		plc->LoadPoint(0).time  = m_fem.m_ftime0;
		plc->LoadPoint(0).value = 0;
		plc->LoadPoint(1).time  = m_fem.m_ftime0 + m_dt*m_ntime;
		plc->LoadPoint(1).value = m_dtmax;
		m_fem.AddLoadCurve(plc);
		m_nmplc = m_fem.m_LC.size()-1;
	}

	// the must point load curve must be evaluated
	// using a step interpolation
	m_fem.m_LC[m_nmplc]->SetInterpolation(FELoadCurve::STEP);

	// activate the boundary conditions
	for (i=0; i<(int) m_BC.size(); ++i) m_BC[i]->Activate();

	// clear the active rigid body BC's
	for (i=0; i<(int) m_fem.m_RB.size(); ++i)
	{
		FERigidBody& RB = m_fem.m_RB[i];
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_fem.GetMaterial(RB.m_mat));
		for (j=0; j<6; ++j)
		{
			if (RB.m_pDC[j])
			{
				RB.m_pDC[j] = 0;
				pm->m_bc[j] = 0;
			}
		}
	}

	// set the active rigid bodies BC's
	for (i=0; i<(int) m_fem.m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *(m_fem.m_RDC[i]);
		FERigidBody& RB = m_fem.m_RB[DC.id];
//		assert(RB.m_pDC[DC.bc] == 0);
		if (RB.m_bActive && DC.IsActive())
		{
			RB.m_pDC[DC.bc] = &DC;
			FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_fem.GetMaterial(RB.m_mat));
			pm->m_bc[DC.bc] = 1;
		}
	}

	// reset nodal ID's
	FEMesh& mesh = m_fem.m_mesh;
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (j=0; j<MAX_NDOFS; ++j)	node.m_ID[j] = node.m_BC[j];
	}

	// set the rigid nodes
	// Note that also the rotational degrees of freedom are fixed
	// for rigid nodes that do not belong to a non-rigid shell element.
	for (i=0; i<(int) m_fem.m_RN.size(); ++i)
	{
		FERigidNode& rn = *m_fem.m_RN[i];
		if (rn.IsActive())
		{
			FENode& node = m_fem.m_mesh.Node(rn.nid);
			node.m_rid = rn.rid;

			// fix degrees of freedom
			node.m_ID[0] = -1;
			node.m_ID[1] = -1;
			node.m_ID[2] = -1;
			if (node.m_bshell == false)
			{
				node.m_ID[3] = -1;
				node.m_ID[4] = -1;
				node.m_ID[5] = -1;
			}
		}
		else 
		{
			FENode& node = m_fem.m_mesh.Node(rn.nid);
			node.m_rid = -1;
		}
	}

	// override prescribed displacements for rigid nodes
	bool bdisp = false;
	for (i=0; i<(int) m_fem.m_DC.size(); ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			// if the node is not free we don't use this prescribed displacement
			// note that we don't do this for prescribed pressures and concentrations
			if ((dc.bc != 6) && (dc.bc != 11))
			{
				FENode& node = m_fem.m_mesh.Node(dc.node);
				if (node.m_rid >= 0) 
				{
					dc.Deactivate();
					bdisp = true;
				}
			}
		}
	}

	if (bdisp) clog.printbox("WARNING", "Rigid degrees of freedom cannot be prescribed.");

	// Sometimes an (ignorant) user might have added a rigid body
	// that is not being used. Since this can cause problems we need
	// to find these rigid bodies.
	vector<int> mec; mec.assign(m_fem.m_nrb, 0);
	for (i=0; i<m_fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = m_fem.m_mesh.Node(i);
		n = node.m_rid;
		if (n >= 0) mec[n]++;
	}

	for (i=0; i<m_fem.m_nrb; ++i)
		if (mec[i] == 0)
		{
			clog.printbox("WARNING", "Rigid body %d is not being used.", m_fem.m_RB[i].m_mat+1);
			m_fem.m_RB[i].m_bActive = false;
		}

	// initialize equations
	// ----->
	// TODO: Should I let the solver take care of this?
	if (m_fem.InitEquations() == false) return false;

	// initialize linear constraints
	// Must be done after equations are initialized
	if (m_fem.InitConstraints() == false) return false;
	// ----->

	// Now we adjust the equation numbers of prescribed dofs according to the above rule
	// Make sure that a prescribed dof has not been fixed
	int ndis = m_fem.m_DC.size();
	for (i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.m_DC[i];
		int nid = DC.node;
		int bc  = DC.bc;

		FENode& node = m_fem.m_mesh.Node(nid); 

		if (DC.IsActive())
		{
			switch (bc)
			{
			case 0: // x-displacement
			case 1: // y-displacement
			case 2: // z-displacement
			case 3: // x-rotation
			case 4: // y-rotation
			case 5: // z-rotation
			case 6: // prescribed pressure
			case 10: // precribed temperature
			case 11: // precribed concentration
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				break;
			case 20: // y-z radial displacement
				n = node.m_ID[1];
				node.m_ID[1] = (n<0?n:-n-2);
				n = node.m_ID[2];
				node.m_ID[2] = (n<0?n:-n-2);
				break;
			}
		}
	}

	// modify the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il)
		{
			list<FELinearConstraint::SlaveDOF>::iterator is = il->slave.begin();
			for (int i=0; i<(int) il->slave.size(); ++i, ++is)
			{
				is->neq = m_fem.m_mesh.Node(is->node).m_ID[is->bc];
			}
		}
	}

	// modify the (aug lag) linear constraints
	if (m_fem.m_LCSet.size())
	{
		int M = m_fem.m_LCSet.size();
		list<FELinearConstraintSet*>::iterator im = m_fem.m_LCSet.begin();
		for (int m=0; m<M; ++m, ++im) (*im)->Init();
	}

	// see if we need to do contact augmentations
	m_baugment = (m_fem.m_nrj > 0 ? true : false);
	for (i=0; i<(int) m_fem.m_CI.size(); ++i)
	{
		FEContactInterface& ci = *m_fem.m_CI[i];
		if (ci.m_blaugon) m_baugment = true;
	}

	// see if we to do incompressible augmentations
	for (i=0; i<(int) m_fem.m_MAT.size(); ++i)
	{
		FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_fem.GetMaterial(i));
		if (pmi && pmi->m_blaugon) m_baugment = true;
	}

	// see if we have to do linear constraint augmentations
	if (m_fem.m_LCSet.size()) m_baugment = true;

	return true;
}

//-----------------------------------------------------------------------------
bool FEAnalysis::Solve()
{
	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		clog.printbox("FATAL ERROR","Initialization has failed.\nAborting run.\n");
		return false;
	}

	// obtain a pointer to the console object. We'll use this to
	// set the title of the console window.
	Console* pShell = Console::GetHandle();

	// convergence flag
	// we initialize it to true so that when a restart is performed after 
	// the last time step we terminate normally.
	bool bconv = true;

	// calculate end time value
	double starttime = m_fem.m_ftime0;
	double endtime = m_fem.m_ftime0 + m_ntime*m_dt0;
	const double eps = endtime*1e-7;

	int nsteps = m_fem.m_Step.size();

	bool bdebug = m_fem.GetDebugFlag();

	if (nsteps > 1)
		pShell->SetTitle("(step %d/%d: %.f%%) %s - %s", m_fem.m_nStep+1, nsteps, (100.f*(m_fem.m_ftime - starttime) / (endtime - starttime)), m_fem.m_szfile_title, (bdebug?"FEBio (debug mode)": "FEBio"));
	else
		pShell->SetTitle("(%.f%%) %s - %s", (100.f*m_fem.m_ftime/endtime), m_fem.m_szfile_title, (bdebug?"FEBio (debug mode)": "FEBio"));

	// keep a stack for push/pop'ing
	stack<FEM> state;

	// print initial progress bar
	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		printf("\nProgress:\n");
		for (int i=0; i<50; ++i) printf("\xB0"); printf("\r");
		clog.SetMode(Logfile::FILE_ONLY);
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

	// repeat for all timesteps
	m_nretries = 0;
	while (endtime - m_fem.m_ftime > eps)
	{
		// keep a copy of the current state, in case
		// we need to retry this time step
		if (m_bautostep) 
		{
			while (!state.empty()) state.pop();
			state.push(m_fem);
		}

		// update time
		m_fem.m_ftime += m_dt;

		int i, j;

		// evaluate load curve values at current time
		for (i=0; i<m_fem.LoadCurves(); ++i) m_fem.GetLoadCurve(i)->Evaluate(m_fem.m_ftime);

		// evaluate parameter lists
		for (i=0; i<(int) m_fem.m_MPL.size(); ++i)
		{
			FEParameterList* pl = m_fem.m_MPL[i];
			list<FEParam>::iterator pi = pl->first();
			for (j=0; j<pl->Parameters(); ++j, ++pi)
			{
				if (pi->m_nlc >= 0)
				{
					double v = m_fem.GetLoadCurve(pi->m_nlc)->Value();
					switch (pi->m_itype)
					{
					case FE_PARAM_INT   : pi->value<int>() = (int) v; break;
					case FE_PARAM_DOUBLE: pi->value<double>() = v; break;
					case FE_PARAM_BOOL  : pi->value<bool>() = (v > 0? true : false); break;
					default: 
						assert(false);
					}
				}
			}
		}

		// solve this timestep,
		try
		{
			int oldmode = 0;
			bconv = m_psolver->SolveStep(m_fem.m_ftime);
		}
		catch (ExitRequest)
		{
			bconv = false;
			clog.printbox("WARNING", "Early termination on user's request");
			break;
		}
		catch (ZeroDiagonal e)
		{
			bconv = false;
			clog.printbox("FATAL ERROR", "%s", e.m_szerr);
			break;
		}
		catch (NANDetected)
		{
			bconv = false;
			clog.printbox("FATAL ERROR", "NAN Detected. Run aborted.");
			break;
		}
		catch (MemException e)
		{
			bconv = false;
			if (e.m_falloc < 1024*1024)
				clog.printbox("FATAL ERROR", "Failed allocating %lg bytes", e.m_falloc);
			else
			{
				double falloc = e.m_falloc / (1024.0*1024.0);
				clog.printbox("FATAL ERROR", "Failed allocating %lg MB", falloc);
			}
			break;
		}
		catch (std::bad_alloc e)
		{
			bconv = false;
			clog.printbox("FATAL ERROR", "A memory allocation failure has occured.\nThe program will now be terminated.");
			break;
		}
		catch (...)
		{
			bconv = false;
			clog.printbox("FATAL ERROR", "An unknown exception has occured.\nThe program will now be terminated.");
			break;
		}

		// update counters
		m_ntotref  += m_psolver->m_ntotref;
		m_ntotiter += m_psolver->m_niter;
		m_ntotrhs  += m_psolver->m_nrhs;

		// see if we have converged
		if (bconv)
		{
			// Yes! We have converged!
			clog.printf("\n\n------- converged at time : %lg\n\n", m_fem.m_ftime);

			// update nr of completed timesteps
			m_ntimesteps++;

			// output results to plot database
			if (m_nplot != FE_PLOT_NEVER)
			{
				if ((m_nplot == FE_PLOT_MUST_POINTS) && (m_nmplc >= 0))
				{
					FELoadCurve& lc = *m_fem.m_LC[m_nmplc];
					if (lc.HasPoint(m_fem.m_ftime)) m_fem.m_plot->Write(m_fem);
				}
				else m_fem.m_plot->Write(m_fem);
			}

			// Dump converged state to the archive
			if (m_bDump)
			{
				DumpFile ar(&m_fem);
				if (ar.Create(m_fem.m_szdump) == false)
				{
					clog.printf("WARNING: Failed creating restart point.\n");
				}
				else 
				{
					m_fem.Serialize(ar);
					clog.printf("\nRestart point created. Archive name is %s\n", m_fem.m_szdump);
				}
			}

			// store additional data to the logfile
			m_fem.m_Data.Write();

			// update time step
			if (m_bautostep && (m_fem.m_ftime + eps < endtime)) AutoTimeStep(m_psolver->m_niter);

			// reset retry counter
			m_nretries = 0;

			// call callback function
			m_fem.DoCallback();
		}
		else 
		{
			// Report the sad news to the user.
			clog.printf("\n\n------- failed to converge at time : %lg\n\n", m_fem.m_ftime);

			// If we have auto time stepping, decrease time step and let's retry
			if (m_bautostep && (m_nretries < m_maxretries))
			{
				// restore the previous state
				m_fem = state.top(); state.pop();
				
				// let's try again
				Retry();
			}
			else 
			{
				// can't retry, so abort
				if (m_nretries >= m_maxretries)	clog.printf("Max. nr of retries reached.\n\n");

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
		clog.flush();

		bool bdebug = m_fem.GetDebugFlag();
		if (nsteps>1)
			pShell->SetTitle("(step %d/%d: %.f%%) %s - %s", m_fem.m_nStep+1, nsteps, (100.f*(m_fem.m_ftime - starttime) / (endtime - starttime)), m_fem.m_szfile_title, (bdebug?"FEBio (debug mode)": "FEBio"));
		else
			pShell->SetTitle("(%.f%%) %s - %s", (100.f*m_fem.m_ftime/endtime), m_fem.m_szfile_title, (bdebug?"FEBio (debug mode)": "FEBio"));
	}

	m_fem.m_ftime0 = m_fem.m_ftime;

	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		clog.SetMode(Logfile::FILE_AND_SCREEN);
	}

	// output report
	clog.printf("\n\nN O N L I N E A R   I T E R A T I O N   I N F O R M A T I O N\n\n");
	clog.printf("\tNumber of time steps completed .................... : %d\n\n", m_ntimesteps);
	clog.printf("\tTotal number of equilibrium iterations ............ : %d\n\n", m_ntotiter);
	clog.printf("\tAverage number of equilibrium iterations .......... : %lg\n\n", (double) m_ntotiter / (double) m_ntimesteps);
	clog.printf("\tTotal number of right hand evaluations ............ : %d\n\n", m_ntotrhs);
	clog.printf("\tTotal number of stiffness reformations ............ : %d\n\n", m_ntotref);

	// get and print elapsed time
	char sztime[64];

	m_psolver->m_SolverTime.time_str(sztime);
	clog.printf("\tTime in solver: %s\n\n", sztime);

	return bconv;
}

//-----------------------------------------------------------------------------
void FEAnalysis::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		// --- analysis data ---
		ar << m_nModule;
		ar << m_nanalysis;
		ar << m_istiffpr;
		ar << m_baugment;
		ar << m_hg;

		// --- Time Step Data ---
		ar << m_ntime;
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
	}
	else
	{
		// --- analysis data ---
		ar >> m_nModule;
		ar >> m_nanalysis;
		ar >> m_istiffpr;
		ar >> m_baugment;
		ar >> m_hg;

		// --- Time Step Data ---
		ar >> m_ntime;
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

		// create a solver
		assert(m_psolver == 0);
		switch (m_nModule)
		{
		case FE_SOLID: 
			m_psolver = new FESolidSolver(m_fem); 
			break;
		case FE_POROELASTIC:
			m_psolver = new FEPoroSolidSolver(m_fem);
			break;
		case FE_POROSOLUTE:
			m_psolver = new FEPoroSoluteSolver(m_fem);
			break;
		case FE_HEAT:
			m_psolver = new FEHeatSolver(m_fem);
			break;
		default:
			throw "Unknown module type in FEAnalysis::Serialize";
		}
	}

	// Seriaize solver data
	m_psolver->Serialize(ar);
}

//-----------------------------------------------------------------------------
//! Restores data for a running restart

void FEAnalysis::Retry()
{
	clog.printf("Retrying time step. Retry attempt %d of max %d\n\n", m_nretries+1, m_maxretries);

	// adjust time step
	static double ddt = 0;
	double dtn;

	if (m_nretries == 0) ddt = (m_dt) / (m_maxretries+1);

	if (m_naggr == 0) dtn = m_dt - ddt;
	else dtn = m_dt*0.5;

	clog.printf("\nAUTO STEPPER: retry step, dt = %lg\n\n", dtn);

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

	// get the must point load curve
	FELoadCurve& lc = *m_fem.m_LC[ m_nmplc ];

	double dtmax = lc.Value(told);

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
			clog.printf("\nAUTO STEPPER: increasing time step, dt = %lg\n\n", dtn);
		else if (dtn < m_dt)
			clog.printf("\nAUTO STEPPER: decreasing time step, dt = %lg\n\n", dtn);
	}

	// check for mustpoints
	double tnew = told + dtn;

	const double eps = m_tend*1e-07;

	int n = lc.FindPoint(told+eps);
	if (n >= 0)
	{
		// TODO: what happens when dtn < dtmin and the next time step fails??
		double tmust = lc.LoadPoint(n).time;
		if (tnew > tmust)
		{
			dtn = tmust - told;
			clog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
		}
		else if (tnew > m_tend)
		{
			dtn = m_tend - told;
			clog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
		}
	}

	m_dt = dtn;
}
