#include "StdAfx.h"
#include "FEAnalysis.h"
#include "fem.h"
#include "console.h"

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

//-----------------------------------------------------------------------------
//! constructor
FEAnalysis::FEAnalysis(FEM& fem) : m_fem(fem)
{
	// --- Analysis data ---
	m_itype = FE_STATIC;	// do quasi-static analysis
	m_istiffpr = 1;			// use pressure stiffness
	m_baugment = false;		// no augmentations

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
	m_psolver = new FESolver(m_fem);
	m_nretries = 0;
	m_maxretries = 5;

	// --- I/O Data ---
	m_bDump = false;
	m_nplot = FE_PLOT_MAJOR_ITRS;
	m_nprint = FE_PRINT_MINOR_ITRS;

	m_BC.setgrowsize(1024);
}

//-----------------------------------------------------------------------------
FEAnalysis::~FEAnalysis(void)
{
}

//-----------------------------------------------------------------------------

void FEAnalysis::Finish()
{
	// deactivate the boundary conditions
	for (int i=0; i<m_BC.size(); ++i) m_BC[i]->Deactivate();
}

//-----------------------------------------------------------------------------
bool FEAnalysis::Init()
{
	int i, j, n;

	// initialize counters
	m_ntotref    = 0;		// total nr of stiffness reformations
	m_ntotiter   = 0;		// total nr of non-linear iterations
	m_ntimesteps = 0;		// time steps completed
	m_ntotrhs    = 0;		// total nr of right hand side evaluations

	// set first time step
	m_dt = m_dt0;

	m_tend = m_fem.m_ftime + m_dt0*m_ntime;

	// init must point curve
	if (m_nmplc < 0)
	{
		FELoadCurve* plc = new FELoadCurve;
		plc->Create(2);
		plc->LoadPoint(0).time  = m_fem.m_ftime;
		plc->LoadPoint(0).value = 0;
		plc->LoadPoint(1).time  = m_fem.m_ftime + m_dt*m_ntime;
		plc->LoadPoint(1).value = m_dtmax;
		m_fem.m_LC.add(plc);
		m_nmplc = m_fem.m_LC.size()-1;
	}

	// the must point load curve must be evaluated
	// using a step interpolation
	m_fem.m_LC[m_nmplc].SetInterpolation(FELoadCurve::STEP);

	// activate the boundary conditions
	for (i=0; i<m_BC.size(); ++i) m_BC[i]->Activate();

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
	for (i=0; i<m_fem.m_RN.size(); ++i)
	{
		FERigidNode& rn = m_fem.m_RN[i];
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
	for (i=0; i<m_fem.m_DC.size(); ++i)
	{
		FENodalDisplacement& dc = m_fem.m_DC[i];

		if (dc.IsActive())
		{
			// if the node is not free we don't use this prescribed displacement
			// note that we don't do this for prescribed pressures
			if (dc.bc != 6)
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
	if (bdisp) m_fem.m_log.printbox("WARNING", "Rigid degrees of freedom cannot be prescribed.");

	// Sometimes an (ignorant) user might have added a rigid body
	// that is not being used. Since this can cause problems we need
	// to find these rigid bodies.
	vector<int> mec(m_fem.m_nrb); mec.zero();
	for (i=0; i<m_fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = m_fem.m_mesh.Node(i);
		n = node.m_rid;
		if (n >= 0) mec[n]++;
	}

	for (i=0; i<m_fem.m_nrb; ++i)
		if (mec[i] == 0)
		{
			m_fem.m_log.printbox("WARNING", "Rigid body %d is not being used.", m_fem.m_RB[i].m_mat+1);
			m_fem.m_RB[i].m_bc[0] = -1;
			m_fem.m_RB[i].m_bc[1] = -1;
			m_fem.m_RB[i].m_bc[2] = -1;
			m_fem.m_RB[i].m_bc[3] = -1;
			m_fem.m_RB[i].m_bc[4] = -1;
			m_fem.m_RB[i].m_bc[5] = -1;
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
		FENodalDisplacement& DC = m_fem.m_DC[i];
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
		for (int l=0; l<m_fem.m_LinC.size(); ++l, ++il)
		{
			list<FELinearConstraint::SlaveDOF>::iterator is = il->slave.begin();
			for (int i=0; i<il->slave.size(); ++i, ++is)
			{
				is->neq = m_fem.m_mesh.Node(is->node).m_ID[is->bc];
			}
		}
	}

	// modify the (aug lag) linear constraints
	if (m_fem.m_LCAL.size())
	{
		list<FEAugLagLinearConstraint>::iterator il = m_fem.m_LCAL.begin();
		FEAugLagLinearConstraint::DOF& ms = il->m_master;
		ms.neq = m_fem.m_mesh.Node(ms.node).m_ID[ms.bc];
		for (int l=0; l<m_fem.m_LCAL.size(); ++l, ++il)
		{
			list<FEAugLagLinearConstraint::SlaveDOF>::iterator is = il->m_slave.begin();
			for (int i=0; i<il->m_slave.size(); ++i, ++is)
			{
				is->neq = m_fem.m_mesh.Node(is->node).m_ID[is->bc];
			}
		}
	}


	// see if we need to do contact augmentations
	m_baugment = (m_fem.m_nrj > 0 ? true : false);
	for (i=0; i<m_fem.m_CI.size(); ++i)
	{
		FESlidingInterface* ps = dynamic_cast<FESlidingInterface*>(&m_fem.m_CI[i]);
		if (ps)
		{
			if (ps->m_blaugon) m_baugment = true;
		}
		else m_baugment = true;
	}

	// see if we to do incompressible augmentations
	for (i=0; i<m_fem.m_MAT.size(); ++i)
	{
		FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(m_fem.GetMaterial(i));
		if (pmi && pmi->m_blaugon) m_baugment = true;
	}


	return true;
}

//-----------------------------------------------------------------------------
bool FEAnalysis::Solve()
{
	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		m_fem.m_log.printbox("FATAL ERROR","Initialization has failed.\nAborting run.\n");
		return false;
	}

	// convergence flag
	// we initialize it to true so that when a restart is performed after 
	// the last time step we terminate normally.
	bool bconv = true;

	// calculate end time value
	double endtime = m_fem.m_ftime + m_ntime*m_dt0;
	const double eps = endtime*1e-7;

	Console::SetTitle("(%.f%%) %s - FEBio", (100.f*m_fem.m_ftime / endtime), m_fem.m_szfile_title);

	// make sure that the timestep is at least the min time step size
	if (m_bautostep) AutoTimeStep(0);

	// keep a stack for push/pop'ing
	stack<FEM> state(1);

	// print initial progress bar
	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		printf("\nProgress:\n");
		for (int i=0; i<50; ++i) printf("\xB0"); printf("\r");
		m_fem.m_log.SetMode(Logfile::FILE_ONLY);
	}

	// repeat for all timesteps
	m_nretries = 0;
	while (endtime - m_fem.m_ftime > eps)
	{
		// keep a copy of the current state, in case
		// we need to retry this time step
		if (m_bautostep) 
		{
			state.Clear();
			state.Push(m_fem);
		}

		// update time
		m_fem.m_ftime += m_dt;

		m_fem.m_log.printf("\n===== beginning time step %d : %lg =====\n", m_ntimesteps+1, m_fem.m_ftime);

		// solve this timestep,
		// that is, use the Newton Raphson method to solve the timestep
		try
		{
			int oldmode = 0;
			bconv = m_psolver->SolveStep(m_fem.m_ftime);
		}
		catch (ExitRequest)
		{
			bconv = false;
			m_fem.m_log.printbox("WARNING", "Early termination on user's request");
			break;
		}
		catch (ZeroDiagonal e)
		{
			bconv = false;
			m_fem.m_log.printbox("FATAL ERROR", "%s", e.m_szerr);
			break;
		}
		catch (NANDetected)
		{
			bconv = false;
			m_fem.m_log.printbox("FATAL ERROR", "NAN Detected. Run aborted.");
			break;
		}
		catch (MemException e)
		{
			bconv = false;
			if (e.m_falloc < 1024*1024)
				m_fem.m_log.printbox("FATAL ERROR", "Failed allocating %lg bytes", e.m_falloc);
			else
			{
				double falloc = e.m_falloc / (1024.0*1024.0);
				m_fem.m_log.printbox("FATAL ERROR", "Failed allocating %lg MB", falloc);
			}
			break;
		}
		catch (...)
		{
			bconv = false;
			m_fem.m_log.printbox("FATAL ERROR", "An unknown exception has occured.\nThe program will now be terminated.");
			break;
		}

		// update counters
		m_ntotref  += m_psolver->m_nref;
		m_ntotiter += m_psolver->m_niter;
		m_ntotrhs  += m_psolver->m_nrhs;

		// see if we have converged
		if (bconv)
		{
			// Yes! We have converged!
			m_fem.m_log.printf("\n\n------- converged at time : %lg\n\n", m_fem.m_ftime);

			// update nr of completed timesteps
			m_ntimesteps++;

			// update time step
			if (m_bautostep && (m_fem.m_ftime + eps < endtime)) AutoTimeStep(m_psolver->m_niter);

			// reset retry counter
			m_nretries = 0;

			// output results to plot database
			if (m_nplot != FE_PLOT_NEVER)
			{
				if ((m_nplot == FE_PLOT_MUST_POINTS) && (m_nmplc >= 0))
				{
					FELoadCurve& lc = m_fem.m_LC[m_nmplc];
					if (lc.HasPoint(m_fem.m_ftime)) m_fem.m_plot.Write(m_fem);
				}
				else m_fem.m_plot.Write(m_fem);
			}

			// Dump converged state to the archive
			if (m_bDump)
			{
				Archive ar;
				if (ar.Create(m_fem.m_szdump) == false)
				{
					m_fem.m_log.printf("WARNING: Failed creating restart point.\n");
				}
				else 
				{
					m_fem.Serialize(ar);
					m_fem.m_log.printf("\nRestart point created. Archive name is %s\n", m_fem.m_szdump);
				}
			}

			// store additional data to the logfile
			m_fem.m_Data.Write();

			// call callback function
			m_fem.DoCallback();
		}
		else 
		{
			// Report the sad news to the user.
			m_fem.m_log.printf("\n\n------- failed to converge at time : %lg\n\n", m_fem.m_ftime);

			// If we have auto time stepping, decrease time step and let's retry
			if (m_bautostep && (m_nretries < m_maxretries))
			{
				// restore the previous state
				state.Pop(m_fem);
				
				// let's try again
				Retry();
			}
			else 
			{
				// can't retry, so abort
				if (m_nretries >= m_maxretries)	m_fem.m_log.printf("Max. nr of retries reached.\n\n");

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
		m_fem.m_log.flush();

		Console::SetTitle("(%.f%%) %s - FEBio", (100.f*m_fem.m_ftime / endtime), m_fem.m_szfile_title);
	}

	if (GetPrintLevel() == FE_PRINT_PROGRESS)
	{
		m_fem.m_log.SetMode(Logfile::FILE_AND_SCREEN);
	}

	// output report
	m_fem.m_log.printf("\n\nN O N L I N E A R   I T E R A T I O N   I N F O R M A T I O N\n\n");
	m_fem.m_log.printf("\tNumber of time steps completed .................... : %d\n\n", m_ntimesteps);
	m_fem.m_log.printf("\tTotal number of equilibrium iterations ............ : %d\n\n", m_ntotiter);
	m_fem.m_log.printf("\tAverage number of equilibrium iterations .......... : %lg\n\n", (double) m_ntotiter / (double) m_ntimesteps);
	m_fem.m_log.printf("\tTotal number of right hand evaluations ............ : %d\n\n", m_ntotrhs);
	m_fem.m_log.printf("\tTotal number of stiffness reformations ............ : %d\n\n", m_ntotref);

	// get and print elapsed time
	char sztime[64];

	m_psolver->m_SolverTime.time_str(sztime);
	m_fem.m_log.printf("\tTime in solver: %s\n\n", sztime);

	return bconv;
}

//-----------------------------------------------------------------------------
void FEAnalysis::Serialize(Archive& ar)
{
	// TODO:	serialize the boundary conditions
	//			not sure how to do this yet.
	if (ar.IsSaving())
	{
		// --- analysis data ---
		ar << m_itype;
		ar << m_istiffpr;
		ar << m_baugment;

		// --- Time Step Data ---
		ar << m_ntime;
		ar << m_dt;
		ar << m_dt0;
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
		ar << m_bDump;
		ar << m_nplot;
		ar << m_nprint;
	}
	else
	{
		// --- analysis data ---
		ar >> m_itype;
		ar >> m_istiffpr;
		ar >> m_baugment;

		// --- Time Step Data ---
		ar >> m_ntime;
		ar >> m_dt;
		ar >> m_dt0;
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
		ar >> m_bDump;
		ar >> m_nplot;
		ar >> m_nprint;
	}

	// serialize solver data
	m_psolver->Serialize(ar);
}

//-----------------------------------------------------------------------------
//! Restores data for a running restart

void FEAnalysis::Retry()
{
	m_fem.m_log.printf("Retrying time step. Retry attempt %d of max %d\n\n", m_nretries+1, m_maxretries);

	// adjust time step
	static double ddt = 0;
	double dtn;

	if (m_nretries == 0) ddt = (m_dt) / (m_maxretries+1);

	if (m_naggr == 0) dtn = m_dt - ddt;
	else dtn = m_dt*0.5;

	m_fem.m_log.printf("\nAUTO STEPPER: retry step, dt = %lg\n\n", dtn);

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
	FELoadCurve& lc = m_fem.m_LC[ m_nmplc ];

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
			m_fem.m_log.printf("\nAUTO STEPPER: increasing time step, dt = %lg\n\n", dtn);
		else if (dtn < m_dt)
			m_fem.m_log.printf("\nAUTO STEPPER: decreasing time step, dt = %lg\n\n", dtn);
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
			m_fem.m_log.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
		}
		else if (tnew > m_tend)
		{
			dtn = m_tend - told;
			m_fem.m_log.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
		}
	}

	m_dt = dtn;
}
