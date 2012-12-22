#include "stdafx.h"
#include "FEAnalysisStep.h"
#include <FECore/FERigid.h>
#include "FEUncoupledMaterial.h"
#include "FESolidSolver.h"
#include "FEHeatSolver.h"
#include "FEBiphasicSolver.h"
#include "FEBiphasicSoluteSolver.h"
#include "FELinearSolidSolver.h"
#include "FECoupledHeatSolidSolver.h"
#include "FEExplicitSolidSolver.h"
#include <FECore/FERigidBody.h>
#include "log.h"

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

//-----------------------------------------------------------------------------
//! constructor
FEAnalysisStep::FEAnalysisStep(FEModel& fem, int ntype) : FEAnalysis(fem, ntype)
{
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
FEAnalysisStep::~FEAnalysisStep(void)
{
}

//-----------------------------------------------------------------------------

void FEAnalysisStep::Finish()
{
	// deactivate the boundary conditions
	for (size_t i=0; i<m_BC.size(); ++i) m_BC[i]->Deactivate();

	// deactivate contact interfaces
	for (size_t i=0; i<m_CI.size(); ++i) m_CI[i]->Deactivate();

	// clean up solver data (i.e. destroy linear solver)
	m_psolver->Clean();
}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FEAnalysisStep::Init()
{
	int i, j, n;

	// set first time step
	// We can't do this since it will mess up the value from a restart
//	m_dt = m_dt0;

	double Dt;
	if (m_ntime == -1) Dt = m_final_time; else Dt = m_dt0*m_ntime;

	m_tend = m_fem.m_ftime0 + Dt;

	// For now, add all domains to the analysis step
	FEMesh& mesh = m_fem.GetMesh();
	int ndom = mesh.Domains();
	ClearDomains();
	for (i=0; i<ndom; ++i) AddDomain(i);

	// activate the boundary conditions
	for (i=0; i<(int) m_BC.size(); ++i) m_BC[i]->Activate();

	// activate contact interface
	for (i=0; i<(int) m_CI.size(); ++i) m_CI[i]->Activate();

	// clear the active rigid body BC's
	int NRB = m_fem.Objects();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.Object(i));
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
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.Object(DC.id));
//		assert(RB.m_pDC[DC.bc] == 0);
		if (RB.IsActive() && DC.IsActive())
		{
			RB.m_pDC[DC.bc] = &DC;
			FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_fem.GetMaterial(RB.m_mat));
			pm->m_bc[DC.bc] = 1;
		}
	}

	// reset nodal ID's
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (j=0; j<MAX_NDOFS; ++j)	node.m_ID[j] = node.m_BC[j];
	}

	// set the rigid nodes
	// Note that also the rotational degrees of freedom are fixed
	// for rigid nodes that do not belong to a non-rigid shell element.
	int nrn = m_fem.RigidNodes();
	for (i=0; i<nrn; ++i)
	{
		FERigidNode& rn = *m_fem.RigidNode(i);
		if (rn.IsActive())
		{
			FENode& node = m_fem.GetMesh().Node(rn.nid);
			node.m_rid = rn.rid;

			// fix degrees of freedom
			node.m_ID[DOF_X] = -1;
			node.m_ID[DOF_Y] = -1;
			node.m_ID[DOF_Z] = -1;
			if (node.m_bshell == false)
			{
				node.m_ID[DOF_U] = -1;
				node.m_ID[DOF_V] = -1;
				node.m_ID[DOF_W] = -1;
			}
		}
		else 
		{
			FENode& node = m_fem.GetMesh().Node(rn.nid);
			node.m_rid = -1;
		}
	}

	// override prescribed displacements for rigid nodes
	bool bdisp = false;
	int nbc = m_fem.PrescribedBCs();
	for (i=0; i<nbc; ++i)
	{
		FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
		if (dc.IsActive())
		{
			// if the node is not free we don't use this prescribed displacement
			// note that we don't do this for prescribed pressures and concentrations
			if ((dc.bc != DOF_P) && (dc.bc < DOF_C))
			{
				FENode& node = m_fem.GetMesh().Node(dc.node);
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
	int nrb = m_fem.Objects();
	vector<int> mec; mec.assign(nrb, 0);
	for (i=0; i<m_fem.GetMesh().Nodes(); ++i)
	{
		FENode& node = m_fem.GetMesh().Node(i);
		n = node.m_rid;
		if (n >= 0) mec[n]++;
	}

	for (i=0; i<nrb; ++i)
		if (mec[i] == 0)
		{
			FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_fem.Object(i));
			clog.printbox("WARNING", "Rigid body %d is not being used.", RB.m_mat+1);
			RB.Activate(false);
		}

	// initialize equations
	// ----->
	if (m_psolver->InitEquations() == false) return false;

	// initialize linear constraints
	// Must be done after equations are initialized
	if (InitConstraints() == false) return false;
	// ----->

	// Now we adjust the equation numbers of prescribed dofs according to the above rule
	// Make sure that a prescribed dof has not been fixed
	// TODO: maybe this can be moved to the FESolver::InitEquations function
	int ndis = m_fem.PrescribedBCs();
	for (i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		int nid = DC.node;
		int bc  = DC.bc;
		bool br = DC.br;

		FENode& node = m_fem.GetMesh().Node(nid); 

		if (DC.IsActive())
		{
			switch (bc)
			{
			case DOF_X: // x-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.x - node.m_r0.x : 0;	// GAA
				break;
			case DOF_Y: // y-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.y - node.m_r0.y : 0;
				break;
			case DOF_Z: // z-displacement
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_rt.z - node.m_r0.z : 0;
				break;
			case DOF_U: // x-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.x - node.m_D0.x : 0;
				break;
			case DOF_V: // y-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.y - node.m_D0.y : 0;
				break;
			case DOF_W: // z-rotation
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_Dt.z - node.m_D0.z : 0;
				break;
			case DOF_P: // prescribed pressure
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_pt - node.m_p0 : 0;
				break;
			case DOF_T: // precribed temperature
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = 0;
				break;
/*			case DOF_C: // precribed concentration
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_ct[0] - node.m_c0[0] : 0;
				break;
			case DOF_C+1: // precribed concentration
				n = node.m_ID[bc];
				node.m_ID[bc] = (n<0?n:-n-2);
				DC.r = br ? node.m_ct[1] - node.m_c0[1] : 0;
				break;*/
				//--> TODO: change bc=20 to something else
			case 20: // y-z radial displacement
				n = node.m_ID[DOF_Y];
				node.m_ID[DOF_Y] = (n<0?n:-n-2);
				n = node.m_ID[DOF_Z];
				node.m_ID[DOF_Z] = (n<0?n:-n-2);
				DC.r = 0;
				break;
			default:	// all prescribed concentrations
				if ((bc >= DOF_C) && (bc < MAX_NDOFS)) {
					n = node.m_ID[bc];
					node.m_ID[bc] = (n<0?n:-n-2);
					int sid = bc - DOF_C;
					DC.r = br ? node.m_ct[sid] - node.m_c0[sid] : 0;
				}
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
				is->neq = m_fem.GetMesh().Node(is->node).m_ID[is->bc];
			}
		}
	}

	// modify the (aug lag) nonlinear constraints
	// TODO: I think this is already done in FEM::Init. Why do I need to do this again?
	int M = m_fem.NonlinearConstraints();
	for (int m=0; m<M; ++m) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(m);
		plc->Init();
	}

	// see if we need to do contact augmentations
	m_baugment = false;
	for (i=0; i<m_fem.ContactInterfaces(); ++i)
	{
		FEContactInterface& ci = *m_fem.ContactInterface(i);
		if (ci.IsActive() && ci.m_blaugon) m_baugment = true;
	}

	// see if we need to do incompressible augmentations
	int nmat = m_fem.Materials();
	for (i=0; i<nmat; ++i)
	{
		FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_fem.GetMaterial(i));
		if (pmi && pmi->m_blaugon) m_baugment = true;
	}

	// see if we have to do nonlinear constraint augmentations
	if (m_fem.NonlinearConstraints() != 0) m_baugment = true;

	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		clog.printbox("FATAL ERROR","Failed to initialize solver.\nAborting run.\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function initializes the linear constraint table (LCT). This table
//! contains for each dof the linear constraint it belongs to. (or -1 if it is
//! not constraint)

bool FEAnalysisStep::InitConstraints()
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
bool FEAnalysisStep::Solve(Progress& prg)
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

	// set progress
	prg.SetProgress(100.f*(m_fem.m_ftime - starttime) / (endtime - starttime));

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
		if (m_bautostep) m_fem.PushState();

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
		for (i=0; i<m_fem.BodyForces(); ++i)
		{
			FEParameterList& pl = m_fem.GetBodyForce(i)->GetParameterList();
			m_fem.EvaluateParameterList(pl);
		}

		// evaluate contact interface parameter lists
		for (i=0; i<m_fem.ContactInterfaces(); ++i)
		{
			FEParameterList& pl = m_fem.ContactInterface(i)->GetParameterList();
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
			// TODO: Fix this feature
			clog.printbox("FATAL ERROR", "Zero diagonal detected. Aborting run.");
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
		catch (FEMultiScaleException)
		{
			bconv = false;
			clog.printbox("FATAL ERROR", "The RVE problem has failed. Aborting macro run.");
			break;
		}
		catch (std::bad_alloc e)
		{
			bconv = false;
			clog.printbox("FATAL ERROR", "A memory allocation failure has occured.\nThe program will now be terminated.");
			break;
		}

// We only catch all exceptions for release versions
#ifndef _DEBUG
		catch (...)
		{
			bconv = false;
			clog.printbox("FATAL ERROR", "An unknown exception has occured.\nThe program will now be terminated.");
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
			clog.printf("\n\n------- converged at time : %lg\n\n", m_fem.m_ftime);

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
				m_fem.PopState();
				
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

		// set progress
		prg.SetProgress(100.f*(m_fem.m_ftime - starttime) / (endtime - starttime));
	}

	if ((m_nplot == FE_PLOT_FINAL) && bconv) m_fem.Write();

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

	FESolver* psolver = dynamic_cast<FESolver*>(m_psolver);
	psolver->m_SolverTime.time_str(sztime);
	clog.printf("\tTime in solver: %s\n\n", sztime);

	return bconv;
}

//-----------------------------------------------------------------------------
void FEAnalysisStep::Serialize(DumpFile& ar)
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

		// create a solver
		assert(m_psolver == 0);
		switch (m_ntype)
		{
		case FE_SOLID         : m_psolver = new FESolidSolver           (m_fem); break;
		case FE_EXPLICIT_SOLID: m_psolver = new FEExplicitSolidSolver   (m_fem); break;
		case FE_BIPHASIC	  : m_psolver = new FEBiphasicSolver        (m_fem); break;
		case FE_POROSOLUTE    : m_psolver = new FEBiphasicSoluteSolver  (m_fem); break;
		case FE_HEAT          : m_psolver = new FEHeatSolver            (m_fem); break;
		case FE_LINEAR_SOLID  : m_psolver = new FELinearSolidSolver     (m_fem); break;
		case FE_HEAT_SOLID    : m_psolver = new FECoupledHeatSolidSolver(m_fem); break;
		default:
			throw "Unknown module type in FEAnalysisStep::Serialize";
		}
	}

	// Seriaize solver data
	m_psolver->Serialize(ar);
}

//-----------------------------------------------------------------------------
//! Restores data for a running restart

void FEAnalysisStep::Retry()
{
	clog.printf("Retrying time step. Retry attempt %d of max %d\n\n", m_nretries+1, m_maxretries);

	// adjust time step
	double dtn;

	if (m_nretries == 0) m_ddt = (m_dt) / (m_maxretries+1);

	if (m_naggr == 0) dtn = m_dt - m_ddt;
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

void FEAnalysisStep::AutoTimeStep(int niter)
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
			clog.printf("\nAUTO STEPPER: increasing time step, dt = %lg\n\n", dtn);
		else if (dtn < m_dt)
			clog.printf("\nAUTO STEPPER: decreasing time step, dt = %lg\n\n", dtn);
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
				clog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
			}
			else if (tnew > m_tend)
			{
				dtn = m_tend - told;
				clog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
			}
		}
	}

	// make sure we are not exceeding the final time
	if (told + dtn > m_tend)
	{
		dtn = m_tend - told;
		clog.printf("MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
	}

	// store time step size
	m_dt = dtn;
}
