#include "stdafx.h"
#include "fem.h"
#include "FECore/FEException.h"
#include "FEBioLib/FESolidSolver.h"
#include "FEBioLib/FETransverselyIsotropic.h"
#include "FEBioLib/FEDiscreteMaterial.h"
#include "FEBioLib/FEElasticSolidDomain.h"
#include "FEBioLib/FEElasticShellDomain.h"
#include "FEBioLib/FEPointBodyForce.h"
#include "FEBioLib/FERigidJoint.h"
#include "FEBioLib/FERigidBody.h"
#include "FEBioLib/FERigid.h"
#include "FEBioLib/FERigidJoint.h"
#include "FEBioLib/FETriphasic.h"
#include "FEBioLib/log.h"
#include "FEBioPlot/LSDYNAPlotFile.h"
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

//-----------------------------------------------------------------------------
// Forward declarations
void Hello(FILE* fp);

//-----------------------------------------------------------------------------
//! This function performs one-time-initialization stuff. All the different 
//! modules are initialized here as well. This routine also performs some
//! data checks

bool FEM::Init()
{
	int i;

	// Open the logfile
	if (!clog.is_valid()) 
	{
		if (clog.open(m_szlog) == false)
		{
			clog.printbox("FATAL ERROR", "Failed creating log file");
			return false;
		}

		// if we don't want to output anything we only output to the logfile
		if (m_pStep->GetPrintLevel() == FE_PRINT_NEVER) clog.SetMode(Logfile::FILE_ONLY);

		// print welcome message to file
		Hello(clog);
	}

	// check step data
	// TODO: should I let the Steps take care of this instead?
	for (i=0; i<(int) m_Step.size(); ++i)
	{
		FEAnalysis& step = *m_Step[i];
		if ((step.m_ntime <= 0) && (step.m_final_time <= 0.0)) { clog.printf("Invalid number of time steps for analysis step %d", i+1); return false; }
		if ((step.m_ntime >  0) && (step.m_final_time >  0.0)) { clog.printf("You must either set the number of time steps or the final time but not both.\n"); return false; }
		if (step.m_dt0   <= 0) { clog.printf("Invalid time step size for analysis step %d", i+1); return false; }
		if (step.m_bautostep)
		{
//			if (m_pStep->m_dtmin <= 0) return err("Invalid minimum time step size");
//			if (m_pStep->m_dtmax <= 0) return err("Invalid maximum time step size");
		}
	}

	// evaluate all loadcurves at the initial time
	for (i=0; i<LoadCurves(); ++i) m_LC[i]->Evaluate(0);

	// if the analysis is run in plain-strain mode we fix all the z-dofs of all nodes
	if (m_nplane_strain >= 0)
	{
		int bc = m_nplane_strain;
		for (int i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[bc] = -1;
	}

	// find and remove isolated vertices
	int ni = m_mesh.RemoveIsolatedVertices();
	if (ni != 0) 
	{
		if (ni == 1)
			clog.printbox("WARNING", "%d isolated vertex removed.", ni);
		else
			clog.printbox("WARNING", "%d isolated vertices removed.", ni);
	}

	// create and initialize the rigid body data
	if (CreateRigidBodies() == false) return false;

	// initialize poroelastic/biphasic and solute data
	if (InitPoroSolute() == false) return false;

	// initialize random number generator
	srand((unsigned) time(NULL));

	// initialize mesh data
	// note that this must be done AFTER the elements have been assigned material point data !
	// this is because the mesh data is reset
	// TODO: perhaps I should not reset the mesh data during the initialization
	if (InitMesh() == false) return false;

	// initialize material data
	if (InitMaterials() == false) return false;

	// initialize contact data
	if (InitContact() == false) return false;

	// intitialize time
	m_ftime = 0;
	m_ftime0 = 0;

	// init some other stuff
	for (i=0; i<(int) m_BF.size(); ++i)
	{
		if (m_BF[i]->Init() == false) return false;
	}

	// initialize nonlinear constraints
	// TODO: This is also initialized in the analysis step. Do I need to do this here?
	for (i=0; i<(int) m_NLC.size(); ++i) m_NLC[i]->Init();

	// open plot database file
	if (m_pStep->m_nplot != FE_PLOT_NEVER)
	{
		if (m_plot == 0) m_plot = new LSDYNAPlotFile;

		if (m_plot->Open(*this, m_szplot) == false)
		{
			clog.printf("ERROR : Failed creating PLOT database\n");
			return false;
		}
	}

	// Since it is assumed that for the first timestep
	// there are no loads or initial displacements, the case n=0 is skipped.
	// Therefor we can output those results here.
	// Offcourse we should actually check if this is indeed
	// the case, otherwise we should also solve for t=0
	if (m_pStep->m_nplot != FE_PLOT_NEVER) m_plot->Write(*this);

	// do the callback
	DoCallback();

	// Alright, all initialization is done, so let's get busy !
	return true;
}

//-----------------------------------------------------------------------------
//! Update the mesh data. This function calculates the initial directors
//! for the shell elements.

// NOTE: This function needs to be called after the rigid bodies have been
// initialized

bool FEM::InitMesh()
{
	int i, j, n, m0, m1, m2, nd;
	int* en;
	vec3d a, b, c;

	int ninverted = 0;

	FEMesh& m = m_mesh;

	// zero initial directors for shell nodes
	for (i=0; i<m.Nodes(); ++i) m.Node(i).m_D0 = vec3d(0,0,0);

	for (nd = 0; nd < m.Domains(); ++nd)
	{
		// check all solid elements to see if they are not initially inverted
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);

				int nint = el.GaussPoints();
				for (int n=0; n<nint; ++n)
				{
					double J0 = pbd->detJ0(el, n);
					if (J0 <= 0)
					{
						clog.printf("**************************** E R R O R ****************************\n");
						clog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.m_nID);
						clog.printf("Jacobian = %lg\n", J0);
						clog.printf("Did you use the right node numbering?\n");
						clog.printf("Nodes:");
						for (int l=0; l<el.Nodes(); ++l)
						{
							clog.printf("%d", el.m_node[l]+1);
							if (l+1 != el.Nodes()) clog.printf(","); else clog.printf("\n");
						}
						clog.printf("*******************************************************************\n\n");
						++ninverted;
					}
				}
			}
		}

		// Calculate the shell directors as the local node normals
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			vec3d r0[4];
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);

				n = el.Nodes();
				en = &el.m_node[0];

				// get the nodes
				for (j=0; j<n; ++j) r0[j] = psd->GetMesh()->Node(en[j]).m_r0;

				for (j=0; j<n; ++j)
				{
					m0 = j;
					m1 = (j+1)%n;
					m2 = (j==0? n-1: j-1);

					a = r0[m0];
					b = r0[m1];
					c = r0[m2];

					m.Node(en[m0]).m_D0 += (b-a)^(c-a);
				}
			}
		}
	}

	// make sure we start with unit directors
	for (i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		node.m_D0.unit();
		node.m_Dt = node.m_D0;
	}

	// Check initially inverted shell elements
	for (nd = 0; nd < m.Domains(); ++nd)
	{
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			// check the connectivity of the shells
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				if (!el.IsRigid())
				{
					int nint = el.GaussPoints();
					for (int n=0; n<nint; ++n)
					{
						double J0 = psd->detJ0(el, n);
						if (J0 <= 0)
						{
							clog.printf("**************************** E R R O R ****************************\n");
							clog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.m_nID);
							clog.printf("Jacobian = %lg\n", J0);
							clog.printf("Did you use the right node numbering?\n");
							clog.printf("Nodes:");
							for (int l=0; l<el.Nodes(); ++l)
							{
								clog.printf("%d", el.m_node[l]+1);
								if (l+1 != el.Nodes()) clog.printf(","); else clog.printf("\n");
							}
							clog.printf("*******************************************************************\n\n");
							++ninverted;
						}
					}
				}
			}
		}
	}

	// report number of inverted elements
	if (ninverted != 0)
	{
		clog.printf("**************************** E R R O R ****************************\n");
		clog.printf(" FEBio found %d initially inverted elements.\n", ninverted);
		clog.printf(" Run will be aborted.\n");
		clog.printf("*******************************************************************\n\n");
		return false;
	}

	// next if a node does not belong to a shell
	// we turn of the rotational degrees of freedom
	// To do this, we first tag all shell nodes
	vector<int> tag(m.Nodes());
	zero(tag);
	for (nd = 0; nd < m.Domains(); ++nd)
	{
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				n = el.Nodes();
				en = &el.m_node[0];
				for (j=0; j<n; ++j) tag[en[j]] = 1;
			}
		}
	}

	// fix rotational degrees of freedom of tagged nodes
	for (i=0; i<m.Nodes(); ++i) 
	{
		FENode& node = m.Node(i);
		if (tag[i] == 0)
		{
			node.m_ID[DOF_U] = -1;
			node.m_ID[DOF_V] = -1;
			node.m_ID[DOF_W] = -1;
		}
	}

	// At this point, the node ID still contains the boundary conditions
	// so we copy that into the m_BC array
	// TODO: perhaps we should put the BC's initially in BC instead of ID.
	for (i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		for (j=0; j<MAX_NDOFS; ++j)	node.m_BC[j] = node.m_ID[j];
	}

	// reset data
	m.Reset();

	// intialize domains
	for (i=0; i<m_mesh.Domains(); ++i) m_mesh.Domain(i).Initialize(*this);

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize material data
bool FEM::InitMaterials()
{
	int i;

	// initialize material data
	for (i=0; i<Materials(); ++i)
	{
		// get the material
		FEMaterial* pmat = GetMaterial(i);

		// initialize material data
		try
		{
			pmat->Init();
		}
		catch (MaterialError e)
		{
			clog.printf("Failed initializing material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			clog.printf("ERROR: %s\n\n", e.Error());
			return false;
		}
		catch (MaterialRangeError e)
		{
			clog.printf("Failed initializing material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			clog.printf("ERROR: parameter \"%s\" out of range ", e.m_szvar);
			if (e.m_bl) clog.printf("["); else clog.printf("(");
			clog.printf("%lg, %lg", e.m_vmin, e.m_vmax);
			if (e.m_br) clog.printf("]"); else clog.printf(")");
			clog.printf("\n\n");
			return false;
		}
		catch (...)
		{
			clog.printf("A fatal error occured during material intialization\n\n");
			return false;
		}
	}

	// initialize discrete materials
	try
	{
		for (i=0; i<(int) m_MAT.size(); ++i)
		{
			FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(m_MAT[i]);
			if (pm)
			{
				if (dynamic_cast<FENonLinearSpring*>(pm))
				{
					FENonLinearSpring* ps = dynamic_cast<FENonLinearSpring*>(pm);
					ps->m_plc = GetLoadCurve(ps->m_nlc);
				}
				pm->Init();
			}
		}
	}
	catch (...)
	{
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function resets the FEM data so that a new run can be done.
//! This routine is called from the optimization routine.

bool FEM::Reset()
{
	int i;

	// initialize materials
	FEMaterial* pmat;
	
	for (i=0; i<Materials(); ++i)
	{
		pmat = GetMaterial(i);
		pmat->Init();
	}

	// reset mesh data
	m_mesh.Reset();

	// reset object data
	int nrb = m_Obj.size();
	for (i=0; i<nrb; ++i) m_Obj[i]->Reset();

	// set up rigid joints
	if (!m_NLC.empty())
	{
		int NC = (int) m_NLC.size();
		for (i=0; i<NC; ++i)
		{
			FENLConstraint* plc = m_NLC[i];
			if (plc->Type() == FE_RIGID_JOINT)
			{
				FERigidJoint& rj = dynamic_cast<FERigidJoint&>(*plc);
				rj.m_F = vec3d(0,0,0);

				FERigidBody& ra = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBa]);
				FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBb]);

				rj.m_qa0 = rj.m_q0 - ra.m_r0;
				rj.m_qb0 = rj.m_q0 - rb.m_r0;
			}
		}
	}

	// set the start time
	m_ftime = 0;
	m_ftime0 = 0;

	// set first time step
	m_pStep->m_dt = m_pStep->m_dt0;
	m_pStep->m_ntotref    = 0;		// total nr of stiffness reformations
	m_pStep->m_ntotiter   = 0;		// total nr of non-linear iterations
	m_pStep->m_ntimesteps = 0;		// time steps completed
	m_pStep->m_ntotrhs    = 0;		// total nr of right hand side evaluations

	// open plot database file
	if (m_pStep->m_nplot != FE_PLOT_NEVER)
	{
		if (m_plot == 0) m_plot = new LSDYNAPlotFile;

		if (m_plot->Open(*this, m_szplot) == false)
		{
			clog.printf("ERROR : Failed creating PLOT database\n");
			return false;
		}
	}

	// Since it is assumed that for the first timestep
	// there are no loads or initial displacements, the case n=0 is skipped.
	// Therefor we can output those results here.
	// Offcourse we should actually check if this is indeed
	// the case, otherwise we should also solve for t=0
	if (m_pStep->m_nplot != FE_PLOT_NEVER) m_plot->Write(*this);
/*
	// reset the log file
	if (!log.is_valid())
	{
		log.open(m_szlog);

		// if we don't want to output anything we only output to the logfile
		if (m_pStep->GetPrintLevel() == FE_PRINT_NEVER) log.SetMode(Logfile::FILE_ONLY);

		// print welcome message to file
		Hello(log);
	}
*/
	// do the callback
	DoCallback();

	// All data is reset successfully
	return true;
}

//-----------------------------------------------------------------------------
//! Initializes contact data
// TODO: I should probably let the analysis step take care of this
bool FEM::InitContact()
{
	// loop over all contact interfaces
	for (int i=0; i<ContactInterfaces(); ++i)
	{
		// get the contact interface
		FEContactInterface& ci = *m_CI[i];

		// initializes contact interface data
		ci.Init();
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize solute-poroelastic data.
//! Find all nodes that are not part of a poro-solute domain and fix the 
//! pressure and concentration DOFS. 
//! \todo This function should probably move to the FEAnalysisStep class.
bool FEM::InitPoroSolute()
{
	int i, j, nd;

	// make sure this is the poro-solute module
	int nstep = m_pStep->GetType();
	bool bporo = (nstep == FE_BIPHASIC  ) || (nstep == FE_POROSOLUTE) || (nstep == FE_TRIPHASIC);
	bool bsolu = (nstep == FE_POROSOLUTE) || (nstep == FE_TRIPHASIC );
	bool btri  = (nstep == FE_TRIPHASIC );
	
	if (!bporo)
	{
		// if there is no poroelasticity
		// we set all pressure degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[DOF_P] = -1;
		
		// also remove prescribed pressures
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& bc   = m_DC[i]->bc;
			if (bc == DOF_P) bc = -1;
		}
	}
	
	if (!bsolu)
	{
		// if there is no solute
		// we set all concentration degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[DOF_C] = -1;
		
		// also remove prescribed concentrations
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& bc   = m_DC[i]->bc;
			if (bc == DOF_C) bc = -1;
		}
	}
	
	if (!btri)
	{
		// if there is no solute
		// we set all concentration degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[DOF_C+1] = -1;
		
		// also remove prescribed concentrations
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& bc   = m_DC[i]->bc;
			if (bc == int(DOF_C+1)) bc = -1;
		}
	}
	
	if (!bporo)
		// let's go back
		return true;
	
	if (bporo)
	{
		// fix all pressure dofs that are not used
		// that is, that are not part of a biphasic, biphasic-solute, or triphasic element
		// this is done in three steps
		// step 1. mark all biphasic nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(pbd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (bm || bsm || btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_P] == 0) m_mesh.Node(n[j]).m_ID[DOF_P] = 1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(psd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (bm || bsm || btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_P] == 0) m_mesh.Node(n[j]).m_ID[DOF_P] = 1;
					}
				}
			}
		}
		
		// step 2. fix pressure dofs of all unmarked nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(pbd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if ((bm == 0) && (bsm == 0) && (btm == 0))
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_P] != 1) m_mesh.Node(n[j]).m_ID[DOF_P] = -1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(psd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if ((bm == 0) && (bsm == 0) && (btm == 0))
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_P] != 1) m_mesh.Node(n[j]).m_ID[DOF_P] = -1;
					}
				}
			}
		}
		
		// step 3. free all marked dofs
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(pbd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (bm || bsm || btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_P] == 1) m_mesh.Node(n[j]).m_ID[DOF_P] = 0;
					}
				}
			}
			
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(psd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (bm || bsm || btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_P] == 1) m_mesh.Node(n[j]).m_ID[DOF_P] = 0;
					}
				}
			}
		}
	}
	
	if (bsolu)
	{
		// fix all neutral/cation concentration dofs that are not used
		// that is, that are not part of a solute-solid element
		// this is done in three steps
		// step 1. mark all solute-solid nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (bsm || btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_C] == 0) m_mesh.Node(n[j]).m_ID[DOF_C] = 1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (bsm || btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_C] == 0) m_mesh.Node(n[j]).m_ID[DOF_C] = 1;
					}
				}
			}
		}
		
		// step 2. fix concentration dofs of all unmarked nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if ((bsm == 0) && (btm == 0))
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C] != 1) m_mesh.Node(n[j]).m_ID[DOF_C] = -1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if ((bsm == 0) && (btm == 0))
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C] != 1) m_mesh.Node(n[j]).m_ID[DOF_C] = -1;
					}
				}
			}
		}
		
		// step 3. free all marked dofs
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (bsm || btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C] == 1) m_mesh.Node(n[j]).m_ID[DOF_C] = 0;
					}
				}
			}
			
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (bsm || btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C] == 1) m_mesh.Node(n[j]).m_ID[DOF_C] = 0;
					}
				}
			}
		}
	}

	if (btri)
	{
		// fix all anion concentration dofs that are not used
		// that is, that are not part of a solute-solid element
		// this is done in three steps
		// step 1. mark all solute-solid nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] == 0) m_mesh.Node(n[j]).m_ID[DOF_C+1] = 1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] == 0) m_mesh.Node(n[j]).m_ID[DOF_C+1] = 1;
					}
				}
			}
		}
		
		// step 2. fix concentration dofs of all unmarked nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm == 0)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] != 1) m_mesh.Node(n[j]).m_ID[DOF_C+1] = -1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm == 0)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] != 1) m_mesh.Node(n[j]).m_ID[DOF_C+1] = -1;
					}
				}
			}
		}
		
		// step 3. free all marked dofs
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] == 1) m_mesh.Node(n[j]).m_ID[DOF_C+1] = 0;
					}
				}
			}
			
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] == 1) m_mesh.Node(n[j]).m_ID[DOF_C+1] = 0;
					}
				}
			}
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! This function creates the rigid bodies by analyzing the rigid materials
//! and the mesh in the model. 
//!
bool FEM::CreateRigidBodies()
{
	int i, j, n, m, nd;
	// count the number of rigid materials
	int nrm = 0;
	for (i=0; i<Materials(); ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(i));
		if (pm) nrm++;
	}
	
	// make sure there are rigid materials
	if (nrm == 0) return true;

	// First we need to figure out how many rigid bodies there are.
	// This is not the same as rigid materials, since a rigid body
	// may be composed of different rigid materials (similarly to a deformable
	// body that may contain different materials). Although there can
	// only be one deformable mesh, there can be several rigid bodies.

	// The mrb array will contain an index to the rigid body the material
	// is attached too.
	vector<int> mrb(Materials());
	n = 0;
	for (i=0; i<Materials(); ++i)
	{
		if (dynamic_cast<FERigidMaterial*>(GetMaterial(i))) mrb[i] = n++;
		else mrb[i] = -1;
	}

	// Next, we assign to all nodes a rigid node number
	// This number is preliminary since rigid materials can be merged
	for (nd = 0; nd < m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					el.m_nrigid = el.GetMatID();
					for (j=0; j<el.Nodes(); ++j)
					{
						n = el.m_node[j];
						FENode& node = m_mesh.Node(n);
						node.m_rid = el.GetMatID();
					}
				}
				else el.m_nrigid = -1;
			}
		}

		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					el.m_nrigid = el.GetMatID();
					for (j=0; j<el.Nodes(); ++j)
					{
						n = el.m_node[j];
						FENode& node = m_mesh.Node(n);
						node.m_rid = el.GetMatID();
					}		
				}
				else el.m_nrigid = -1;
			}
		}
	}

	// now we can merge rigid materials
	// if a rigid element has two nodes that connect to two different
	// rigid materials we need to merge. 
	bool bdone;
	do
	{
		bdone = true;
		for (nd=0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (el.m_nrigid >= 0)
					{
						m = m_mesh.Node(el.m_node[0]).m_rid;
						for (j=1; j<el.Nodes(); ++j)
						{
							n = m_mesh.Node(el.m_node[j]).m_rid;
							if (mrb[n] != mrb[m])
							{
								if (mrb[n]<mrb[m]) mrb[m] = mrb[n]; else mrb[n] = mrb[m];
								bdone = false;
							}
						}
					}
				}
			}

			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (el.m_nrigid >= 0)
					{
						m = m_mesh.Node(el.m_node[0]).m_rid;
						for (j=1; j<el.Nodes(); ++j)
						{
							n = m_mesh.Node(el.m_node[j]).m_rid;
							if (mrb[n] != mrb[m])
							{
								if (mrb[n]<mrb[m]) mrb[m] = mrb[n]; else mrb[n] = mrb[m];
								bdone = false;
							}
						}
					}
				}
			}
		}
	}
	while (!bdone);

	// since we may have lost a rigid body in the merge process
	// we reindex the RB's.
	int nmat = Materials();
	vector<int> mrc; mrc.assign(nmat, -1);
	for (i=0; i<nmat; ++i) if (mrb[i] >= 0) mrc[mrb[i]] = 0;
	int nrb = 0;
	for (i=0; i<nmat; ++i)
	{
		if (mrc[i] == 0) mrc[i] = nrb++;
	}

	for (i=0; i<nmat; ++i) 
	{
		if (mrb[i] >= 0) mrb[i] = mrc[mrb[i]];
	}

	// set rigid body index for materials
	for (i=0; i<Materials(); ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(i));
		if (pm)	
		{
			pm->m_nRB = mrb[i];
		}
	}

	// assign rigid body index to rigid elements
	for (nd=0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(el.GetMatID()));
				if (pm)
					el.m_nrigid = pm->m_nRB;
				else
					el.m_nrigid = -1;
			}
		}

		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(el.GetMatID()));
				if (pm)
					el.m_nrigid = pm->m_nRB;
				else
					el.m_nrigid = -1;
			}
		}
	}

	// assign rigid body index to nodes
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0) node.m_rid = mrb[ node.m_rid ];
	}

	// Ok, we now know how many rigid bodies there are
	// so let's create them
	m_Obj.clear();
	for (i=0; i<nrb; ++i)
	{
		// create a new rigid body
		FERigidBody* prb = new FERigidBody(this);
		prb->m_nID = i;

		// Since a rigid body may contain several rigid materials
		// we find the first material that this body has and use
		// that materials data to set up the rigid body data
		FERigidMaterial* pm = 0;
		for (j=0; j<Materials(); ++j)
		{
			pm = dynamic_cast<FERigidMaterial*> (GetMaterial(j));

			if (pm && (pm->m_nRB == i))	break;
		}
		assert(j<Materials());
		prb->m_mat = j;

		// initialize center of mass
		if (pm->m_com == 1)
		{
			// grab the com from the material
			prb->m_r0 = prb->m_rt = pm->m_rc;
		}
		else
		{
			// calculate the com
			prb->UpdateCOM();
		}

		// add it to the pile
		m_Obj.push_back(prb);
	}

	// set up rigid joints
	if (!m_NLC.empty())
	{
		FERigidMaterial* pm;
		int NC = m_NLC.size();
		for (i=0; i<NC; ++i)
		{
			FENLConstraint* plc = m_NLC[i];
			if (plc->Type() == FE_RIGID_JOINT)
			{
				FERigidJoint& rj = dynamic_cast<FERigidJoint&>(*plc);
				rj.m_F = vec3d(0,0,0);

				pm = dynamic_cast<FERigidMaterial*> (GetMaterial(rj.m_nRBa));
				if (pm == 0)
				{
					clog.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", i+1);
					return false;
				}
				rj.m_nRBa = pm->m_nRB;

				pm = dynamic_cast<FERigidMaterial*> (GetMaterial(rj.m_nRBb));
				if (pm == 0)
				{
					clog.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", i+1);
					return false;
				}
				rj.m_nRBb = pm->m_nRB;

				FERigidBody& ra = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBa]);
				FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBb]);

				rj.m_qa0 = rj.m_q0 - ra.m_r0;
				rj.m_qb0 = rj.m_q0 - rb.m_r0;
			}
		}
	}

	// overwrite rigid nodes degrees of freedom
	// We do this so that these dofs do not
	// get equation numbers assigned to them. Later we'll assign
	// the rigid dofs equations numbers to these nodes
	for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_bshell = false;
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				if (el.m_nrigid < 0)
				{
					int n = el.Nodes();
					for (j=0; j<n; ++j) m_mesh.Node(el.m_node[j]).m_bshell = true;
				}
			}
		}
	}

	// The following fixes the degrees of freedom for rigid nodes.
	// Note that also the rotational degrees of freedom are fixed
	// for rigid nodes that do not belong to a non-rigid shell element.
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0)
		{
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
	}

	// assign correct rigid body ID's to rigid nodes
	for (i=0; i<(int) m_RN.size(); ++i)
	{
		FERigidNode& rn = *m_RN[i];
		rn.rid = mrb[rn.rid];
	}

	// let's find all rigid surface elements
	// a surface element is rigid when it has no free nodes
	for (int is = 0; is < (int) m_SL.size(); ++is)
	{
		FESurfaceLoad* ps = m_SL[is];
		for (i=0; i<ps->Surface().Elements(); ++i)
		{
			FESurfaceElement& el = ps->Surface().Element(i);
			int N = el.Nodes();
			el.m_nrigid = 0;
			for (j=0; j<N; ++j) 
			{
				FENode& node = m_mesh.Node(el.m_node[j]);
				if (node.m_rid < 0) 
				{
					el.m_nrigid = -1;
					break;
				}
			}
		}
	}

	// the rigid body constraints are still associated with the rigid materials
	// so we now associate them with the rigid bodies
	for (i=0; i<(int) m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *m_RDC[i];
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(DC.id-1));
		DC.id = pm->m_nRB;
	}
	for (i=0; i<(int) m_RFC.size(); ++i)
	{
		FERigidBodyForce& FC = *m_RFC[i];
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(FC.id-1));
		FC.id = pm->m_nRB;
	}

	// set the rigid body parents
	for (i=0; i<(int) m_Obj.size(); ++i)
	{
		FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[i]);
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_MAT[rb.m_mat]);
		assert(pm);
		if (pm->m_pmid > -1)
		{
			FERigidMaterial* ppm = dynamic_cast<FERigidMaterial*>(m_MAT[pm->m_pmid-1]);
			assert(ppm);
			FERigidBody& prb = dynamic_cast<FERigidBody&>(*m_Obj[ppm->m_nRB]);
			rb.m_prb = &prb;

			// we also need to open up all the RB's degree of freedoms
			pm->m_bc[0] = 1;
			pm->m_bc[1] = 1;
			pm->m_bc[2] = 1;
			pm->m_bc[3] = 1;
			pm->m_bc[4] = 1;
			pm->m_bc[5] = 1;
		}
	}

	return true;
}
