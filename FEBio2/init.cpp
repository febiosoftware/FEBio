#include "stdafx.h"
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include "fem.h"
#include "FECore/FEException.h"
#include "FEBioLib/log.h"
#include "FESolidSolver.h"
#include "LSDYNAPlotFile.h"
#include "FEBioLib/FETransverselyIsotropic.h"
#include "FEBioLib/FEDiscreteMaterial.h"
#include "FEBioLib/FEElasticSolidDomain.h"
#include "FEBioLib/FEElasticShellDomain.h"
#include "FEBioLib/FEPointBodyForce.h"

// Forward declarations
void Hello(FILE* fp);

bool err(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// print to file
	va_start(args, sz);
	vfprintf(stderr, sz, args);
	va_end(args);

	return false;
}

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
		if ((step.m_ntime <= 0) && (step.m_final_time <= 0.0)) return err("Invalid number of time steps for analysis step %d", i+1);
		if ((step.m_ntime >  0) && (step.m_final_time >  0.0)) return err("You must either set the number of time steps or the final time but not both.\n");
		if (step.m_dt0   <= 0) return err("Invalid time step size for analysis step %d", i+1);
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
						fprintf(stderr, "**************************** E R R O R ****************************\n");
						fprintf(stderr, "Negative jacobian detected at integration point %d of element %d\n", n+1, el.m_nID);
						fprintf(stderr, "Jacobian = %lg\n", J0);
						fprintf(stderr, "Did you use the right node numbering?\n");
						fprintf(stderr, "Nodes:");
						for (int l=0; l<el.Nodes(); ++l)
						{
							fprintf(stderr, "%d", el.m_node[l]+1);
							if (l+1 != el.Nodes()) fprintf(stderr, ","); else fprintf(stderr, "\n");
						}
						fprintf(stderr, "*******************************************************************\n\n");
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
							fprintf(stderr, "**************************** E R R O R ****************************\n");
							fprintf(stderr, "Negative jacobian detected at integration point %d of element %d\n", n+1, el.m_nID);
							fprintf(stderr, "Jacobian = %lg\n", J0);
							fprintf(stderr, "Did you use the right node numbering?\n");
							fprintf(stderr, "Nodes:");
							for (int l=0; l<el.Nodes(); ++l)
							{
								fprintf(stderr, "%d", el.m_node[l]+1);
								if (l+1 != el.Nodes()) fprintf(stderr, ","); else fprintf(stderr, "\n");
							}
							fprintf(stderr, "*******************************************************************\n\n");
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
		fprintf(stderr, "**************************** E R R O R ****************************\n");
		fprintf(stderr, " FEBio found %d initially inverted elements.\n", ninverted);
		fprintf(stderr, " Run will be aborted.\n");
		fprintf(stderr, "*******************************************************************\n\n");
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
//! This function initializes the linear constraint table (LCT). This table
//! contains for each dof the linear constraint it belongs to. (or -1 if it is
//! not constraint)

bool FEM::InitConstraints()
{
	int nlin = m_LinC.size();
	if (nlin == 0) return true;

	int i;

	// set the equation numbers for the linear constraints
	list<FELinearConstraint>::iterator it = m_LinC.begin();
	for (i=0; i<nlin; ++i, ++it)
	{
		FELinearConstraint& lc = *it;
		lc.master.neq = m_mesh.Node(lc.master.node).m_ID[lc.master.bc];

		// make sure the master did not get assigned an equation
		assert(lc.master.neq == -1);

		// set the slave equation numbers
		list<FELinearConstraint::SlaveDOF>::iterator is = lc.slave.begin();
		int nn = lc.slave.size();
		for (int n=0; n<nn; ++n, ++is)
		{
			FELinearConstraint::SlaveDOF& sn = *is;
			sn.neq = m_mesh.Node(sn.node).m_ID[sn.bc];
		}		
	}

	// create the linear constraint table
	m_LCT.assign(m_mesh.Nodes()*MAX_NDOFS, -1);

	list<FELinearConstraint>::iterator ic = m_LinC.begin();
	for (i=0; i<nlin; ++i, ++ic)
	{
		FELinearConstraint& lc = *ic;
		int n = lc.master.node;
		int m = lc.master.bc;
		
		m_LCT[n*MAX_NDOFS+m] = i;
	}

	// to simplify accessing the linear constraint data
	// we store all pointers in an array
	// TODO: perhaps I should store the linear constraints that way
	// anyways and get rid of the list
	m_LCA.resize(nlin);
	ic = m_LinC.begin();
	for (i=0; i<nlin; ++i, ++ic) m_LCA[i] = &(*ic);

	// let's do the aug lag linear constraints
	// TODO: This is also done in FEM::Init and FEAnalysis::Init. Where do I really need to do this?
	if (m_NLC.size() > 0)
	{
		int N = m_NLC.size();
		for (i=0; i<N; ++i) m_NLC[i]->Init();
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
