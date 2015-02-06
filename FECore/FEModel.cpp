#include "stdafx.h"
#include "FEModel.h"
#include "FEShellDomain.h"
#include "FERigidBody.h"
#include <string>
#include "log.h"
#include "FECoreKernel.h"
using namespace std;

//-----------------------------------------------------------------------------
FEModel::FEModel(void)
{
	m_sztitle[0] = 0;

	// --- Analysis Data ---
	m_pStep = 0;
	m_nStep = -1;
	m_nplane_strain = -1;	// don't use plain strain mode
	m_ftime = 0;
	m_ftime0 = 0;
	m_bwopt = 0;

	// additional data
	m_ut4_alpha = 0.05;
	m_ut4_bdev = false;
	m_udghex_hg = 1.0;
}

//-----------------------------------------------------------------------------
//! Delete all dynamically allocated data
FEModel::~FEModel(void)
{
	Clear();
}

//-----------------------------------------------------------------------------
void FEModel::Clear()
{
	size_t i;
	for (i=0; i<m_Step.size(); ++i) delete m_Step[i]; m_Step.clear();
	for (i=0; i<m_CI.size  (); ++i) delete m_CI [i] ; m_CI.clear  ();
	for (i=0; i<m_MAT.size (); ++i) delete m_MAT[i] ; m_MAT.clear ();
	for (i=0; i<m_LC.size  (); ++i) delete m_LC [i] ; m_LC.clear  ();
	for (i=0; i<m_BL.size  (); ++i) delete m_BL [i] ; m_BL.clear  ();
	for (i=0; i<m_DC.size  (); ++i) delete m_DC [i] ; m_DC.clear  ();
	for (i=0; i<m_FC.size  (); ++i) delete m_FC [i] ; m_FC.clear  ();
	for (i=0; i<m_SL.size  (); ++i) delete m_SL [i] ; m_SL.clear  ();
	for (i=0; i<m_RDC.size (); ++i) delete m_RDC[i] ; m_RDC.clear ();
	for (i=0; i<m_RFC.size (); ++i) delete m_RFC[i] ; m_RFC.clear ();
	for (i=0; i<m_RBV.size (); ++i) delete m_RBV[i] ; m_RBV.clear ();
	for (i=0; i<m_RBW.size (); ++i) delete m_RBW[i] ; m_RBW.clear ();
	for (i=0; i<m_RN.size  (); ++i) delete m_RN [i] ; m_RN.clear  ();
	for (i=0; i<m_NLC.size (); ++i) delete m_NLC[i] ; m_NLC.clear ();
	for (i=0; i<m_Obj.size (); ++i) delete m_Obj[i] ; m_Obj.clear ();

	// global data
	for (i=0; i<m_GD.size(); ++i) delete m_GD[i]; m_GD.clear();
	m_Const.clear();
}

//-----------------------------------------------------------------------------
void FEModel::ClearBCs()
{
	for (size_t i=0; i<m_DC.size  (); ++i) delete m_DC[i];
	m_DC.clear();
}

//-----------------------------------------------------------------------------
bool FEModel::Init()
{
	// intitialize time
	m_ftime = 0;
	m_ftime0 = 0;

	// check step data
	// TODO: should I let the Steps take care of this instead?
	for (int i=0; i<(int) m_Step.size(); ++i)
	{
		FEAnalysis& step = *m_Step[i];
		if ((step.m_ntime <= 0) && (step.m_final_time <= 0.0)) { felog.printf("Invalid number of time steps for analysis step %d", i+1); return false; }
		if ((step.m_ntime >  0) && (step.m_final_time >  0.0)) { felog.printf("You must either set the number of time steps or the final time but not both.\n"); return false; }
		if (step.m_dt0   <= 0) { felog.printf("Invalid time step size for analysis step %d", i+1); return false; }
		if (step.m_bautostep)
		{
//			if (m_pStep->m_dtmin <= 0) return err("Invalid minimum time step size");
//			if (m_pStep->m_dtmax <= 0) return err("Invalid maximum time step size");
		}
	}

	// evaluate all loadcurves at the initial time
	for (int i=0; i<LoadCurves(); ++i) m_LC[i]->Evaluate(0);
    
    // evaluate all parameter lists
    EvaluateAllParameterLists();

	// if the analysis is run in plane-strain mode we fix all the z-dofs of all nodes
	if (m_nplane_strain >= 0)
	{
		int bc = m_nplane_strain;
		for (int i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_BC[bc] = -1;
	}

	// validate BC's
	if (InitBCs() == false) return false;

	// create and initialize the rigid body data
	if (InitObjects() == false) return false;

	// initialize material data
	// NOTE: call this before InitMesh since we need to initialize the FECoordSysMap
	// before we can calculate the local element coordinate systems.
	if (InitMaterials() == false) return false;

	// initialize mesh data
	// NOTE: this must be done AFTER the elements have been assigned material point data !
	// this is because the mesh data is reset
	// TODO: perhaps I should not reset the mesh data during the initialization
	if (InitMesh() == false) return false;

	// initialize contact data
	if (InitContact() == false) return false;

	// init body loads
	if (InitBodyLoads() == false) return false;

	// initialize nonlinear constraints
	// TODO: This is also initialized in the analysis step. Do I need to do this here?
	if (InitConstraints() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! See if the BC's are setup correctly.
bool FEModel::InitBCs()
{
	// get the number of loadcurves
	int NLC = LoadCurves();

	// check the prescribed BC's
	int NBC = PrescribedBCs();
	for (int i=0; i<NBC; ++i)
	{
		FEPrescribedBC* pbc = PrescribedBC(i);
		if ((pbc->lc < 0)||(pbc->lc >= NLC))
		{
			felog.printf("ERROR: Invalid loadcurve in prescribed BC %d\n", i+1);
			return false;
		}
	}

	// check the nodal loads
	int NNL = NodalLoads();
	for (int i=0; i<NNL; ++i)
	{
		FENodalForce* pbc = NodalLoad(i);
		if ((pbc->lc < 0)||(pbc->lc >= NLC))
		{
			felog.printf("ERROR: Invalid loadcurve in nodal load %d\n", i+1);
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Update the mesh data. This function calculates the initial directors
//! for the shell elements.

// NOTE: This function needs to be called after the rigid bodies have been
// initialized

bool FEModel::InitMesh()
{
	// Initialize mesh
	if (m_mesh.Init() == false) return false;

	// intialize domain data
	// TODO: I'd like the mesh to take care of this, but I need to pass a FEModel pointer for some reason.
	for (int i=0; i<m_mesh.Domains(); ++i) m_mesh.Domain(i).Initialize(*this);

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize material data
bool FEModel::InitMaterials()
{
	// initialize material data
	for (int i=0; i<Materials(); ++i)
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
			felog.printf("Failed initializing material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			felog.printf("ERROR: %s\n\n", e.Error());
			return false;
		}
		catch (MaterialRangeError e)
		{
			felog.printf("Failed initializing material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			felog.printf("ERROR: parameter \"%s\" out of range ", e.m_szvar);
			if (e.m_bl) felog.printf("["); else felog.printf("(");
			felog.printf("%lg, %lg", e.m_vmin, e.m_vmax);
			if (e.m_br) felog.printf("]"); else felog.printf(")");
			felog.printf("\n\n");
			return false;
		}
		catch (...)
		{
			felog.printf("A fatal error occured during material intialization\n\n");
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function creates and initializes the rigid bodies, the only non-deformable
//! "objects" supported so far.
bool FEModel::InitObjects()
{
	// count the number of rigid materials
	int nrm = 0;
	for (int i=0; i<Materials(); ++i)
	{
		if (GetMaterial(i)->IsRigid()) nrm++;
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
	int n = 0;
	for (int i=0; i<Materials(); ++i)
	{
		if (GetMaterial(i)->IsRigid()) mrb[i] = n++;
		else mrb[i] = -1;
	}

	// Next, we assign to all nodes a rigid node number
	// This number is preliminary since rigid materials can be merged
	for (int nd = 0; nd < m_mesh.Domains(); ++nd)
	{
		FEDomain& dom = m_mesh.Domain(nd);
		FEMaterial* pmat = dom.GetMaterial();
		if (pmat->IsRigid())
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				el.m_nrigid = pmat->GetID() - 1;
				for (int j=0; j<el.Nodes(); ++j)
				{
					int n = el.m_node[j];
					FENode& node = m_mesh.Node(n);
					node.m_rid = pmat->GetID() - 1;
				}
			}
		}
		else 
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				el.m_nrigid = -1;
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
		for (int nd=0; nd<m_mesh.Domains(); ++nd)
		{
			FEDomain& dom = m_mesh.Domain(nd);
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				if (el.m_nrigid >= 0)
				{
					int m = m_mesh.Node(el.m_node[0]).m_rid;
					for (int j=1; j<el.Nodes(); ++j)
					{
						int n = m_mesh.Node(el.m_node[j]).m_rid;
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
	while (!bdone);

	// since we may have lost a rigid body in the merge process
	// we reindex the RB's.
	int nmat = Materials();
	vector<int> mrc; mrc.assign(nmat, -1);
	for (int i=0; i<nmat; ++i) if (mrb[i] >= 0) mrc[mrb[i]] = 0;
	int nrb = 0;
	for (int i=0; i<nmat; ++i)
	{
		if (mrc[i] == 0) mrc[i] = nrb++;
	}

	for (int i=0; i<nmat; ++i) 
	{
		if (mrb[i] >= 0) mrb[i] = mrc[mrb[i]];
	}

	// set rigid body index for materials
	for (int i=0; i<Materials(); ++i)
	{
		FEMaterial* pm = GetMaterial(i);
		if (pm->IsRigid())	
		{
			pm->SetRigidBodyID(mrb[i]);
		}
	}

	// assign rigid body index to rigid elements
	for (int nd=0; nd<m_mesh.Domains(); ++nd)
	{
		FEDomain& dom = m_mesh.Domain(nd);
		FEMaterial* pm = dom.GetMaterial();
		for (int i=0; i<dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			if (pm->IsRigid())
				el.m_nrigid = pm->GetRigidBodyID();
			else
				el.m_nrigid = -1;
		}
	}

	// assign rigid body index to nodes
	for (int i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0) node.m_rid = mrb[ node.m_rid ];
	}

	// Ok, we now know how many rigid bodies there are
	// so let's create them
	m_Obj.clear();
	for (int i=0; i<nrb; ++i)
	{
		// create a new rigid body
		FERigidBody* prb = new FERigidBody(this);
		prb->m_nID = i;

		// Since a rigid body may contain several rigid materials
		// we find the first material that this body has and use
		// that materials data to set up the rigid body data
		int j;
		FEMaterial* pm = 0;
		for (j=0; j<Materials(); ++j)
		{
			pm = GetMaterial(j);
			if (pm && (pm->GetRigidBodyID() == i))	break;
		}
		assert(j<Materials());
		prb->m_mat = j;

		// add it to the pile
		m_Obj.push_back(prb);
	}

	// overwrite rigid nodes degrees of freedom
	// We do this so that these dofs do not
	// get equation numbers assigned to them. Later we'll assign
	// the rigid dofs equations numbers to these nodes
	for (int i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_bshell = false;
	for (int nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEDomain& dom = m_mesh.Domain(nd);
		if (dom.Class() == FE_DOMAIN_SHELL)
		{
			int N = dom.Elements();
			for (int i=0; i<N; ++i)
			{
				FEElement& el = dom.ElementRef(i);
				if (el.m_nrigid < 0)
				{
					int n = el.Nodes();
					for (int j=0; j<n; ++j) m_mesh.Node(el.m_node[j]).m_bshell = true;
				}
			}
		}
	}

	// The following fixes the degrees of freedom for rigid nodes.
	// Note that also the rotational degrees of freedom are fixed
	// for rigid nodes that do not belong to a non-rigid shell element.
	for (int i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0)
		{
			node.m_BC[DOF_X] = -1;
			node.m_BC[DOF_Y] = -1;
			node.m_BC[DOF_Z] = -1;
			if (node.m_bshell == false)
			{
				node.m_BC[DOF_U] = -1;
				node.m_BC[DOF_V] = -1;
				node.m_BC[DOF_W] = -1;
			}
		}
	}

	// assign correct rigid body ID's to rigid nodes
	for (int i=0; i<(int) m_RN.size(); ++i)
	{
		FERigidNode& rn = *m_RN[i];
		rn.rid = mrb[rn.rid];
	}

	// let's find all rigid surface elements
	// a surface element is rigid when it has no free nodes
	for (int is = 0; is < (int) m_SL.size(); ++is)
	{
		FESurfaceLoad* ps = m_SL[is];
		for (int i=0; i<ps->Surface().Elements(); ++i)
		{
			FESurfaceElement& el = ps->Surface().Element(i);
			int N = el.Nodes();
			el.m_nrigid = 0;
			for (int j=0; j<N; ++j) 
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
	for (int i=0; i<(int) m_RBC.size(); ++i)
	{
		FERigidBodyFixedBC& BC = *m_RBC[i];
		FEMaterial* pm = GetMaterial(BC.id-1);
		BC.id = pm->GetRigidBodyID(); assert(BC.id >= 0);
	}
	for (int i=0; i<(int) m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *m_RDC[i];
		FEMaterial* pm = GetMaterial(DC.id-1);
		DC.id = pm->GetRigidBodyID(); assert(DC.id >= 0);
	}
	for (int i=0; i<(int) m_RFC.size(); ++i)
	{
		FERigidBodyForce& FC = *m_RFC[i];
		FEMaterial* pm = GetMaterial(FC.id-1);
		FC.id = pm->GetRigidBodyID(); assert(FC.id >= 0);
	}
	for (int i=0; i<(int) m_RBV.size(); ++i)
	{
		FERigidBodyVelocity& RV = *m_RBV[i];
		FEMaterial* pm = GetMaterial(RV.id-1);
		RV.id = pm->GetRigidBodyID(); assert(RV.id>=0);
	}
	for (int i=0; i<(int) m_RBW.size(); ++i)
	{
		FERigidBodyAngularVelocity& RW = *m_RBW[i];
		FEMaterial* pm = GetMaterial(RW.id-1);
		RW.id = pm->GetRigidBodyID(); assert(RW.id>=0);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Initializes contact data
//! \todo Contact interfaces have two initialization functions: Init() and
//!   Activate(). The Init member is called here and allocates the required memory
//!   for each interface. The Activate() member is called during the step initialization
//!   which is called later during the solve phase. However, for global interfaces (i.e.
//!   interfaces that are active during the entire simulation), the Activate() member is
//!   not called in the solve phase. That is why we have to call it here. Global interfaces
//!   can be idenfitifed since they are active during the initialization. 
//!   I am not entirely a fan of this approach but it does solve the problem that contact
//!   interface shoulds only do work (e.g. update projection status) when they are active, but
//!   have to allocate memory during the initialization fase.
//!
bool FEModel::InitContact()
{
	// loop over all contact interfaces
	for (int i=0; i<SurfacePairInteractions(); ++i)
	{
		// get the contact interface
		FESurfacePairInteraction& ci = *m_CI[i];

		// initializes contact interface data
		if (ci.Init() == false) return false;

		// If the contact interface is active
		// we have to call the Activate() member. 
		if (ci.IsActive()) ci.Activate();
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEModel::InitConstraints()
{
	for (int i=0; i<(int) m_NLC.size(); ++i)
	{
		if (m_NLC[i]->Init() == false) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEModel::InitBodyLoads()
{
	for (int i=0; i<(int) m_BL.size(); ++i)
	{
		if (m_BL[i]->Init() == false) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
//! This function solves the FE problem by calling the solve method for each step.
bool FEModel::Solve()
{
	// convergence flag
	bool bconv = true;

	// loop over all analysis steps
	// Note that we don't necessarily from step 0.
	// This is because the user could have restarted
	// the analysis. 
	for (size_t nstep=m_nStep; nstep < m_Step.size(); ++nstep)
	{
		// set the current analysis step
		m_nStep = nstep;
		m_pStep = m_Step[nstep];

		// intitialize step data
		if (m_pStep->Init() == false)
		{
			bconv = false;
			break;
		}

		// solve the analaysis step
		bconv = m_pStep->Solve();

		// break if the step has failed
		if (bconv == false) break;

		// wrap it up
		m_pStep->Finish();
	}

	return bconv;
}

//-----------------------------------------------------------------------------
bool FEModel::Reset()
{
	// initialize materials
	for (int i=0; i<Materials(); ++i)
	{
		FEMaterial* pmat = GetMaterial(i);
		pmat->Init();
	}

	// reset mesh data
	m_mesh.Reset();

	// reset object data
	int nrb = m_Obj.size();
	for (int i=0; i<nrb; ++i) m_Obj[i]->Reset();

	// set up rigid joints
	if (!m_NLC.empty())
	{
		int NC = (int) m_NLC.size();
		for (int i=0; i<NC; ++i)
		{
			FENLConstraint* plc = m_NLC[i];
			plc->Reset();
		}
	}

	// set the start time
	m_ftime = 0;
	m_ftime0 = 0;

	// set first time step
	m_pStep = m_Step[0];
	m_nStep = 0;
	for (int i=0; i<(int)m_Step.size(); ++i) m_Step[i]->Reset();

	// reset contact data
	// TODO: I just call Init which I think is okay
	InitContact();

	return true;
}

//-----------------------------------------------------------------------------
//! Get the current time information.
FETimePoint FEModel::GetCurrentTime()
{
	FETimePoint pt;
	pt.t = m_ftime;
	pt.dt = GetCurrentStep()->m_dt;
	return pt;
}

//-----------------------------------------------------------------------------
//! This function is used when pushing the FEM state data. Since we don't need
//! to copy all the data, this function only copies the data that needs to be 
//! restored for a running restart.
//!
//! \todo Shallow copy nonlinear constraints
void FEModel::ShallowCopy(DumpStream& dmp, bool bsave)
{
	// stream model data
	if (bsave)
	{
		dmp << m_ftime;
	}
	else
	{
		dmp >> m_ftime;
	}

	// stream mesh
	m_mesh.ShallowCopy(dmp, bsave);

	// stream rigid body data
	for (int i=0; i<(int) m_Obj.size(); ++i) m_Obj[i]->ShallowCopy(dmp, bsave);

	// stream contact data
	for (int i=0; i<SurfacePairInteractions(); ++i) m_CI[i]->ShallowCopy(dmp, bsave);
}

//=============================================================================
//    P A R A M E T E R   F U N C T I O N S
//=============================================================================

//-----------------------------------------------------------------------------
//! Return a pointer to the named variable

//! This function returns a pointer to a named variable. Currently, we only
//! support names of the form:
//!		material_name.parameter_name
//!		material_name.elastic.parameter_name (nested material)
//!		material_name.solid_name.parameter_name (solid mixture)
//!		material_name.solid.parameter_name (biphasic material)
//!		material_name.permeability.parameter_name (biphasic material)
//!		material_name.solid.solid_name.parameter_name (biphasic material with solid mixture)
//! The 'material_name' is a user defined name for a material.
//! The 'parameter_name' is the predefined name of the variable.
//! The keywords 'elastic', 'solid', and 'permeability' must appear as shown.
//! \todo perhaps I should use XPath to refer to material parameters ?

double* FEModel::FindParameter(const char* szparam)
{
	char szname[256];
	strcpy(szname, szparam);

	// get the material and parameter name
	char* ch = strchr((char*)szname, '.');
	if (ch == 0) return 0;
	*ch = 0;
	const char* szmat = szname;
	const char* szvar = ch+1;

	// find the material with the same name
	FEMaterial* pmat = 0;
	int nmat = -1;
	for (int i=0; i<Materials(); ++i)
	{
		pmat = GetMaterial(i);
		nmat = i;
		if (strcmp(szmat, pmat->GetName()) == 0) break;
		pmat = 0;
	}

	// make sure we found a material with the same name
	if (pmat == 0) return 0;

	// if the variable is a vector, then we require an index
	char* szarg = strchr((char*) szvar, '[');
	int index = 0;
	if (szarg)
	{
		*szarg = 0; szarg++;
		const char* ch = strchr(szarg, ']');
		assert(ch);
		index = atoi(szarg);
	}

	// find the material parameter
	ParamString sz(szvar);
	FEParam* pp = pmat->GetParameter(sz);
	if (pp) return pp->pvalue<double>(index);

	// the rigid bodies are dealt with differently
	int nrb = m_Obj.size();
	for (int i=0; i<nrb; ++i)
	{
		FEObject& ob = *m_Obj[i];
		if (ob.GetMaterialID() == nmat)
		{
			FEParam* pp = ob.GetParameter(sz);
			if (pp) return pp->pvalue<double>(index);
		}
	}

	// oh, oh, we didn't find it
	return 0;
}

//-----------------------------------------------------------------------------
void FEModel::EvaluateAllParameterLists()
{
	// evaluate material parameter lists
	for (int i=0; i<Materials(); ++i)
	{
		// get the material
		FEMaterial* pm = GetMaterial(i);

		// evaluate its parameter list
		EvaluateParameterList(pm);
	}

	// evaluate surface load parameter lists
	for (int i=0; i<SurfaceLoads(); ++i)
	{
		FEParameterList& pl = SurfaceLoad(i)->GetParameterList();
		EvaluateParameterList(pl);
	}

	// evaluate body load parameter lists
	for (int i=0; i<BodyLoads(); ++i)
	{
		FEParameterList& pl = GetBodyLoad(i)->GetParameterList();
		EvaluateParameterList(pl);
	}

	// evaluate contact interface parameter lists
	for (int i=0; i<SurfacePairInteractions(); ++i)
	{
		FEParameterList& pl = SurfacePairInteraction(i)->GetParameterList();
		EvaluateParameterList(pl);
	}

	// evaluate constraint parameter lists
	for (int i=0; i<NonlinearConstraints(); ++i)
	{
		FEParameterList& pl = NonlinearConstraint(i)->GetParameterList();
		EvaluateParameterList(pl);
	}
}

//-----------------------------------------------------------------------------
//! Evaluate a parameter list
void FEModel::EvaluateParameterList(FEParameterList &pl)
{
	list<FEParam>::iterator pi = pl.first();
	for (int j=0; j<pl.Parameters(); ++j, ++pi)
	{
		if (pi->m_nlc >= 0)
		{
			double v = GetLoadCurve(pi->m_nlc)->Value();
			switch (pi->m_itype)
			{
			case FE_PARAM_INT   : pi->value<int>() = (int) v; break;
			case FE_PARAM_DOUBLE: pi->value<double>() = pi->m_scl*v; break;
			case FE_PARAM_BOOL  : pi->value<bool>() = (v > 0? true : false); break;
			case FE_PARAM_VEC3D : pi->value<vec3d>() = pi->m_vscl*v; break;
			default: 
				assert(false);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function evaluates parameter lists. First the FECoreBase's parameter
//! list is evaluated. Then, the parameter lists of all the properties are 
//! evaluated recursively.
void FEModel::EvaluateParameterList(FECoreBase* pc)
{
	// evaluate the component's parameter list
	EvaluateParameterList(pc->GetParameterList());

	// evaluate the properties' parameter lists
	int N = pc->Properties();
	for (int i=0; i<N; ++i)
	{
		FECoreBase* pci = pc->GetProperty(i);
		if (pci) EvaluateParameterList(pci);
	}
}

//-----------------------------------------------------------------------------
// This function adds a callback routine
//
void FEModel::AddCallback(FECORE_CB_FNC pcb, unsigned int nwhen, void *pd)
{
	FECORE_CALLBACK cb;
	cb.m_pcb = pcb;
	cb.m_pd = pd;
	cb.m_nwhen = nwhen;

	m_pcb.push_back(cb);
}

//-----------------------------------------------------------------------------
// Call the callback function if there is one defined
//
void FEModel::DoCallback(unsigned int nevent)
{
	list<FECORE_CALLBACK>::iterator it = m_pcb.begin();
	for (int i=0; i<(int) m_pcb.size(); ++i, ++it)
	{
		// call the callback function
		if (it->m_nwhen & nevent) (it->m_pcb)(this, it->m_pd);
	}
}

//-----------------------------------------------------------------------------
void FEModel::SetGlobalConstant(const string& s, double v)
{
	m_Const[s] = v;
	return;
}

//-----------------------------------------------------------------------------
double FEModel::GetGlobalConstant(const string& s)
{
	return (m_Const.count(s) ? m_Const.find(s)->second : 0);
}

//-----------------------------------------------------------------------------
void FEModel::AddGlobalData(FEGlobalData* psd)
{
	m_GD.push_back(psd);
}

//-----------------------------------------------------------------------------
FEGlobalData* FEModel::GetGlobalData(int i)
{
	return m_GD[i];
}

//-----------------------------------------------------------------------------
int FEModel::GlobalDataItems()
{
	return (int) m_GD.size();
}

//-----------------------------------------------------------------------------
//! Set the title of the model
void FEModel::SetTitle(const char* sz)
{ 
	strcpy(m_sztitle, sz); 
}

//-----------------------------------------------------------------------------
//! Return the title of the model
const char* FEModel::GetTitle()
{ 
	return m_sztitle; 
}

//-----------------------------------------------------------------------------
//! Find a BC based on its ID. This is needed for restarts.
FEBoundaryCondition* FEModel::FindBC(int nid)
{
	int i;
	for (i=0; i<(int) m_DC.size(); ++i) if (m_DC[i]->GetID() == nid) return m_DC[i];

	for (i=0; i<(int) m_FC.size(); ++i) if (m_FC[i]->GetID() == nid) return m_FC[i];

	for (i=0; i<(int) m_SL.size(); ++i) if (m_SL[i]->GetID() == nid) return m_SL[i];

	for (i=0; i<(int) m_RBC.size(); ++i) if (m_RBC[i]->GetID() == nid) return m_RBC[i];

	for (i=0; i<(int) m_RDC.size(); ++i) if (m_RDC[i]->GetID() == nid) return m_RDC[i];

	for (i=0; i<(int) m_RFC.size(); ++i) if (m_RFC[i]->GetID() == nid) return m_RFC[i];

	for (i=0; i<(int) m_RN.size(); ++i) if (m_RN[i]->GetID() == nid) return m_RN[i];

	return 0;
}

//-----------------------------------------------------------------------------
//! Find a Contact Interface based on its ID. This is needed for restarts.
FESurfacePairInteraction* FEModel::FindCI(int nid)
{
	for (int i=0; i<(int) m_CI.size(); ++i) if (m_CI[i]->GetID() == nid) return m_CI[i];

	return 0;
}

//-----------------------------------------------------------------------------
//! Find a Nonlinear Constraint based on its ID. This is needed for restarts.
FENLConstraint* FEModel::FindNLC(int nid)
{
	for (int i=0; i<(int) m_NLC.size(); ++i) if (m_NLC[i]->GetID() == nid) return m_NLC[i];

	return 0;
}

//-----------------------------------------------------------------------------
//! This function copies the model data from the fem object. Note that it only copies
//! the model definition, i.e. mesh, bc's, contact interfaces, etc..
void FEModel::CopyFrom(FEModel& fem)
{
	// clear the current model data
	Clear();

	// --- Parameters ---

	// copy parameters (not sure if I need/want to copy all of these)
	m_nsolver = fem.m_nsolver;
	m_bwopt   = fem.m_bwopt;
	m_nStep   = fem.m_nStep;
	m_ftime   = fem.m_ftime;
	m_ftime0  = fem.m_ftime0;
	m_nplane_strain = fem.m_nplane_strain;
	m_debug   = fem.m_debug;
	m_ut4_alpha = fem.m_ut4_alpha;
	m_ut4_bdev  = fem.m_ut4_bdev;
	m_udghex_hg = fem.m_udghex_hg;
	m_pStep = 0;

	// --- Steps ---

	// copy the steps
	// NOTE: This does not copy the boundary conditions of the steps
	int NSTEP = fem.Steps();
	for (int i=0; i<NSTEP; ++i)
	{
		// get the type string
		FEAnalysis* ps = fem.GetStep(i);
		const char* sztype = ps->GetTypeStr();

		// create a new step
		FEAnalysis* pnew = fecore_new<FEAnalysis>(FEANALYSIS_ID, sztype, this);
		assert(pnew);

		// copy parameter data
		pnew->GetParameterList() = ps->GetParameterList();

		// copy additional info
		pnew->m_nanalysis = ps->m_nanalysis;
		pnew->m_istiffpr = ps->m_istiffpr;
		pnew->m_baugment = ps->m_baugment;

		pnew->m_ntime		= ps->m_ntime;
		pnew->m_final_time	= ps->m_final_time;
		pnew->m_dt			= ps->m_dt;
		pnew->m_dt0			= ps->m_dt0;
		pnew->m_tstart		= ps->m_tstart;
		pnew->m_tend		= ps->m_tend;
		pnew->m_bautostep	= ps->m_bautostep;
		pnew->m_iteopt		= ps->m_iteopt;
		pnew->m_dtmin		= ps->m_dtmin;
		pnew->m_dtmax		= ps->m_dtmax;
		pnew->m_ddt			= ps->m_ddt;
		pnew->m_nmplc		= ps->m_nmplc;
		pnew->m_naggr		= ps->m_naggr;

		// copy the solver
		FESolver* psolver = ps->m_psolver;
		sztype = psolver->GetTypeStr();

		// create a new solver
		FESolver* pnew_solver = fecore_new<FESolver>(FESOLVER_ID, sztype, this);
		assert(pnew_solver);
		pnew->m_psolver = pnew_solver;

		// add the step
		AddStep(pnew);
	}

	// --- Materials ---

	// copy material data
	int NMAT = fem.Materials();
	for (int i=0; i<NMAT; ++i)
	{
		// get the type info from the old material
		FEMaterial* pmat = fem.GetMaterial(i);
		const char* sztype = pmat->GetTypeStr();

		// create a new material
		FEMaterial* pnew = fecore_new<FEMaterial>(FEMATERIAL_ID, sztype, this);
		assert(pnew);

		// copy material data
		// we only copy material parameters
		pnew->GetParameterList() = pmat->GetParameterList();

		// add the material
		AddMaterial(pnew);

	}
	assert(m_MAT.size() == fem.m_MAT.size());

	// --- Mesh ---

	// copy the mesh data
	// NOTE: This will not assign materials to the new domains
	m_mesh.CopyFrom(fem.m_mesh);
	assert(m_mesh.Domains()==fem.m_mesh.Domains());

	// next, we need to assign the new materials to the new domains
	// let's first create a table of material indices for the old domains
	int NDOM = fem.m_mesh.Domains();
	vector<int> LUT(NDOM);
	for (int i=0; i<NDOM; ++i)
	{
		FEMaterial* pm = fem.m_mesh.Domain(i).GetMaterial();
		for (int j=0; j<NMAT; ++j)
		{
			if (pm == fem.m_MAT[j])
			{
				LUT[i] = j;
				break;
			}
		}
	}

	// since both the domains and the materials are created in the same order
	// we can use the lookup table to assign materials to domains
	int ND = m_mesh.Domains();
	for (int i=0; i<ND; ++i)
	{
		FEDomain& dom = m_mesh.Domain(i);
		dom.SetMaterial(m_MAT[LUT[i]]);
		assert(dom.GetMaterial());
	}

	// --- boundary conditions ---

	int NDC = fem.PrescribedBCs();
	for (int i=0; i<NDC; ++i)
	{
		FEPrescribedBC* pbc = fem.PrescribedBC(i);
		FEPrescribedBC* pnew = new FEPrescribedBC(this);

		// copy data
		pnew->bc = pbc->bc;
		pnew->br = pbc->br;
		pnew->lc = pbc->lc;
		pnew->node = pbc->node;
		pnew->r = pbc->r;
		pnew->s = pbc->s;

		// add to model
		AddPrescribedBC(pnew);
	}

	// --- Load curves ---

	// copy load curves
	int NLC = fem.LoadCurves();
	for (int i=0; i<NLC; ++i)
	{
		FELoadCurve* plc = new FELoadCurve(*fem.m_LC[i]);
		m_LC.push_back(plc);
	}
}
