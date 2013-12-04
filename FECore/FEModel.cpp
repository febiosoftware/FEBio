#include "stdafx.h"
#include "FEModel.h"
#include "FEShellDomain.h"
#include "FERigidBody.h"
#include <string>
#include "log.h"
using namespace std;

//-----------------------------------------------------------------------------
FEModel::FEModel(void)
{
	// --- Analysis Data ---
	m_pStep = 0;
	m_nStep = -1;
	m_nplane_strain = -1;	// don't use plain strain mode
	m_ftime = 0;
	m_ftime0 = 0;
	m_bwopt = 0;
}

//-----------------------------------------------------------------------------
//! Delete all dynamically allocated data
FEModel::~FEModel(void)
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
	for (i=0; i<m_RN.size  (); ++i) delete m_RN [i] ; m_RN.clear  ();
	for (i=0; i<m_NLC.size (); ++i) delete m_NLC[i] ; m_NLC.clear ();
	for (i=0; i<m_Obj.size (); ++i) delete m_Obj[i] ; m_Obj.clear ();

	// global data
	for (i=0; i<m_SD.size(); ++i) delete m_SD[i]; m_SD.clear();
	for (i=0; i<m_SBM.size(); ++i) delete m_SBM[i]; m_SBM.clear();
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

	// if the analysis is run in plain-strain mode we fix all the z-dofs of all nodes
	if (m_nplane_strain >= 0)
	{
		int bc = m_nplane_strain;
		for (int i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_BC[bc] = -1;
	}

	// find and remove isolated vertices
	int ni = m_mesh.RemoveIsolatedVertices();
	if (ni != 0) 
	{
		if (ni == 1)
			felog.printbox("WARNING", "%d isolated vertex removed.", ni);
		else
			felog.printbox("WARNING", "%d isolated vertices removed.", ni);
	}

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
	InitConstraints();

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
		for (int i=0; i<dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			if (dom.GetMaterial()->IsRigid())
			{
				el.m_nrigid = el.GetMatID();
				for (int j=0; j<el.Nodes(); ++j)
				{
					int n = el.m_node[j];
					FENode& node = m_mesh.Node(n);
					node.m_rid = el.GetMatID();
				}
			}
			else el.m_nrigid = -1;
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
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (int i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
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
	if (pmat == 0) return false;

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
			default: 
				assert(false);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function evaluates material parameter lists. Since some of the materials
//! can have other materials as sub-componenents, we need to set up a recursive
//! call to evaluate the parameter lists of the sub-materials.
void FEModel::EvaluateMaterialParameters(FEMaterial* pm)
{
	// evaluate the materials' parameter list
	EvaluateParameterList(pm->GetParameterList());

	// evaluate the material properties
	int N = pm->Properties();
	for (int i=0; i<N; ++i)
	{
		FEMaterial* pmi = pm->GetProperty(i);
		if (pmi) EvaluateMaterialParameters(pmi);
	}
}

//-----------------------------------------------------------------------------
// This function adds a callback routine
//
void FEModel::AddCallback(FEBIO_CB_FNC pcb, unsigned int nwhen, void *pd)
{
	FEBIO_CALLBACK cb;
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
	list<FEBIO_CALLBACK>::iterator it = m_pcb.begin();
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
void FEModel::AddSoluteData(FESoluteData* psd)
{
	m_SD.push_back(psd);
}

//-----------------------------------------------------------------------------
FESoluteData* FEModel::FindSoluteData(int nid)
{
	int i;
	for (i=0; i<(int) m_SD.size(); ++i) if (m_SD[i]->m_nID == nid) return m_SD[i];
	
	return 0;
}

//-----------------------------------------------------------------------------
void FEModel::AddSBMData(FESBMData* psd)
{
	m_SBM.push_back(psd);
}

//-----------------------------------------------------------------------------
FESBMData* FEModel::FindSBMData(int nid)
{
	int i;
	for (i=0; i<(int) m_SBM.size(); ++i) if (m_SBM[i]->m_nID == nid) return m_SBM[i];
	
	return 0;
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

	for (i=0; i<(int) m_RDC.size(); ++i) if (m_RDC[i]->GetID() == nid) return m_RDC[i];

	for (i=0; i<(int) m_RFC.size(); ++i) if (m_RFC[i]->GetID() == nid) return m_RFC[i];

	for (i=0; i<(int) m_RN.size(); ++i) if (m_RN[i]->GetID() == nid) return m_RN[i];

	return 0;
}
