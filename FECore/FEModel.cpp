#include "stdafx.h"
#include "FEModel.h"
#include "FEShellDomain.h"
#include "FERigidBody.h"
#include <string>
#include "log.h"
#include "FECoreKernel.h"
using namespace std;

// set the default linear solver (0 is equivalent to skyline solver)
int FEModel::m_ndefault_solver = 0;

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
	m_nsolver = m_ndefault_solver;

	// create a rigid system
	// TODO: Perhaps I can set it inially to zero and only allocate it if a rigid body is defined
	m_prs = new FERigidSystem(this);
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
	for (i=0; i<m_BC.size  (); ++i) delete m_BC [i] ; m_BC.clear  ();
	for (i=0; i<m_DC.size  (); ++i) delete m_DC [i] ; m_DC.clear  ();
	for (i=0; i<m_IC.size  (); ++i) delete m_IC [i] ; m_IC.clear  ();
	for (i=0; i<m_FC.size  (); ++i) delete m_FC [i] ; m_FC.clear  ();
	for (i=0; i<m_SL.size  (); ++i) delete m_SL [i] ; m_SL.clear  ();
	for (i=0; i<m_ML.size  (); ++i) delete m_ML [i] ; m_ML.clear  ();
	for (i=0; i<m_RDC.size (); ++i) delete m_RDC[i] ; m_RDC.clear ();
	for (i=0; i<m_RBV.size (); ++i) delete m_RBV[i] ; m_RBV.clear ();
	for (i=0; i<m_RBW.size (); ++i) delete m_RBW[i] ; m_RBW.clear ();
	for (i=0; i<m_RN.size  (); ++i) delete m_RN [i] ; m_RN.clear  ();
	for (i=0; i<m_NLC.size (); ++i) delete m_NLC[i] ; m_NLC.clear ();

	// global data
	for (i=0; i<m_GD.size(); ++i) delete m_GD[i]; m_GD.clear();
	m_Const.clear();

	// clear the rigid system (if there is one)
	if (m_prs) m_prs->Clear();
}

//-----------------------------------------------------------------------------
void FEModel::AddFixedBC(int node, int bc)
{
	FEFixedBC* pbc = dynamic_cast<FEFixedBC*>(fecore_new<FEBoundaryCondition>(FEBC_ID, "fix", this));
	pbc->SetDOF(bc);
	pbc->AddNode(node);
	AddFixedBC(pbc);
}

//-----------------------------------------------------------------------------
void FEModel::ClearBCs()
{
	for (size_t i=0; i<m_DC.size  (); ++i) delete m_DC[i];
	for (size_t i=0; i<m_BC.size  (); ++i) delete m_BC[i];
	m_DC.clear();
	m_BC.clear();
}

//-----------------------------------------------------------------------------
bool FEModel::Init()
{
	// intitialize time
	m_ftime = 0;
	m_ftime0 = 0;

	// check step data
	for (int i=0; i<(int) m_Step.size(); ++i)
	{
		FEAnalysis& step = *m_Step[i];
		if (step.Init() == false) return false;
	}

	// evaluate all loadcurves at the initial time
	for (int i=0; i<LoadCurves(); ++i) m_LC[i]->Evaluate(0);
    
    // evaluate all parameter lists
    if (EvaluateAllParameterLists() == false) return false;

	// create and initialize the rigid body data
	// NOTE: Do this first, since some BC's look at the nodes' rigid id.
	if (m_prs && (m_prs->Init() == false)) return false;

	// validate BC's
	if (InitBCs() == false) return false;

	// initialize material data
	// NOTE: call this before InitMesh since we need to initialize the FECoordSysMap
	//       before we can calculate the local element coordinate systems.
	// NOTE: This must be called after the rigid system is initialiazed since the rigid materials will
	//       reference the rigid bodies
	if (InitMaterials() == false) return false;

	// initialize model loads
	// NOTE: This must be called after the InitMaterials since the COM of the rigid bodies
	//       are set in that function. 
	if (InitModelLoads() == false) return false;

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
	if (InitConstraints() == false) return false;

	// activate all permanent dofs
	Activate();

	// do the callback
	DoCallback(CB_INIT);

	return true;
}

//-----------------------------------------------------------------------------
//! See if the BC's are setup correctly.
bool FEModel::InitBCs()
{
	// if the analysis is run in plane-strain mode we fix all the z-dofs of all nodes
//--> I want to remove this. Perhaps I can make a new class for this. (derived from FEFixedBC?)
	if (m_nplane_strain >= 0)
	{
		int bc = m_nplane_strain;
		for (int i=0; i<m_mesh.Nodes(); ++i) AddFixedBC(i, bc);
	}
//-->

	// get the number of loadcurves
	int NLC = LoadCurves();

	// check the prescribed BC's
	int NBC = PrescribedBCs();
	for (int i=0; i<NBC; ++i)
	{
		FEPrescribedBC* pbc = PrescribedBC(i);
		if (pbc->Init() == false) return false;
	}

	// check the nodal loads
	int NNL = NodalLoads();
	for (int i=0; i<NNL; ++i)
	{
		FENodalLoad* pbc = NodalLoad(i);
		if (pbc->Init() == false) return false;
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
	for (int i=0; i<m_mesh.Domains(); ++i)
	{
		if (m_mesh.Domain(i).Initialize(*this) == false) return false;
	}

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
//! Initialize rigid force data
bool FEModel::InitModelLoads()
{
	// call the Init() function of all rigid forces
	for (int i=0; i<(int) m_ML.size(); ++i)
	{
		FEModelLoad& FC = *m_ML[i];
		if (FC.Init() == false) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Initializes contact data
bool FEModel::InitContact()
{
	// loop over all contact interfaces
	for (int i=0; i<SurfacePairInteractions(); ++i)
	{
		// get the contact interface
		FESurfacePairInteraction& ci = *m_CI[i];

		// initializes contact interface data
		if (ci.Init() == false) return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize the nonlinear constraints.
//! This function is called during model initialization (\sa FEModel::Init)
bool FEModel::InitConstraints()
{
	for (int i=0; i<(int) m_NLC.size(); ++i)
	{
		FENLConstraint* plc = m_NLC[i];

		// initialize
		if (plc->Init() == false) return false;
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
		if (m_pStep->Activate() == false)
		{
			bconv = false;
			break;
		}

		// solve the analaysis step
		bconv = m_pStep->Solve();

		// break if the step has failed
		if (bconv == false) break;

		// wrap it up
		m_pStep->Deactivate();
	}

	return bconv;
}

//-----------------------------------------------------------------------------
// Model activation.
// BC's that are not assigned to a step will not have their Activate() member called
// so we do it here. This function is called in Init() and Reset()
void FEModel::Activate()
{
	// fixed dofs
	for (int i=0; i<(int) m_BC.size(); ++i)
	{
		FEFixedBC& bc = *m_BC[i];
		if (bc.IsActive()) bc.Activate();
	}

	// initial conditions
	// Must be activated before prescribed BC's
	// since relative prescribed BC's use the initial values
	for (int i=0; i<(int) m_IC.size(); ++i)
	{
		FEInitialCondition& ic = *m_IC[i];
		if (ic.IsActive()) ic.Activate();
	}

	// prescribed dofs
	for (int i=0; i<(int) m_DC.size(); ++i)
	{
		FEPrescribedBC& dc = *m_DC[i];
		if (dc.IsActive()) dc.Activate();
	}

	// model loads
	for (int i=0; i<(int) m_ML.size(); ++i)
	{
		FEModelLoad& FC = *m_ML[i];
		if (FC.IsActive()) FC.Activate();
	}

	// nonlinear constraints
	for (int i=0; i<(int) m_NLC.size(); ++i)
	{
		FENLConstraint* plc = m_NLC[i];
		if (plc->IsActive()) plc->Activate();
	}

	// contact interfaces
	for (int i=0; i<SurfacePairInteractions(); ++i)
	{
		FESurfacePairInteraction& ci = *m_CI[i];
		if (ci.IsActive()) ci.Activate();
	}

	// rigid nodes
	for (int i=0; i<(int) m_RN.size(); ++i)
	{
		FERigidNode& rn = *m_RN[i];
		if (rn.IsActive()) rn.Activate();
	}

	// rigid body displacements
	for (int i=0; i<(int)m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& rc = *m_RDC[i];
		if (rc.IsActive()) rc.Activate();
	}

	// fixed rigid body dofs
	for (int i=0; i<(int)m_RBC.size(); ++i)
	{
		FERigidBodyFixedBC& rc = *m_RBC[i];
		if (rc.IsActive()) rc.Activate();
	}

	// initial rigid velocity
	for (int i=0; i<(int) m_RBV.size(); ++i)
	{
		FERigidBodyVelocity& RV = *m_RBV[i];
		if (RV.IsActive()) RV.Activate();
	}

	// initial rigid angular velocity
	for (int i=0; i<(int) m_RBW.size(); ++i)
	{
		FERigidBodyAngularVelocity& RW = *m_RBW[i];
		if (RW.IsActive()) RW.Activate();
	}
}

//-----------------------------------------------------------------------------
//! \todo Do I really need this function. I think calling FEModel::Init achieves the
//! same effect.
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
	if (m_prs && (m_prs->Reset() == false)) return false;

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

	// Call Activate() to activate all permanent BC's
	Activate();

	return true;
}

//-----------------------------------------------------------------------------
//! Get the current time information.
FETimePoint FEModel::GetTime()
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
	if (m_prs) m_prs->ShallowCopy(dmp, bsave);

	// stream contact data
	for (int i=0; i<SurfacePairInteractions(); ++i) m_CI[i]->ShallowCopy(dmp, bsave);

	// stream nonlinear constraints
	for (int i=0; i<NonlinearConstraints(); ++i) m_NLC[i]->ShallowCopy(dmp, bsave);
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

	// if we get here we didn't find
	// Let's try the rigid system
	if (m_prs) return m_prs->FindParameter(nmat, sz, index);

	// oh, oh, we didn't find it
	return 0;
}

//-----------------------------------------------------------------------------
//! Evaluates all load curves at the specified time
void FEModel::EvaluateLoadCurves(double time)
{
	const int NLC = LoadCurves();
	for (int i=0; i<NLC; ++i) GetLoadCurve(i)->Evaluate(time);
}

//-----------------------------------------------------------------------------
bool FEModel::EvaluateAllParameterLists()
{
	// evaluate material parameter lists
	for (int i=0; i<Materials(); ++i)
	{
		// get the material
		FEMaterial* pm = GetMaterial(i);

		// evaluate its parameter list
		if (EvaluateParameterList(pm) == false) return false;
	}

	// evaluate surface load parameter lists
	for (int i=0; i<SurfaceLoads(); ++i)
	{
		FEParameterList& pl = SurfaceLoad(i)->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
	}

	// evaluate body load parameter lists
	for (int i=0; i<BodyLoads(); ++i)
	{
		FEParameterList& pl = GetBodyLoad(i)->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
	}

	// evaluate contact interface parameter lists
	for (int i=0; i<SurfacePairInteractions(); ++i)
	{
		FEParameterList& pl = SurfacePairInteraction(i)->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
	}

	// evaluate constraint parameter lists
	for (int i=0; i<NonlinearConstraints(); ++i)
	{
		FEParameterList& pl = NonlinearConstraint(i)->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
	}

	// evaluate model loads
	for (int i=0; i<(int)m_ML.size(); ++i)
	{
		FEParameterList& pl = m_ML[i]->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Evaluate a parameter list
bool FEModel::EvaluateParameterList(FEParameterList &pl)
{
	const int NLC = LoadCurves();
	list<FEParam>::iterator pi = pl.first();
	for (int j=0; j<pl.Parameters(); ++j, ++pi)
	{
		int nlc = pi->m_nlc;
		if (pi->m_nlc >= 0)
		{
			if (nlc >= NLC) return false;

			double v = GetLoadCurve(nlc)->Value();
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

	return true;
}

//-----------------------------------------------------------------------------
//! This function evaluates parameter lists. First the FECoreBase's parameter
//! list is evaluated. Then, the parameter lists of all the properties are 
//! evaluated recursively.
bool FEModel::EvaluateParameterList(FECoreBase* pc)
{
	// evaluate the component's parameter list
	if (EvaluateParameterList(pc->GetParameterList()) == false) return false;

	// evaluate the properties' parameter lists
	int N = pc->Properties();
	for (int i=0; i<N; ++i)
	{
		FECoreBase* pci = pc->GetProperty(i);
		if (pci)
		{
			if (EvaluateParameterList(pci) == false) return false;
		}
	}

	return true;
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
		if (it->m_nwhen & nevent) (it->m_pcb)(this, nevent, it->m_pd);
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
FEModelComponent* FEModel::FindModelComponent(int nid)
{
	int i;
	for (i=0; i<(int) m_BC.size (); ++i) if (m_BC [i]->GetClassID() == nid) return m_BC [i];
	for (i=0; i<(int) m_DC.size (); ++i) if (m_DC [i]->GetClassID() == nid) return m_DC [i];
	for (i=0; i<(int) m_IC.size (); ++i) if (m_IC [i]->GetClassID() == nid) return m_IC [i];
	for (i=0; i<(int) m_FC.size (); ++i) if (m_FC [i]->GetClassID() == nid) return m_FC [i];
	for (i=0; i<(int) m_SL.size (); ++i) if (m_SL [i]->GetClassID() == nid) return m_SL [i];
	for (i=0; i<(int) m_ML.size (); ++i) if (m_ML [i]->GetClassID() == nid) return m_ML [i];
	for (i=0; i<(int) m_RBC.size(); ++i) if (m_RBC[i]->GetClassID() == nid) return m_RBC[i];
	for (i=0; i<(int) m_RDC.size(); ++i) if (m_RDC[i]->GetClassID() == nid) return m_RDC[i];
	for (i=0; i<(int) m_RN.size (); ++i) if (m_RN [i]->GetClassID() == nid) return m_RN [i];
	for (i=0; i<(int) m_CI.size (); ++i) if (m_CI [i]->GetClassID() == nid) return m_CI [i];
	for (i=0; i<(int) m_NLC.size(); ++i) if (m_NLC[i]->GetClassID() == nid) return m_NLC[i];
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

		pnew->SetID(pmat->GetID());

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
		FEPrescribedBC* pnew = new FEPrescribedBC(this, *pbc);

		// add to model
		AddPrescribedBC(pnew);
	}

	// --- contact interfaces ---
	int NCI = fem.SurfacePairInteractions();
	for (int i=0; i<NCI; ++i)
	{
		// get the next interaction
		FESurfacePairInteraction* pci = fem.SurfacePairInteraction(i);
		const char* sztype = pci->GetTypeStr();

		// create a new contact interface
		FESurfacePairInteraction* pnew = fecore_new<FESurfacePairInteraction>(FESURFACEPAIRINTERACTION_ID, sztype, this);
		assert(pnew);

		// create a copy
		pnew->CopyFrom(pci);

		// add the new interface
		AddSurfacePairInteraction(pnew);

		// add the surfaces to the surface list
		m_mesh.AddSurface(pnew->GetMasterSurface());
		m_mesh.AddSurface(pnew->GetSlaveSurface ());
	}

	// --- nonlinear constraints ---
	int NLC = fem.NonlinearConstraints();
	for (int i=0; i<NLC; ++i)
	{
		// get the next constraint
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		const char* sztype = plc->GetTypeStr();

		// create a new nonlinear constraint
		FENLConstraint* plc_new = fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, sztype, this);
		assert(plc_new);

		// create a copy
		plc_new->CopyFrom(plc);

		// add the nonlinear constraint
		AddNonlinearConstraint(plc_new);

		// add the surface to the mesh (if any)
		FESurface* ps = plc_new->GetSurface(0);
		if (ps) m_mesh.AddSurface(ps);
	}

	// --- Load curves ---

	// copy load curves
	int NLD = fem.LoadCurves();
	for (int i=0; i<NLD; ++i)
	{
		FELoadCurve* plc = new FELoadCurve(*fem.m_LC[i]);
		m_LC.push_back(plc);
	}
}
