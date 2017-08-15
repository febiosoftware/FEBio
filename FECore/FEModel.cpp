#include "stdafx.h"
#include "FEModel.h"
#include "FEDataLoadCurve.h"
#include "FEMaterial.h"
#include "FERigidSystem.h"
#include "FEModelLoad.h"
#include "BC.h"
#include "FESurfaceLoad.h"
#include "FEEdgeLoad.h"
#include "FEBodyLoad.h"
#include "FEInitialCondition.h"
#include "FESurfacePairInteraction.h"
#include "FENLConstraint.h"
#include "FEAnalysis.h"
#include "FEGlobalData.h"
#include "FECoreKernel.h"
#include "FELinearConstraintManager.h"
#include "log.h"
#include "FERigidBody.h"
#include "FEModelData.h"
#include "FEDataArray.h"
#include <string>
#include <map>
using namespace std;

//-----------------------------------------------------------------------------
// Implementation class for the FEModel class
class FEModel::Implementation
{
public:
	Implementation(FEModel* fem) : m_fem(fem)
	{
		m_sztitle[0] = 0;

		// --- Analysis Data ---
		m_pStep = 0;
		m_nStep = -1;
		m_ftime = 0;
		m_ftime0 = 0;
		m_bwopt = 0;

		// additional data
		m_linearSolver = FECoreKernel::m_ndefault_solver;

		// create a rigid system
		m_prs = new FERigidSystem(fem);

		// create the linear constraint manager
		m_LCM = new FELinearConstraintManager(fem);
	}

public:
	// helper functions for serialization
	void SerializeLoadData    (DumpStream& ar);
	void SerializeGlobals     (DumpStream& ar);
	void SerializeMaterials   (DumpStream& ar);
	void SerializeGeometry    (DumpStream& ar);
	void SerializeContactData (DumpStream& ar);
	void SerializeBoundaryData(DumpStream& ar);
	void SerializeAnalysisData(DumpStream& ar);


public: // TODO: Find a better place for these parameters
	// I want to make this parameter part of the FEAnalysis, since 
	// it could be different analysis steps (in multi-step problems) may
	// require different solvers.
	int		m_linearSolver;			//!< type of (linear) solver selected

	int		m_bwopt;			//!< bandwidth optimization flag
	double	m_ftime;			//!< current time value
	double	m_ftime0;			//!< start time of current step

public:
	FEVecPropertyT<FEMaterial>					m_MAT;		//!< array of materials
	FEVecPropertyT<FEFixedBC>					m_BC;		//!< fixed constraints
	FEVecPropertyT<FEPrescribedBC>				m_DC;		//!< prescribed constraints
	FEVecPropertyT<FENodalLoad>					m_FC;		//!< concentrated nodal loads
	FEVecPropertyT<FESurfaceLoad>				m_SL;		//!< surface loads
	FEVecPropertyT<FEEdgeLoad>					m_EL;		//!< edge loads
	FEVecPropertyT<FEBodyLoad>					m_BL;		//!< body load data
	FEVecPropertyT<FEInitialCondition>			m_IC;		//!< initial conditions
	FEVecPropertyT<FESurfacePairInteraction>	m_CI;		//!< contact interface array
	FEVecPropertyT<FENLConstraint>				m_NLC;		//!< nonlinear constraints
	FEVecPropertyT<FEModelLoad>					m_ML;		//!< model loads
	FEVecPropertyT<FELoadCurve>					m_LC;		//!< load curve data
	FEVecPropertyT<FEAnalysis>					m_Step;		//!< array of analysis steps
	FEVecPropertyT<FEModelData>					m_Data;		//!< the model output data

public:
	FEAnalysis*		m_pStep;	//!< pointer to current analysis step
	int				m_nStep;	//!< current index of analysis step

public:
	// The model
	FEModel*	m_fem;

	// DOFS data
	DOFS	m_dofs;				//!< list of degree of freedoms in this model

	// Geometry data
	FEMesh		m_mesh;			//!< the one and only FE mesh

	// the rigid body system
	FERigidSystem*		m_prs;	//!< the rigid body system manages rigid bodies

	// linear constraint data
	FELinearConstraintManager*	m_LCM;

public:
	char	m_sztitle[MAX_STRING];	//!< problem title

public: // Global Data
	std::map<string, double> m_Const;	//!< Global model constants
	vector<FEGlobalData*>	m_GD;		//!< global data structures

public:
	vector<pair<string, FEDataArray*> >	m_DataArray;
};

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEModel, FECoreBase)
	ADD_PARAMETER(m_imp->m_ftime, FE_PARAM_DOUBLE, "time");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEModel::FEModel(void) : FECoreBase(FEMODEL_ID), m_imp(new FEModel::Implementation(this))
{
	m_ut4_alpha = 0.05;
	m_ut4_bdev = false;
	m_udghex_hg = 1.0;

	// set the name
	SetName("fem");

	// Add all properties
	AddProperty(&m_imp->m_MAT, "material");
	AddProperty(&m_imp->m_BC , "bc_fixed");
	AddProperty(&m_imp->m_DC , "bc_prescribed");
	AddProperty(&m_imp->m_FC , "nodal_load");
	AddProperty(&m_imp->m_SL , "surface_load");
	AddProperty(&m_imp->m_EL , "edge_load");
	AddProperty(&m_imp->m_BL , "body_load");
	AddProperty(&m_imp->m_IC , "initial" );
	AddProperty(&m_imp->m_CI , "contact");
	AddProperty(&m_imp->m_NLC, "constraint");
	AddProperty(&m_imp->m_ML , "model_load");
	AddProperty(&m_imp->m_LC , "loadcurve");
	AddProperty(&m_imp->m_Step, "step");
	AddProperty(&m_imp->m_Data, "data");
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
	// clear dofs
	m_imp->m_dofs.Reset();

	// clear all properties
	m_imp->m_MAT.Clear();
	m_imp->m_BC.Clear();
	m_imp->m_DC.Clear();
	m_imp->m_FC.Clear();
	m_imp->m_SL.Clear();
	m_imp->m_EL.Clear();
	m_imp->m_BL.Clear();
	m_imp->m_IC.Clear();
	m_imp->m_CI.Clear();
	m_imp->m_NLC.Clear();
	m_imp->m_ML.Clear();
	m_imp->m_LC.Clear();
	m_imp->m_Step.Clear();

	// global data
	for (size_t i = 0; i<m_imp->m_GD.size(); ++i) delete m_imp->m_GD[i]; m_imp->m_GD.clear();
	m_imp->m_Const.clear();

	// clear the rigid system (if there is one)
	if (m_imp->m_prs) m_imp->m_prs->Clear();

	// clear the linear constraints
	if (m_imp->m_LCM) m_imp->m_LCM->Clear();

	// clear the mesh
	m_imp->m_mesh.Clear();
}

//-----------------------------------------------------------------------------
int FEModel::GetLinearSolverType() const
{
	return m_imp->m_linearSolver;
}

//-----------------------------------------------------------------------------
//! set the linear solver type
void FEModel::SetLinearSolverType(int ntype)
{
	m_imp->m_linearSolver = ntype;
}

//-----------------------------------------------------------------------------
//! see if we need to optimize bandwidth of linear system
bool FEModel::OptimizeBandwidth() const
{
	return (m_imp->m_bwopt == 1);
}

//-----------------------------------------------------------------------------
void FEModel::SetOptimizeBandwidth(bool b)
{
	m_imp->m_bwopt = b;
}

//-----------------------------------------------------------------------------
int FEModel::FixedBCs() { return (int)m_imp->m_BC.size(); }

//-----------------------------------------------------------------------------
FEFixedBC* FEModel::FixedBC(int i) { return m_imp->m_BC[i]; }

//-----------------------------------------------------------------------------
void FEModel::AddFixedBC(FEFixedBC* pbc) { m_imp->m_BC.AddProperty(pbc); }

//-----------------------------------------------------------------------------
int FEModel::PrescribedBCs() { return (int)m_imp->m_DC.size(); }

//-----------------------------------------------------------------------------
FEPrescribedBC* FEModel::PrescribedBC(int i) { return m_imp->m_DC[i]; }

//-----------------------------------------------------------------------------
void FEModel::AddPrescribedBC(FEPrescribedBC* pbc) { m_imp->m_DC.AddProperty(pbc); }

//-----------------------------------------------------------------------------
int FEModel::InitialConditions() { return (int)m_imp->m_IC.size(); }

//-----------------------------------------------------------------------------
FEInitialCondition* FEModel::InitialCondition(int i) { return m_imp->m_IC[i]; }

//-----------------------------------------------------------------------------
void FEModel::AddInitialCondition(FEInitialCondition* pbc) { m_imp->m_IC.AddProperty(pbc); }

//-----------------------------------------------------------------------------
int FEModel::NodalLoads() { return (int)m_imp->m_FC.size(); }

//-----------------------------------------------------------------------------
FENodalLoad* FEModel::NodalLoad(int i) { return m_imp->m_FC[i]; }

//-----------------------------------------------------------------------------
void FEModel::AddNodalLoad(FENodalLoad* pfc) { m_imp->m_FC.AddProperty(pfc); }

//-----------------------------------------------------------------------------
int FEModel::SurfaceLoads() { return (int)m_imp->m_SL.size(); }

//-----------------------------------------------------------------------------
FESurfaceLoad* FEModel::SurfaceLoad(int i) { return m_imp->m_SL[i]; }

//-----------------------------------------------------------------------------
void FEModel::AddSurfaceLoad(FESurfaceLoad* psl) { m_imp->m_SL.AddProperty(psl); }

//-----------------------------------------------------------------------------
int FEModel::EdgeLoads() { return (int)m_imp->m_EL.size(); }

//-----------------------------------------------------------------------------
FEEdgeLoad* FEModel::EdgeLoad(int i) { return m_imp->m_EL[i]; }

//-----------------------------------------------------------------------------
void FEModel::AddEdgeLoad(FEEdgeLoad* psl) { m_imp->m_EL.AddProperty(psl); }

//-----------------------------------------------------------------------------
//! Add a body load to the model
void FEModel::AddBodyLoad(FEBodyLoad* pf) { m_imp->m_BL.AddProperty(pf); }

//-----------------------------------------------------------------------------
//! get the number of body loads
int FEModel::BodyLoads() { return (int)m_imp->m_BL.size(); }

//-----------------------------------------------------------------------------
//! return a pointer to a body load
FEBodyLoad* FEModel::GetBodyLoad(int i) { return m_imp->m_BL[i]; }

//-----------------------------------------------------------------------------
//! retrieve the number of steps
int FEModel::Steps() { return (int)m_imp->m_Step.size(); }

//-----------------------------------------------------------------------------
//! clear the steps
void FEModel::ClearSteps() { m_imp->m_Step.Clear(); }

//-----------------------------------------------------------------------------
//! Add an analysis step
void FEModel::AddStep(FEAnalysis* pstep) { m_imp->m_Step.AddProperty(pstep); }

//-----------------------------------------------------------------------------
//! Get a particular step
FEAnalysis* FEModel::GetStep(int i) { return m_imp->m_Step[i]; }

//-----------------------------------------------------------------------------
//! Get the current step
FEAnalysis* FEModel::GetCurrentStep() { return m_imp->m_pStep; }

//-----------------------------------------------------------------------------
//! Set the current step
void FEModel::SetCurrentStep(FEAnalysis* pstep) { m_imp->m_pStep = pstep; }

//-----------------------------------------------------------------------------
//! Set the current step index
int FEModel::GetCurrentStepIndex() const
{
	return m_imp->m_nStep;
}

//-----------------------------------------------------------------------------
//! Set the current step index
void FEModel::SetCurrentStepIndex(int n)
{
	m_imp->m_nStep = n;
}

//-----------------------------------------------------------------------------
//! return number of surface pair interactions
int FEModel::SurfacePairInteractions() { return (int)m_imp->m_CI.size(); }

//-----------------------------------------------------------------------------
//! retrive a surface pair interaction
FESurfacePairInteraction* FEModel::SurfacePairInteraction(int i) { return m_imp->m_CI[i]; }

//-----------------------------------------------------------------------------
//! Add a surface pair interaction
void FEModel::AddSurfacePairInteraction(FESurfacePairInteraction* pci) { m_imp->m_CI.AddProperty(pci); }

//-----------------------------------------------------------------------------
//! return number of nonlinear constraints
int FEModel::NonlinearConstraints() { return (int)m_imp->m_NLC.size(); }

//-----------------------------------------------------------------------------
//! retrieve a nonlinear constraint
FENLConstraint* FEModel::NonlinearConstraint(int i) { return m_imp->m_NLC[i]; }

//-----------------------------------------------------------------------------
//! add a nonlinear constraint
void FEModel::AddNonlinearConstraint(FENLConstraint* pnlc) { m_imp->m_NLC.AddProperty(pnlc); }

//-----------------------------------------------------------------------------
//! return the number of model loads
int FEModel::ModelLoads() { return (int)m_imp->m_ML.size(); }

//-----------------------------------------------------------------------------
//! retrieve a model load
FEModelLoad* FEModel::ModelLoad(int i) { return m_imp->m_ML[i]; }

//-----------------------------------------------------------------------------
//! Add a model load
void FEModel::AddModelLoad(FEModelLoad* pml) { m_imp->m_ML.AddProperty(pml); }

//-----------------------------------------------------------------------------
// get the FE mesh
FEMesh& FEModel::GetMesh() { return m_imp->m_mesh; }

//-----------------------------------------------------------------------------
// get the rigid system
FERigidSystem* FEModel::GetRigidSystem() { return m_imp->m_prs; }

//-----------------------------------------------------------------------------
FELinearConstraintManager& FEModel::GetLinearConstraintManager() { return *m_imp->m_LCM; }

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
	m_imp->m_DC.Clear();
	m_imp->m_BC.Clear();
}

//-----------------------------------------------------------------------------
bool FEModel::Init()
{
	// intitialize time
	m_imp->m_ftime = 0;
	m_imp->m_ftime0 = 0;

	// initialize global data
	// TODO: I'd like to do this here for consistency, but
	//       the problem is that solute dofs (i.e. concentration dofs) have
	//       to be allocated before the materials are read in.
	//       So right now the Init function is called when the solute data is created.
/*	for (int i=0; i<(int) m_GD.size(); ++i)
	{
		FEGlobalData* pd = m_GD[i]; assert(pd);
		if (pd->Init() == false) return false;
	}
*/
	// check step data
	for (int i = 0; i<(int)m_imp->m_Step.size(); ++i)
	{
		FEAnalysis& step = *m_imp->m_Step[i];
		if (step.Init() == false) return false;
	}

	// evaluate all loadcurves at the initial time
	for (int i = 0; i<LoadCurves(); ++i) m_imp->m_LC[i]->Evaluate(0);
    
    // evaluate all parameter lists
    if (EvaluateAllParameterLists() == false) return false;

	// create and initialize the rigid body data
	// NOTE: Do this first, since some BC's look at the nodes' rigid id.
	if (m_imp->m_prs && (m_imp->m_prs->Init() == false)) return false;

	// validate BC's
	if (InitBCs() == false) return false;

	// initialize material data
	// NOTE: call this before FEMesh::Init() since we need to initialize the FECoordSysMap
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
	if (m_imp->m_mesh.Init() == false) return false;

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

    // check the surface loads
    int NSL = SurfaceLoads();
    for (int i=0; i<NSL; ++i)
    {
        FESurfaceLoad* pbc = SurfaceLoad(i);
        if (pbc->Init() == false) return false;
    }
    
    // check the edge loads
    int NEL = EdgeLoads();
    for (int i=0; i<NEL; ++i)
    {
        FEEdgeLoad* pbc = EdgeLoad(i);
        if (pbc->Init() == false) return false;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEModel::AddMaterial(FEMaterial* pm) 
{ 
	m_imp->m_MAT.AddProperty(pm);
}

//-----------------------------------------------------------------------------
//! get the number of materials
int FEModel::Materials() { return (int)m_imp->m_MAT.size(); }

//-----------------------------------------------------------------------------
//! return a pointer to a material
FEMaterial* FEModel::GetMaterial(int i) { return m_imp->m_MAT[i]; }

//-----------------------------------------------------------------------------
FEMaterial* FEModel::FindMaterial(int nid)
{
	for (int i = 0; i<Materials(); ++i)
	{
		FEMaterial* pm = GetMaterial(i);
		if (pm->GetID() == nid) return pm;
	}
	return 0;
}

//-----------------------------------------------------------------------------
FEMaterial* FEModel::FindMaterial(const char* sz)
{
	if (sz == 0) return 0;
	for (int i = 0; i<Materials(); ++i)
	{
		FEMaterial* pm = GetMaterial(i);
		const char* szname = pm->GetName();
		if (szname && (strcmp(szname, sz) == 0)) return pm;
	}
	return 0;
}

//-----------------------------------------------------------------------------
FESurfaceLoad* FEModel::FindSurfaceLoad(const char* sz)
{
	if (sz == 0) return 0;
	for (int i = 0; i<SurfaceLoads(); ++i)
	{
		FESurfaceLoad* pl = SurfaceLoad(i);
		const char* szname = pl->GetName();
		if (szname && (strcmp(szname, sz) == 0)) return pl;
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! Initialize material data (This also does an initial validation).
bool FEModel::InitMaterials()
{
	// initialize material data
	for (int i=0; i<Materials(); ++i)
	{
		// get the material
		FEMaterial* pmat = GetMaterial(i);

		// initialize material data
		if (pmat->Init() == false)
		{
			const char* szerr = fecore_get_error_string();
			if (szerr == 0) szerr = "unknown error";
			felog.printf("Failed initializing material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			felog.printf("ERROR: %s\n\n", szerr);
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! validate material data
bool FEModel::ValidateMaterials()
{
	// initialize material data
	for (int i=0; i<Materials(); ++i)
	{
		// get the material
		FEMaterial* pmat = GetMaterial(i);

		// initialize material data
		if (pmat->Validate() == false)
		{
			const char* szerr = fecore_get_error_string();
			if (szerr == 0) szerr = "unknown error";
			felog.printf("Failed validating material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			felog.printf("ERROR: %s\n\n", szerr);
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Add a loadcurve to the model
void FEModel::AddLoadCurve(FELoadCurve* plc) 
{ 
	m_imp->m_LC.AddProperty(plc) ;
}

//-----------------------------------------------------------------------------
//! get a loadcurve
FELoadCurve* FEModel::GetLoadCurve(int i)
{ 
	return m_imp->m_LC[i];
}

//-----------------------------------------------------------------------------
//! get the number of loadcurves
int FEModel::LoadCurves() const 
{ 
	return (int)m_imp->m_LC.size();
}

//-----------------------------------------------------------------------------
//! Initialize rigid force data
bool FEModel::InitModelLoads()
{
	// call the Init() function of all rigid forces
	for (int i=0; i<ModelLoads(); ++i)
	{
		FEModelLoad& FC = *ModelLoad(i);
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
		FESurfacePairInteraction& ci = *SurfacePairInteraction(i);

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
	for (int i=0; i<NonlinearConstraints(); ++i)
	{
		FENLConstraint* plc = NonlinearConstraint(i);

		// initialize
		if (plc->Init() == false) return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEModel::InitBodyLoads()
{
	for (int i=0; i<BodyLoads(); ++i)
	{
		if (GetBodyLoad(i)->Init() == false) return false;
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
	for (size_t nstep = m_imp->m_nStep; nstep < Steps(); ++nstep)
	{
		// set the current analysis step
		m_imp->m_nStep = (int)nstep;
		m_imp->m_pStep = m_imp->m_Step[(int)nstep];

		// intitialize step data
		if (m_imp->m_pStep->Activate() == false)
		{
			bconv = false;
			break;
		}

		DoCallback(CB_STEP_ACTIVE);

		// solve the analaysis step
		bconv = m_imp->m_pStep->Solve();

		// break if the step has failed
		if (bconv == false) break;

		// do callbacks
		DoCallback(CB_STEP_SOLVED);

		// wrap it up
		m_imp->m_pStep->Deactivate();
	}

	// do the callbacks
	DoCallback(CB_SOLVED);

	return bconv;
}

//-----------------------------------------------------------------------------
// Model activation.
// BC's that are not assigned to a step will not have their Activate() member called
// so we do it here. This function is called in Init() and Reset()
void FEModel::Activate()
{
	// fixed dofs
	for (int i=0; i<FixedBCs(); ++i)
	{
		FEFixedBC& bc = *FixedBC(i);
		if (bc.IsActive()) bc.Activate();
	}

	// initial conditions
	// Must be activated before prescribed BC's
	// since relative prescribed BC's use the initial values
	for (int i=0; i<InitialConditions(); ++i)
	{
		FEInitialCondition& ic = *InitialCondition(i);
		if (ic.IsActive()) ic.Activate();
	}

	// prescribed dofs
	for (int i=0; i<PrescribedBCs(); ++i)
	{
		FEPrescribedBC& dc = *PrescribedBC(i);
		if (dc.IsActive()) dc.Activate();
	}

	// model loads
	for (int i=0; i<ModelLoads(); ++i)
	{
		FEModelLoad& FC = *ModelLoad(i);
		if (FC.IsActive()) FC.Activate();
	}

	// nonlinear constraints
	for (int i=0; i<NonlinearConstraints(); ++i)
	{
		FENLConstraint* plc = NonlinearConstraint(i);
		if (plc->IsActive()) plc->Activate();
	}

	// contact interfaces
	for (int i=0; i<SurfacePairInteractions(); ++i)
	{
		FESurfacePairInteraction& ci = *SurfacePairInteraction(i);
		if (ci.IsActive()) ci.Activate();
	}

	// activate rigid components
	if (m_imp->m_prs) m_imp->m_prs->Activate();

	// activate linear constraints
	if (m_imp->m_LCM) m_imp->m_LCM->Activate();
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
	m_imp->m_mesh.Reset();

	// reset object data
	if (m_imp->m_prs && (m_imp->m_prs->Reset() == false)) return false;

	// set up rigid joints
	if (m_imp->m_NLC.size() > 0)
	{
		int NC = (int)m_imp->m_NLC.size();
		for (int i=0; i<NC; ++i)
		{
			FENLConstraint* plc = m_imp->m_NLC[i];
			plc->Reset();
		}
	}

	// set the start time
	m_imp->m_ftime = 0;
	m_imp->m_ftime0 = 0;

	// set first time step
	m_imp->m_pStep = m_imp->m_Step[0];
	m_imp->m_nStep = 0;
	for (int i=0; i<Steps(); ++i) GetStep(i)->Reset();

	// reset contact data
	// TODO: I just call Init which I think is okay
	InitContact();

	// Call Activate() to activate all permanent BC's
	Activate();

	return true;
}

//-----------------------------------------------------------------------------
//! Get the current time information.
FETimeInfo FEModel::GetTime()
{
	return FETimeInfo(m_imp->m_ftime, GetCurrentStep()->m_dt);
}

//-----------------------------------------------------------------------------
double FEModel::GetStartTime() const { return m_imp->m_ftime0; }

//-----------------------------------------------------------------------------
void FEModel::SetStartTime(double t) { m_imp->m_ftime0 = t; }

//-----------------------------------------------------------------------------
double FEModel::GetCurrentTime() const { return m_imp->m_ftime; }

//-----------------------------------------------------------------------------
void FEModel::SetCurrentTime(double t) { m_imp->m_ftime = t; }

//=============================================================================
//    P A R A M E T E R   F U N C T I O N S
//=============================================================================

//-----------------------------------------------------------------------------
// helper functions for accessing components of parameters via parameter strings
FEParamValue GetParameterComponent(const ParamString& paramName, FEParam* param)
{
	// make sure we have something to do
	if (param == 0) return FEParamValue();

	if (param->type() == FE_PARAM_DOUBLE)
	{
		int lc = param->GetLoadCurve();
		if (lc == -1) return param->paramValue();
		else return param->GetScale();
	}
	else if (param->type() == FE_PARAM_VEC3D)
	{
		vec3d* v = param->pvalue<vec3d>(0);
		assert(v);
		if (v)
		{
			if      (paramName == "x") return FEParamValue(v->x);
			else if (paramName == "y") return FEParamValue(v->y);
			else if (paramName == "z") return FEParamValue(v->z);
			else return FEParamValue();
		}
		else return FEParamValue();
	}
	else if (param->type() == FE_PARAM_BOOL)
	{
		return param->paramValue();
	}
	else if (param->type() == FE_PARAM_STD_VECTOR_DOUBLE)
	{
		vector<double>& data = param->value< vector<double> >();
		int index = paramName.Index();
		if ((index >= 0) && (index < data.size()))
		return FEParamValue(data[index]);
	}

	return FEParamValue();
}

//-----------------------------------------------------------------------------
//! Return a pointer to the named variable
//! This function returns a pointer to a named variable.

FEParamValue FEModel::FindParameter(const ParamString& paramString)
{
	// make sure it starts with the name of this model
	if (paramString != GetName()) return FEParamValue();

	// see what the next reference is
	ParamString next = paramString.next();

	FEParam* param = GetParameter(next);
	if (param)
	{
		ParamString paramComp = next.last();
		return GetParameterComponent(paramComp, param);
	}

	// if we get here, handle some special cases
	if (next == "mesh")
	{
		ParamString nodeString = next.next();
		if (nodeString == "node")
		{
			FEMesh& mesh = GetMesh();
			FENode* node = 0;
			int nid = nodeString.Index();
			if (nid < 0)
			{
				ParamString fnc = nodeString.next();
				if (fnc == "fromId")
				{
					nid = fnc.Index();
					node = mesh.FindNodeFromID(nid);
					nodeString = nodeString.next();
				}
				else return FEParamValue();
			}
			else if ((nid >=0) && (nid < mesh.Nodes()))
			{
				node = &mesh.Node(nid);
			}

			if (node)
			{
				ParamString paramString = nodeString.next();
				if (paramString == "position")
				{
					vec3d& rt = node->m_rt;
					ParamString c = paramString.next();
					if (c == "x") return FEParamValue(rt.x);
					if (c == "y") return FEParamValue(rt.y);
					if (c == "z") return FEParamValue(rt.z);
					return FEParamValue();
				}
				else return FEParamValue();
			}
			else return FEParamValue();
		}
		else return FEParamValue();
	}

	if ((next == "rigidbody") && m_imp->m_prs)
	{
		FEMaterial* mat = 0;
		if (next.IDString()) mat = FindMaterial(next.IDString());
		if ((mat != 0) && (mat->IsRigid()))
		{
			ParamString paramName = next.next();

			// the rigid bodies are dealt with differently
			int nmat = mat->GetID() - 1;
			int NRB = m_imp->m_prs->Objects();
			for (int i = 0; i<NRB; ++i)
			{
				FERigidBody* ob = m_imp->m_prs->Object(i);
				if (ob && (ob->GetMaterialID() == nmat))
				{
					FEParam* pp = ob->GetParameter(paramName);
					return GetParameterComponent(paramName.last(), pp);
				}
			}
		}
	}

	// oh, oh, we didn't find it
	return FEParamValue();
}

//-----------------------------------------------------------------------------
FECoreBase* FEModel::FindComponent(const ParamString& prop)
{
	// make sure it starts with the name of this model
	if (prop != GetName()) return 0;

	// see what the next reference is
	ParamString next = prop.next();

	// next, find the property
	FECoreBase* pc = GetProperty(next);

	return pc;
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

	// prescribed displacements
	for (int i=0; i<PrescribedBCs(); ++i)
	{
		FEParameterList& pl = PrescribedBC(i)->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
	}

	// evaluate nodal loads
	for (int i=0; i<NodalLoads(); ++i)
	{
		FEParameterList& pl = NodalLoad(i)->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
	}

	// evaluate edge load parameter lists
	for (int i=0; i<EdgeLoads(); ++i)
	{
		FEParameterList& pl = EdgeLoad(i)->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
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
	for (int i=0; i<ModelLoads(); ++i)
	{
		FEParameterList& pl = ModelLoad(i)->GetParameterList();
		if (EvaluateParameterList(pl) == false) return false;
	}

	// give the rigid system a chance
	if (m_imp->m_prs->EvaluateParameterLists() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Evaluate a parameter list
bool FEModel::EvaluateParameterList(FEParameterList &pl)
{
	const int NLC = LoadCurves();
	FEParamIterator pi = pl.first();
	for (int j=0; j<pl.Parameters(); ++j, ++pi)
	{
		int nlc = pi->GetLoadCurve();
		if (nlc >= 0)
		{
			if (nlc >= NLC) return fecore_error("Invalid load curve ID");

			double v = GetLoadCurve(nlc)->Value();
			switch (pi->type())
			{
			case FE_PARAM_INT   : pi->value<int>() = (int) v; break;
			case FE_PARAM_DOUBLE: pi->value<double>() = pi->GetScaleDouble()*v; break;
			case FE_PARAM_BOOL  : pi->value<bool>() = (v > 0? true : false); break;
			case FE_PARAM_VEC3D : pi->value<vec3d>() = pi->GetScaleVec3d()*v; break;
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
DOFS& FEModel::GetDOFS()
{
	return m_imp->m_dofs;
}

//-----------------------------------------------------------------------------
int FEModel::GetDOFIndex(const char* sz)
{
	return m_imp->m_dofs.GetDOF(sz);
}

//-----------------------------------------------------------------------------
int FEModel::GetDOFIndex(const char* szvar, int n)
{
	return m_imp->m_dofs.GetDOF(szvar, n);
}

//-----------------------------------------------------------------------------
// Call the callback function if there is one defined
//
bool FEModel::DoCallback(unsigned int nevent)
{
	try
	{
		// do the callbacks
		bool bret = CallbackHandler::DoCallback(this, nevent);
		return bret;
	}
	catch (ExitRequest)
	{
		if (nevent == CB_MAJOR_ITERS) return false;
		else throw;
	}
	catch (ForceConversion)
	{
		throw;
	}
	catch (IterationFailure)
	{
		throw;
	}
	catch (...)
	{
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEModel::SetGlobalConstant(const string& s, double v)
{
	m_imp->m_Const[s] = v;
	return;
}

//-----------------------------------------------------------------------------
double FEModel::GetGlobalConstant(const string& s)
{
	return (m_imp->m_Const.count(s) ? m_imp->m_Const.find(s)->second : 0);
}

//-----------------------------------------------------------------------------
void FEModel::AddGlobalData(FEGlobalData* psd)
{
	m_imp->m_GD.push_back(psd);
}

//-----------------------------------------------------------------------------
FEGlobalData* FEModel::GetGlobalData(int i)
{
	return m_imp->m_GD[i];
}

//-----------------------------------------------------------------------------
int FEModel::GlobalDataItems()
{
	return (int)m_imp->m_GD.size();
}

//-----------------------------------------------------------------------------
void FEModel::AddModelData(FEModelData* data)
{
	m_imp->m_Data.AddProperty(data);
}

//-----------------------------------------------------------------------------
FEModelData* FEModel::GetModelData(int i)
{
	return m_imp->m_Data[i];
}

//-----------------------------------------------------------------------------
int FEModel::ModelDataItems() const
{
	return (int) m_imp->m_Data.size();
}

//-----------------------------------------------------------------------------
void FEModel::UpdateModelData()
{
	for (int i=0; i<ModelDataItems(); ++i)
	{
		FEModelData* data = GetModelData(i);
		data->Update();
	}
}

//-----------------------------------------------------------------------------
//! Set the title of the model
void FEModel::SetTitle(const char* sz)
{ 
	strcpy(m_imp->m_sztitle, sz);
}

//-----------------------------------------------------------------------------
//! Return the title of the model
const char* FEModel::GetTitle()
{ 
	return m_imp->m_sztitle;
}

//-----------------------------------------------------------------------------
//! Find a BC based on its ID. This is needed for restarts.
FEModelComponent* FEModel::FindModelComponent(int nid)
{
	int i;
	for (i=0; i<(int) m_imp->m_BC.size (); ++i) if (m_imp->m_BC [i]->GetClassID() == nid) return m_imp->m_BC [i];
	for (i=0; i<(int) m_imp->m_DC.size (); ++i) if (m_imp->m_DC [i]->GetClassID() == nid) return m_imp->m_DC [i];
	for (i=0; i<(int) m_imp->m_IC.size (); ++i) if (m_imp->m_IC [i]->GetClassID() == nid) return m_imp->m_IC [i];
	for (i=0; i<(int) m_imp->m_FC.size (); ++i) if (m_imp->m_FC [i]->GetClassID() == nid) return m_imp->m_FC [i];
	for (i=0; i<(int) m_imp->m_SL.size (); ++i) if (m_imp->m_SL [i]->GetClassID() == nid) return m_imp->m_SL [i];
	for (i=0; i<(int) m_imp->m_ML.size (); ++i) if (m_imp->m_ML [i]->GetClassID() == nid) return m_imp->m_ML [i];
	for (i=0; i<(int) m_imp->m_CI.size (); ++i) if (m_imp->m_CI [i]->GetClassID() == nid) return m_imp->m_CI [i];
	for (i=0; i<(int) m_imp->m_NLC.size(); ++i) if (m_imp->m_NLC[i]->GetClassID() == nid) return m_imp->m_NLC[i];
	return (m_imp->m_prs ? m_imp->m_prs->FindModelComponent(nid) : 0);
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
	m_imp->m_linearSolver = fem.m_imp->m_linearSolver;
	m_imp->m_bwopt = fem.m_imp->m_bwopt;
	m_imp->m_nStep = fem.m_imp->m_nStep;
	m_imp->m_ftime = fem.m_imp->m_ftime;
	m_imp->m_ftime0 = fem.m_imp->m_ftime0;
	m_ut4_alpha = fem.m_ut4_alpha;
	m_ut4_bdev = fem.m_ut4_bdev;
	m_udghex_hg = fem.m_udghex_hg;
	m_imp->m_pStep = 0;

	// --- Steps ---

	// copy the steps
	// NOTE: This does not copy the boundary conditions of the steps
	int NSTEP = fem.Steps();
	for (int i=0; i<NSTEP; ++i)
	{
		// get the type string
		FEAnalysis* ps = fem.GetStep(i);

		// create a new step
		FEAnalysis* pnew = new FEAnalysis(this);

		// copy additional info
		pnew->m_nanalysis = ps->m_nanalysis;

		pnew->m_ntime		= ps->m_ntime;
		pnew->m_final_time	= ps->m_final_time;
		pnew->m_dt			= ps->m_dt;
		pnew->m_dt0			= ps->m_dt0;
		pnew->m_dtp         = ps->m_dtp;
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
		FESolver* psolver = ps->GetFESolver();
		const char* sztype = psolver->GetTypeStr();

		// create a new solver
		FESolver* pnew_solver = fecore_new<FESolver>(FESOLVER_ID, sztype, this);
		assert(pnew_solver);
		pnew->SetFESolver(pnew_solver);

		// copy parameters
		pnew_solver->GetParameterList() = psolver->GetParameterList();

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
	assert(m_imp->m_MAT.size() == fem.m_imp->m_MAT.size());

	// --- Mesh ---

	// copy the mesh data
	// NOTE: This will not assign materials to the new domains
	// A. copy nodes
	FEMesh& sourceMesh = fem.GetMesh();
	FEMesh& mesh = GetMesh();
	int N = sourceMesh.Nodes();
	mesh.CreateNodes(N);
	for (int i=0; i<N; ++i)
	{
		mesh.Node(i) = sourceMesh.Node(i);
	}

	// B. domains
	// let's first create a table of material indices for the old domains
	int NDOM = sourceMesh.Domains();
	vector<int> LUT(NDOM);
	for (int i=0; i<NDOM; ++i)
	{
		FEMaterial* pm = sourceMesh.Domain(i).GetMaterial();
		for (int j=0; j<NMAT; ++j)
		{
			if (pm == fem.GetMaterial(j))
			{
				LUT[i] = j;
				break;
			}
		}
	}

	// now allocate domains
	for (int i=0; i<NDOM; ++i)
	{
		FEDomain& dom = sourceMesh.Domain(i);
		const char* sz = dom.GetTypeStr();

		// create a new domain
		FEDomain* pd = fecore_new<FEDomain>(FEDOMAIN_ID, sz, this);
		assert(pd);
		pd->SetMaterial(GetMaterial(LUT[i]));

		// copy domain data
		pd->CopyFrom(&dom);

		// add it to the mesh
		mesh.AddDomain(pd);
	}

	// --- boundary conditions ---

	int NDC = fem.PrescribedBCs();
	for (int i=0; i<NDC; ++i)
	{
		FEPrescribedBC* pbc = fem.PrescribedBC(i);
		const char* sz = pbc->GetTypeStr();

		FEPrescribedBC* pnew = fecore_new<FEPrescribedBC>(FEBC_ID, sz, this);
		assert(pnew);

		pnew->CopyFrom(pbc);

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
		mesh.AddSurface(pnew->GetMasterSurface());
		mesh.AddSurface(pnew->GetSlaveSurface ());
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
		if (ps) mesh.AddSurface(ps);
	}

	// --- Load curves ---
	// copy load curves
	int NLD = fem.LoadCurves();
	for (int i = 0; i<NLD; ++i)
	{
		FELoadCurve* lc = fem.GetLoadCurve(i);
		FELoadCurve* newlc = fecore_new<FELoadCurve>(FELOADCURVE_ID, lc->GetTypeStr(), this);
		if (newlc == 0)
		{
			// we can get here if the load curve was not created with fecore_new
			// do it the hard way
			if (dynamic_cast<FELinearRamp*>(lc)) newlc = new FELinearRamp(this);
			else if (dynamic_cast<FEDataLoadCurve*>(lc)) newlc = new FEDataLoadCurve(this);
		}
		bool b = newlc->CopyFrom(lc); assert(b);
		AddLoadCurve(newlc);
	}

	// copy linear constraints
	if (fem.m_imp->m_LCM)
	{
		if (m_imp->m_LCM) delete m_imp->m_LCM;
		m_imp->m_LCM = new FELinearConstraintManager(this);
		m_imp->m_LCM->CopyFrom(*fem.m_imp->m_LCM);
	}

	// TODO: copy all the properties
//	assert(false);
}

//-----------------------------------------------------------------------------
// This function serializes data to a stream.
// This is used for running and cold restarts.
void FEModel::Serialize(DumpStream& ar)
{
	if (ar.IsShallow())
	{
		// stream model data
		if (ar.IsSaving())
		{
			ar << m_imp->m_ftime;
		}
		else
		{
			ar >> m_imp->m_ftime;
		}
		ar.check();

		// stream mesh
		m_imp->m_mesh.Serialize(ar);
		ar.check();

		// stream rigid body data
		if (m_imp->m_prs) m_imp->m_prs->Serialize(ar);
		ar.check();

		// stream contact data
		for (int i = 0; i<SurfacePairInteractions(); ++i) m_imp->m_CI[i]->Serialize(ar);
		ar.check();

		// stream nonlinear constraints
		for (int i = 0; i<NonlinearConstraints(); ++i) m_imp->m_NLC[i]->Serialize(ar);
		ar.check();
	}
	else
	{
		if (ar.IsSaving() == false) Clear();

		if (ar.IsSaving())
		{
			ar << m_imp->m_sztitle;
		}
		else
		{
			ar >> m_imp->m_sztitle;
		}

		m_imp->m_dofs.Serialize(ar);
		m_imp->SerializeLoadData(ar);
		m_imp->SerializeGlobals(ar);
		m_imp->SerializeMaterials(ar);
		m_imp->SerializeGeometry(ar);
		m_imp->SerializeContactData(ar);
		m_imp->SerializeBoundaryData(ar);
		m_imp->SerializeAnalysisData(ar);
	}
}

//-----------------------------------------------------------------------------
//! Serialize load curves
void FEModel::Implementation::SerializeLoadData(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		// save curve data
		ar << m_fem->LoadCurves();
		for (int i = 0; i<m_fem->LoadCurves(); ++i)
		{
			FELoadCurve* lc = m_fem->GetLoadCurve(i);
			ar << lc->GetTypeStr();
			
			lc->Serialize(ar);
		}
	}
	else
	{
		// loadcurve data
		char szlc[256] = { 0 };
		int nlc;
		ar >> nlc;
		m_LC.Clear();
		for (int i=0; i<nlc; ++i)
		{
			ar >> szlc;
			FELoadCurve* plc = fecore_new<FELoadCurve>(FELOADCURVE_ID, szlc, m_fem);
			plc->Serialize(ar);
			m_fem->AddLoadCurve(plc);
		}
	}
}

//-----------------------------------------------------------------------------
//! Serialize global data
void FEModel::Implementation::SerializeGlobals(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		int NC = (int)m_Const.size();
		ar << NC;
		if (NC > 0)
		{
			char sz[256] = {0};
			map<string, double>::iterator it;
			for (it = m_Const.begin(); it != m_Const.end(); ++it)
			{
				strcpy(sz, it->first.c_str());
				ar << sz;
				ar << it->second;
			}
		}
		int nGD = m_fem->GlobalDataItems();
		ar << nGD;
		for (int i=0; i<nGD; i++) 
		{
			FEGlobalData* pgd = m_fem->GetGlobalData(i);
			ar << pgd->GetTypeStr();
			pgd->Serialize(ar);
		}
	}
	else
	{
		char sz[256] = {0};
		double v;
		int NC;
		ar >> NC;
		m_Const.clear();
		for (int i=0; i<NC; ++i)
		{
			ar >> sz >> v;
			m_fem->SetGlobalConstant(string(sz), v);
		}
		int nGD;
		ar >> nGD;
		if (nGD) 
		{
			char sztype[256];
			for (int i=0; i<nGD; ++i)
			{
				ar >> sztype;
				FEGlobalData* pgd = fecore_new<FEGlobalData>(FEGLOBALDATA_ID, sztype, m_fem);
				pgd->Serialize(ar);
				m_fem->AddGlobalData(pgd);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! serialize material data
void FEModel::Implementation::SerializeMaterials(DumpStream& ar)
{
	FECoreKernel& febio = FECoreKernel::GetInstance();

	if (ar.IsSaving())
	{
		// store the nr of materials
		ar << m_fem->Materials();

		// store the materials
		for (int i = 0; i<m_fem->Materials(); ++i)
		{
			FEMaterial* pmat = m_fem->GetMaterial(i);

			// store the type string
			ar << pmat->GetTypeStr();

			// store the name
			ar << pmat->GetName();

			// store material parameters
			pmat->Serialize(ar);
		}
	}
	else
	{
		// read the number of materials
		int nmat;
		ar >> nmat;

		// read the material data
		char szmat[256] = {0}, szvar[256] = {0};
		for (int i=0; i<nmat; ++i)
		{
			// read the type string
			ar >> szmat;

			// create a material
			FEMaterial* pmat = fecore_new<FEMaterial>(FEMATERIAL_ID, szmat, m_fem);
			assert(pmat);

			// read the name
			ar >> szmat;
			pmat->SetName(szmat);

			// read all parameters
			pmat->Serialize(ar);

			// Add material and parameter list to FEM
			m_fem->AddMaterial(pmat);
		}
	}
}

//-----------------------------------------------------------------------------
void FEModel::Implementation::SerializeGeometry(DumpStream &ar)
{
	// serialize the mesh first 
	m_mesh.Serialize(ar);

	// serialize the rigid system
	if (m_prs) m_prs->Serialize(ar);
}

//-----------------------------------------------------------------------------
//! serialize contact data
void FEModel::Implementation::SerializeContactData(DumpStream &ar)
{
	FECoreKernel& febio = FECoreKernel::GetInstance();

	if (ar.IsSaving())
	{
		ar << m_fem->SurfacePairInteractions();
		for (int i = 0; i<m_fem->SurfacePairInteractions(); ++i)
		{
			FESurfacePairInteraction* pci = m_fem->SurfacePairInteraction(i);

			// store the type string
			ar << pci->GetTypeStr();

			pci->Serialize(ar);
		}
	}
	else
	{
		int numci;
		ar >> numci;

		char szci[256] = {0};
		for (int i=0; i<numci; ++i)
		{
			// get the interface type
			ar >> szci;

			// create a new interface
			FESurfacePairInteraction* pci = fecore_new<FESurfacePairInteraction>(FESURFACEPAIRINTERACTION_ID, szci, m_fem);

			// serialize interface data from archive
			pci->Serialize(ar);

			// add interface to list
			m_fem->AddSurfacePairInteraction(pci);

			// add surfaces to mesh
			FEMesh& m = m_mesh;
			if (pci->GetMasterSurface()) m.AddSurface(pci->GetMasterSurface());
			m.AddSurface(pci->GetSlaveSurface());
		}	
	}
}

//-----------------------------------------------------------------------------
//! \todo Do we need to store the m_bActive flag of the boundary conditions?
void FEModel::Implementation::SerializeBoundaryData(DumpStream& ar)
{
	FECoreKernel& febio = FECoreKernel::GetInstance();

	if (ar.IsSaving())
	{
		// fixed bc's
		ar << (int)m_BC.size();
		for (int i = 0; i<(int)m_BC.size(); ++i)
		{
			FEFixedBC& bc = *m_BC[i];
			bc.Serialize(ar);
		}

		// displacements
		ar << (int)m_DC.size();
		for (int i = 0; i<(int)m_DC.size(); ++i)
		{
			FEPrescribedBC& dc = *m_DC[i];
			dc.Serialize(ar);
		}

		// initial conditions
		ar << (int)m_IC.size();
		for (int i = 0; i<(int)m_IC.size(); ++i)
		{
			FEInitialCondition& ic = *m_IC[i];
			ar << ic.GetTypeStr();
			ic.Serialize(ar);
		}

		// nodal loads
		ar << (int)m_FC.size();
		for (int i = 0; i<(int)m_FC.size(); ++i)
		{
			FENodalLoad& fc = *m_FC[i];
			fc.Serialize(ar);
		}

		// surface loads
		ar << (int)m_SL.size();
		for (int i = 0; i<(int)m_SL.size(); ++i)
		{
			FESurfaceLoad* psl = m_SL[i];

			// get the surface
			FESurface& s = psl->GetSurface();
			s.Serialize(ar);

			// save the load data
			ar << psl->GetTypeStr();
			psl->Serialize(ar);
		}

		// edge loads
		ar << (int)m_EL.size();
		for (int i = 0; i<(int)m_EL.size(); ++i)
		{
			FEEdgeLoad* pel = m_EL[i];

			// get the edge
			FEEdge& e = pel->Edge();
			e.Serialize(ar);

			// save the load data
			ar << pel->GetTypeStr();
			pel->Serialize(ar);
		}

		// body loads
		ar << (int)m_BL.size();
		for (int i = 0; i<(int)m_BL.size(); ++i)
		{
			FEBodyLoad* pbl = m_BL[i];
			ar << pbl->GetTypeStr();
			pbl->Serialize(ar);
		}

		// model loads
		ar << (int)m_ML.size();
		for (int i = 0; i<(int)m_ML.size(); ++i)
		{
			FEModelLoad& ml = *m_ML[i];
			ar << ml.GetTypeStr();
			ml.Serialize(ar);
		}

		// nonlinear constraints
		int n = (int)m_NLC.size();
		ar << n;
		if (n) 
		{
			for (int i=0; i<n; ++i) 
			{
				FENLConstraint& ci = *m_NLC[i];
				ar << ci.GetTypeStr();
				ci.Serialize(ar);
			}
		}
	}
	else
	{
		int n;
		char sz[256] = {0};

		// fixed bc's
		// NOTE: I think this may create a memory leak if m_BC is not empty
		ar >> n;
		m_BC.Clear();
		for (int i=0; i<n; ++i) 
		{
			FEFixedBC* pbc = new FEFixedBC(m_fem);
			pbc->Serialize(ar);
			m_fem->AddFixedBC(pbc);
		}

		// displacements
		ar >> n;
		m_DC.Clear();
		for (int i=0; i<n; ++i) 
		{
			FEPrescribedDOF* pdc = fecore_new<FEPrescribedDOF>(FEBC_ID, "prescribe", m_fem);
			pdc->Serialize(ar);
			m_fem->AddPrescribedBC(pdc);
		}

		// initial conditions
		ar >> n;
		m_IC.Clear();
		for (int i=0; i<n; ++i) 
		{
			ar >> sz;
			FEInitialCondition* pic = fecore_new<FEInitialCondition>(FEIC_ID, sz, m_fem);
			assert(pic);
			pic->Serialize(ar);
			m_fem->AddInitialCondition(pic);
		}

		// nodal loads
		ar >> n;
		m_FC.Clear();
		for (int i=0; i<n; ++i)
		{
			FENodalLoad* pfc = new FENodalLoad(m_fem);
			pfc->Serialize(ar);
			m_fem->AddNodalLoad(pfc);
		}

		// surface loads
		ar >> n;
		m_SL.Clear();
		for (int i=0; i<n; ++i)
		{
			// create a new surface
			FESurface* psurf = new FESurface(&m_mesh);
			psurf->Serialize(ar);

			// read load data
			char sztype[256] = {0};
			ar >> sztype;
			FESurfaceLoad* ps = fecore_new<FESurfaceLoad>(FESURFACELOAD_ID, sztype, m_fem);
			assert(ps);
			ps->SetSurface(psurf);

			ps->Serialize(ar);

			m_SL.AddProperty(ps);
			m_mesh.AddSurface(psurf);
		}

		// edge loads
		ar >> n;
		m_EL.Clear();
		for (int i=0; i<n; ++i)
		{
			// create a new edge
			FEEdge* pedge = new FEEdge(&m_mesh);
			pedge->Serialize(ar);

			// read load data
			char sztype[256] = {0};
			ar >> sztype;
			FEEdgeLoad* pel = fecore_new<FEEdgeLoad>(FEEDGELOAD_ID, sztype, m_fem);
			assert(pel);
			pel->SetEdge(pedge);

			pel->Serialize(ar);

			m_EL.AddProperty(pel);
			m_mesh.AddEdge(pedge);
		}

		// body loads
		int nbl;
		ar >> nbl;
		m_BL.Clear();
		char szbl[256] = {0};
		for (int i=0; i<nbl; ++i)
		{
			ar >> szbl;
			FEBodyLoad* pbl = fecore_new<FEBodyLoad>(FEBODYLOAD_ID, szbl, m_fem);
			assert(pbl);

			pbl->Serialize(ar);
			m_BL.AddProperty(pbl);
		}

		// model loads
		ar >> n;
		m_ML.Clear();
		for (int i=0; i<n; ++i)
		{
			// read load data
			char sztype[256] = {0};
			ar >> sztype;
			FEModelLoad* pml = fecore_new<FEModelLoad>(FEBC_ID, sztype, m_fem);
			assert(pml);

			pml->Serialize(ar);
			m_fem->AddModelLoad(pml);
		}

		// non-linear constraints
		ar >> n;
		m_NLC.Clear();
		for (int i=0; i<n; ++i)
		{
			char sztype[256] = { 0 };
			ar >> sztype;
			FENLConstraint* pc = fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, sztype, m_fem);
			assert(pc);

			pc->Serialize(ar);
			m_fem->AddNonlinearConstraint(pc);
		}
	}

	// serialize rigid stuff
	if (m_prs) m_prs->Serialize(ar);

	// serialize linear constraints
	if (m_LCM) m_LCM->Serialize(ar);
}

//-----------------------------------------------------------------------------
//! Serialize analysis data
void FEModel::Implementation::SerializeAnalysisData(DumpStream &ar)
{
	if (ar.IsSaving())
	{
		// analysis steps
		ar << (int)m_Step.size();
		for (int i = 0; i<(int)m_Step.size(); ++i)
		{
			m_Step[i]->Serialize(ar);
		}

		ar << m_nStep;
		ar << m_ftime << m_ftime0;

		// direct solver data
		ar << m_linearSolver;
		ar << m_bwopt;
	}
	else
	{
		m_Step.Clear();

		char sztype[256] = {0};

		// analysis steps
		int nsteps;
		ar >> nsteps;
		for (int i=0; i<nsteps; ++i)
		{
			FEAnalysis* pstep = new FEAnalysis(m_fem); assert(pstep);
			pstep->Serialize(ar);
			m_fem->AddStep(pstep);
		}
		ar >> m_nStep;
		ar >> m_ftime >> m_ftime0;

		// direct solver data
		ar >> m_linearSolver;
		ar >> m_bwopt;

		// set the correct step
		m_pStep = m_Step[m_nStep];
	}
}

//-----------------------------------------------------------------------------
void FEModel::BuildMatrixProfile(FEGlobalMatrix& G, bool breset)
{
	FEAnalysis* pstep = GetCurrentStep();
	FEMesh& mesh = GetMesh();
	FERigidSystem& rigid = *GetRigidSystem();
    DOFS& fedofs = GetDOFS();
    int MAX_NDOFS = fedofs.GetTotalDOFS();

	// when reset is true we build the entire matrix profile
	// (otherwise we only build the "dynamic" profile)
	if (breset)
	{
		vector<int> elm;

		// Add all elements to the profile
		// Loop over all active domains
		for (int nd=0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& d = mesh.Domain(nd);
			d.BuildMatrixProfile(G);
		}

		// Add rigid bodies to the profile
		rigid.BuildMatrixProfile(G);

		// linear constraints
		if (m_imp->m_LCM) m_imp->m_LCM->BuildMatrixProfile(G);
	}
	else
	{
		// Do the "dynamic" profile. That is the part of the profile that always changes
		// This is mostly contact
		// do the nonlinear constraints
		int M = NonlinearConstraints();
		for (int m=0; m<M; ++m)
		{
			FENLConstraint* pnlc = NonlinearConstraint(m);
			if (pnlc->IsActive()) pnlc->BuildMatrixProfile(G);
		}	

		// All following "elements" are nonstatic. That is, they can change
		// connectivity between calls to this function. All of these elements
		// are related to contact analysis (at this point).
		if (SurfacePairInteractions() > 0)
		{
			// Add all contact interface elements
			for (int i=0; i<SurfacePairInteractions(); ++i)
			{
				FESurfacePairInteraction* pci = SurfacePairInteraction(i);
				if (pci->IsActive()) pci->BuildMatrixProfile(G);
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEModel::GetNodeData(int ndof, vector<double>& data)
{
	// get the dofs
	DOFS& dofs = GetDOFS();

	// make sure the dof index is valid
	if ((ndof < 0) || (ndof >= dofs.GetTotalDOFS())) return false;

	// get the mesh and number of nodes
	FEMesh& mesh = GetMesh();
	int N = mesh.Nodes();

	// make sure data is correct size
	data.resize(N, 0.0);

	// loop over all nodes
	for (int i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		data[i] = node.get(ndof);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEModel::ClearDataArrays()
{
	// clear the surface maps
	for (int i = 0; i<(int)m_imp->m_DataArray.size(); ++i) delete m_imp->m_DataArray[i].second;
	m_imp->m_DataArray.clear();
}

//-----------------------------------------------------------------------------
void FEModel::AddDataArray(const char* szname, FEDataArray* map)
{
	m_imp->m_DataArray.push_back(pair<string, FEDataArray*>(string(szname), map));
}

//-----------------------------------------------------------------------------
FEDataArray* FEModel::FindDataArray(const char* szmap)
{
	string name(szmap);
	for (int i = 0; i<(int)m_imp->m_DataArray.size(); ++i)
	{
		if (m_imp->m_DataArray[i].first == name) return m_imp->m_DataArray[i].second;
	}
	return 0;
}
