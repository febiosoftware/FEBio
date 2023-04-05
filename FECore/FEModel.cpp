/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEModel.h"
#include "FELoadController.h"
#include "FEMaterial.h"
#include "FEModelLoad.h"
#include "FEBoundaryCondition.h"
#include "FENodalLoad.h"
#include "FESurfaceLoad.h"
#include "FEEdgeLoad.h"
#include "FEBodyLoad.h"
#include "FEInitialCondition.h"
#include "FESurfacePairConstraint.h"
#include "FENLConstraint.h"
#include "FEAnalysis.h"
#include "FEGlobalData.h"
#include "FECoreKernel.h"
#include "FELinearConstraintManager.h"
#include "log.h"
#include "FEDataArray.h"
#include "FESurfaceConstraint.h"
#include "FEModelParam.h"
#include "FEShellDomain.h"
#include "FEEdge.h"
#include "FEMeshAdaptor.h"
#include <string>
#include <map>
#include "DumpStream.h"
#include "LinearSolver.h"
#include "FETimeStepController.h"
#include "Timer.h"
#include "DumpMemStream.h"
#include "FEPlotDataStore.h"
#include "FESolidDomain.h"
#include "FEShellDomain.h"
#include "FETrussDomain.h"
#include "FEDomain2D.h"
#include "FEDiscreteDomain.h"
#include "FEDataGenerator.h"
#include "FEModule.h"
#include <stdarg.h>
using namespace std;

//-----------------------------------------------------------------------------
// Implementation class for the FEModel class
class FEModel::Implementation
{
public:
	struct LoadParam
	{
		FEParam*			param;
		int					lc;

		double		m_scl;
		vec3d		m_vscl;

		LoadParam()
		{
			m_scl = 1.0;
			m_vscl = vec3d(0, 0, 0);
		}

		void Serialize(DumpStream& ar)
		{
			ar & lc;
			ar & m_scl & m_vscl;

			if (ar.IsShallow() == false)
			{
				// we can't save the FEParam* directly, so we need to store meta data and try to find it on loading
				if (ar.IsSaving())
				{
					FECoreBase* pc = dynamic_cast<FECoreBase*>(param->parent()); assert(pc);
					ar << pc;
					ar << param->name();
				}
				else
				{
					FECoreBase* pc = nullptr;
					ar >> pc; assert(pc);
					
					char name[256] = { 0 };
					ar >> name;

					param = pc->FindParameter(name); assert(param);
				}
			}
			else param = nullptr;
		}
	};

public:
	Implementation(FEModel* fem) : m_fem(fem), m_mesh(fem), m_dmp(*fem)
	{
		// --- Analysis Data ---
		m_pStep = 0;
		m_nStep = -1;
		m_ftime0 = 0;

		m_nupdates = 0;

		m_bsolved = false;

		m_block_log = false;

		m_printParams = false;

		m_meshUpdate = false;

		// create the linear constraint manager
		m_LCM = new FELinearConstraintManager(fem);

		// allocate timers
		// Make sure enough timers are allocated for all the TimerIds!
		m_timers.resize(7);
	}

	void Serialize(DumpStream& ar);

	void PushState()
	{
		DumpMemStream& ar = m_dmp;
		ar.clear(); // this also prepares the stream for writing
		m_fem->Serialize(ar);
	}

	bool PopState()
	{
		// get the dump stream
		DumpMemStream& ar = m_dmp;

		// make sure we have data to rewind
		if (ar.size() == 0) return false;

		// prepare the archive for reading
		ar.Open(false, true);

		// restore the previous state
		m_fem->Serialize(m_dmp);

		return true;
	}

public: // TODO: Find a better place for these parameters
	FETimeInfo	m_timeInfo;			//!< current time value
	double		m_ftime0;			//!< start time of current step

	bool	m_block_log;

	int		m_nupdates;	//!< number of calls to FEModel::Update

public:
	std::vector<FEMaterial*>				m_MAT;		//!< array of materials
	std::vector<FEBoundaryCondition*>		m_BC;		//!< boundary conditions
	std::vector<FEModelLoad*>				m_ML;		//!< model loads
	std::vector<FEInitialCondition*>		m_IC;		//!< initial conditions
    std::vector<FESurfacePairConstraint*>   m_CI;       //!< contact interface array
	std::vector<FENLConstraint*>			m_NLC;		//!< nonlinear constraints
	std::vector<FELoadController*>			m_LC;		//!< load controller data
	std::vector<FEAnalysis*>				m_Step;		//!< array of analysis steps
	std::vector<FEMeshAdaptor*>				m_MA;		//!< mesh adaptors
	std::vector<FEMeshDataGenerator*>		m_MD;		//!< mesh data generators

	std::vector<LoadParam>		m_Param;	//!< list of parameters controller by load controllers
	std::vector<Timer>			m_timers;	// list of timers

public:
	FEAnalysis*		m_pStep;	//!< pointer to current analysis step
	int				m_nStep;	//!< current index of analysis step
	bool			m_printParams;	//!< print parameters
	bool			m_meshUpdate;	//!< mesh update flag

	std::string		m_units;	// units string

public:
	// The model
	FEModel*	m_fem;

	// module name
	std::string		m_moduleName;

	bool	m_bsolved;	// solved flag

	// DOFS data
	DOFS	m_dofs;				//!< list of degree of freedoms in this model

	// Geometry data
	FEMesh		m_mesh;			//!< the one and only FE mesh

	// linear constraint data
	FELinearConstraintManager*	m_LCM;

	DataStore	m_dataStore;			//!< the data store used for data logging

	FEPlotDataStore	m_plotData;		//!< Output request for plot file

	DumpMemStream	m_dmp;	// only used by incremental solver

public: // Global Data
	std::map<string, double> m_Const;	//!< Global model constants
	vector<FEGlobalData*>	m_GD;		//!< global data structures
	std::vector<FEGlobalVariable*>	m_Var;

	FEMODEL_MEMORY_STATS	m_memstats;
};

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEModel, FECoreBase)

	// model parameters
	ADD_PARAMETER(m_imp->m_timeInfo.currentTime, "time");

	// model properties
	ADD_PROPERTY(m_imp->m_MAT , "material"       );
	ADD_PROPERTY(m_imp->m_BC  , "bc"             );
	ADD_PROPERTY(m_imp->m_ML  , "load"           );
	ADD_PROPERTY(m_imp->m_IC  , "initial"        );
	ADD_PROPERTY(m_imp->m_CI  , "contact"        );
	ADD_PROPERTY(m_imp->m_NLC , "constraint"     );
	ADD_PROPERTY(m_imp->m_MA  , "mesh_adaptor"   );
	ADD_PROPERTY(m_imp->m_LC  , "load_controller");
	ADD_PROPERTY(m_imp->m_MD  , "mesh_data"      );
	ADD_PROPERTY(m_imp->m_Step, "step"           );

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEModel::FEModel(void) : FECoreBase(this), m_imp(new FEModel::Implementation(this))
{
	// set the name
	SetName("fem");

	// reset all timers
	ResetAllTimers();
}

//-----------------------------------------------------------------------------
//! Delete all dynamically allocated data
FEModel::~FEModel(void)
{
	Clear();
}

//-----------------------------------------------------------------------------
//! return the data store
DataStore& FEModel::GetDataStore()
{
	return m_imp->m_dataStore;
}

//-----------------------------------------------------------------------------
FEPlotDataStore& FEModel::GetPlotDataStore() { return m_imp->m_plotData; }

//-----------------------------------------------------------------------------
const FEPlotDataStore& FEModel::GetPlotDataStore() const { return m_imp->m_plotData; }

//-----------------------------------------------------------------------------
//! will return true if the model solved succussfully
bool FEModel::IsSolved() const
{
	return m_imp->m_bsolved;
}

//-----------------------------------------------------------------------------
// call this function to set the mesh's update flag
void FEModel::SetMeshUpdateFlag(bool b)
{
	m_imp->m_meshUpdate = b;
}

//-----------------------------------------------------------------------------
void FEModel::Clear()
{
	// clear dofs
	m_imp->m_dofs.Reset();

	// clear all properties
	for (FEMaterial* mat             : m_imp->m_MAT ) delete  mat; m_imp->m_MAT.clear();
	for (FEBoundaryCondition* bc     : m_imp->m_BC  ) delete   bc; m_imp->m_BC.clear();
	for (FEModelLoad* ml             : m_imp->m_ML  ) delete   ml; m_imp->m_ML.clear();
	for (FEInitialCondition* ic      : m_imp->m_IC  ) delete   ic; m_imp->m_IC.clear();
	for (FESurfacePairConstraint* ci : m_imp->m_CI  ) delete   ci; m_imp->m_CI.clear();
	for (FENLConstraint* nlc         : m_imp->m_NLC ) delete   nlc; m_imp->m_NLC.clear();
	for (FELoadController* lc        : m_imp->m_LC  ) delete   lc; m_imp->m_LC.clear();
	for (FEMeshDataGenerator* md     : m_imp->m_MD  ) delete   md; m_imp->m_MD.clear();
	for (FEAnalysis* step            : m_imp->m_Step) delete step; m_imp->m_Step.clear();

	// global data
	for (size_t i = 0; i<m_imp->m_GD.size(); ++i) delete m_imp->m_GD[i]; m_imp->m_GD.clear();
	m_imp->m_Const.clear();

	// global variables (TODO: Should I delete the corresponding parameters?)
	for (size_t i = 0; i < m_imp->m_Var.size(); ++i) delete m_imp->m_Var[i];
	m_imp->m_Var.clear();

	// clear the linear constraints
	if (m_imp->m_LCM) m_imp->m_LCM->Clear();

	// clear the mesh
	m_imp->m_mesh.Clear();

	// clear load parameters
	m_imp->m_Param.clear();
}

//-----------------------------------------------------------------------------
//! set the module name
void FEModel::SetActiveModule(const std::string& moduleName)
{
	m_imp->m_moduleName = moduleName;
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	fecore.SetActiveModule(moduleName.c_str());
	FEModule* pmod = fecore.GetActiveModule();
	pmod->InitModel(this);
}

//-----------------------------------------------------------------------------
//! get the module name
string FEModel::GetModuleName() const
{
	return m_imp->m_moduleName;
}

//-----------------------------------------------------------------------------
int FEModel::BoundaryConditions() const { return (int)m_imp->m_BC.size(); }

//-----------------------------------------------------------------------------
FEBoundaryCondition* FEModel::BoundaryCondition(int i) { return m_imp->m_BC[i]; }

//-----------------------------------------------------------------------------
void FEModel::AddBoundaryCondition(FEBoundaryCondition* pbc) { m_imp->m_BC.push_back(pbc); }

//-----------------------------------------------------------------------------
void FEModel::ClearBoundaryConditions() { m_imp->m_BC.clear(); }

//-----------------------------------------------------------------------------
int FEModel::InitialConditions() { return (int)m_imp->m_IC.size(); }

//-----------------------------------------------------------------------------
FEInitialCondition* FEModel::InitialCondition(int i) { return m_imp->m_IC[i]; }

//-----------------------------------------------------------------------------
void FEModel::AddInitialCondition(FEInitialCondition* pbc) { m_imp->m_IC.push_back(pbc); }

//-----------------------------------------------------------------------------
//! retrieve the number of steps
int FEModel::Steps() { return (int)m_imp->m_Step.size(); }

//-----------------------------------------------------------------------------
//! clear the steps
void FEModel::ClearSteps() { m_imp->m_Step.clear(); }

//-----------------------------------------------------------------------------
//! Add an analysis step
void FEModel::AddStep(FEAnalysis* pstep) { m_imp->m_Step.push_back(pstep); }

//-----------------------------------------------------------------------------
//! Get a particular step
FEAnalysis* FEModel::GetStep(int i) { return m_imp->m_Step[i]; }

//-----------------------------------------------------------------------------
//! Get the current step
FEAnalysis* FEModel::GetCurrentStep() { return m_imp->m_pStep; }
const FEAnalysis* FEModel::GetCurrentStep() const { return m_imp->m_pStep; }

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
int FEModel::SurfacePairConstraints() { return (int)m_imp->m_CI.size(); }

//-----------------------------------------------------------------------------
//! retrive a surface pair interaction
FESurfacePairConstraint* FEModel::SurfacePairConstraint(int i) { return m_imp->m_CI[i]; }

//-----------------------------------------------------------------------------
//! Add a surface pair interaction
void FEModel::AddSurfacePairConstraint(FESurfacePairConstraint* pci) { m_imp->m_CI.push_back(pci); }

//-----------------------------------------------------------------------------
//! return number of nonlinear constraints
int FEModel::NonlinearConstraints() { return (int)m_imp->m_NLC.size(); }

//-----------------------------------------------------------------------------
//! retrieve a nonlinear constraint
FENLConstraint* FEModel::NonlinearConstraint(int i) { return m_imp->m_NLC[i]; }

//-----------------------------------------------------------------------------
//! add a nonlinear constraint
void FEModel::AddNonlinearConstraint(FENLConstraint* pnlc) { m_imp->m_NLC.push_back(pnlc); }

//-----------------------------------------------------------------------------
//! return the number of model loads
int FEModel::ModelLoads() { return (int)m_imp->m_ML.size(); }

//-----------------------------------------------------------------------------
//! retrieve a model load
FEModelLoad* FEModel::ModelLoad(int i) { return m_imp->m_ML[i]; }

//-----------------------------------------------------------------------------
//! Add a model load
void FEModel::AddModelLoad(FEModelLoad* pml) { m_imp->m_ML.push_back(pml); }

//-----------------------------------------------------------------------------
// get the FE mesh
FEMesh& FEModel::GetMesh() { return m_imp->m_mesh; }

//-----------------------------------------------------------------------------
FELinearConstraintManager& FEModel::GetLinearConstraintManager() { return *m_imp->m_LCM; }

//-----------------------------------------------------------------------------
bool FEModel::Init()
{
	// make sure there is something to do
	if (m_imp->m_Step.size() == 0) return false;

	// intitialize time
	FETimeInfo& tp = GetTime();
	tp.currentTime = 0;
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
	// create and initialize the rigid body data
	// NOTE: Do this first, since some BC's look at the nodes' rigid id.
	if (InitRigidSystem() == false) return false;

	// evaluate all load controllers at the initial time
	for (int i = 0; i < LoadControllers(); ++i)
	{
		FELoadController* plc = m_imp->m_LC[i];
		if (plc->Init() == false)
		{
			std::string s = plc->GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Load controller %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
		plc->Evaluate(0);
	}

	// evaluate all mesh data generators
	for (int i = 0; i < MeshDataGenerators(); ++i)
	{
		FEMeshDataGenerator* pmd = m_imp->m_MD[i];
		if (pmd->Init() == false)
		{
			std::string s = pmd->GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Node data generator %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
		pmd->Evaluate(0);
	}

	// check step data
	for (int i = 0; i<(int)m_imp->m_Step.size(); ++i)
	{
		FEAnalysis& step = *m_imp->m_Step[i];
		if (step.Init() == false)
		{
			std::string s = step.GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Step %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
	}

	// validate BC's
	if (InitBCs() == false) return false;

	// initialize material data
	// NOTE: This must be called after the rigid system is initialiazed since the rigid materials will
	//       reference the rigid bodies
	if (InitMaterials() == false) return false;

	// initialize mesh data
	// NOTE: this must be done AFTER the elements have been assigned material point data !
	// this is because the mesh data is reset
	// TODO: perhaps I should not reset the mesh data during the initialization
	if (InitMesh() == false) return false;

	// initialize model loads
	// NOTE: This must be called after the InitMaterials since the COM of the rigid bodies
	//       are set in that function. 
	if (InitModelLoads() == false) return false;

	// initialize contact data
	if (InitContact() == false) return false;

	// initialize nonlinear constraints
	if (InitConstraints() == false) return false;

	// initialize mesh adaptors
	for (int i = 0; i < MeshAdaptors(); ++i)
	{
		FEMeshAdaptor* ma = MeshAdaptor(i);
		if (ma->Init() == false)
		{
			std::string s = ma->GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Mesh adaptor %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
	}

	// evaluate all load parameters
	// Do this last in case any model components redefined their load curves.
	if (EvaluateLoadParameters() == false) return false;

	// activate all permanent dofs
	Activate();

	// check if all load curves are being used
	int NLC = LoadControllers();
	vector<int> tag(NLC, 0);
	for (int i = 0; i < m_imp->m_Param.size(); ++i)
	{
		int lc = m_imp->m_Param[i].lc;
		tag[lc]++;
	}
	for (int i = 0; i < m_imp->m_Step.size(); ++i)
	{
		FEAnalysis* step = m_imp->m_Step[i];
		if (step->m_timeController)
		{
			int lc = step->m_timeController->m_nmplc;
			if (lc >= 0) tag[lc]++;
		}
	}
	int unused = 0;
	for (int i = 0; i < NLC; ++i) if (tag[i] == 0) unused++;
	if (unused != 0)
	{
		feLogWarning("Model has %d unreferenced load controllers.", unused);
	}

	bool ret = false;
	try
	{
		ret = DoCallback(CB_INIT);
	}
	catch (std::exception c)
	{
		ret = false;
		feLogError(c.what());
	}

	// do the callback
	return ret;
}

//-----------------------------------------------------------------------------
// get the number of calls to Update()
int FEModel::UpdateCounter() const
{
	return m_imp->m_nupdates;
}

//-----------------------------------------------------------------------------
void FEModel::IncrementUpdateCounter()
{
	m_imp->m_nupdates++;
}

//-----------------------------------------------------------------------------
void FEModel::Update()
{
	TRACK_TIME(TimerID::Timer_Update);

	// update model counter
	m_imp->m_nupdates++;
	
	// update mesh
	FEMesh& mesh = GetMesh();
	const FETimeInfo& tp = GetTime();
	mesh.Update(tp);

	// set the mesh update flag to false
	// If any load sets this to true, the
	// mesh will also be update after the loads are updated
	m_imp->m_meshUpdate = false;

	int nvel = BoundaryConditions();
	for (int i = 0; i < nvel; ++i)
	{
		FEBoundaryCondition& bc = *BoundaryCondition(i);
		if (bc.IsActive()) bc.UpdateModel();
	}

	// update all model loads
	for (int i = 0; i < ModelLoads(); ++i)
	{
		FEModelLoad* pml = ModelLoad(i);
		if (pml && pml->IsActive()) pml->Update();
	}

	// update all paired-interfaces
	for (int i = 0; i < SurfacePairConstraints(); ++i)
	{
		FESurfacePairConstraint* psc = SurfacePairConstraint(i);
		if (psc && psc->IsActive()) psc->Update();
	}

	// update all constraints
	for (int i = 0; i < NonlinearConstraints(); ++i)
	{
		FENLConstraint* pc = NonlinearConstraint(i);
		if (pc && pc->IsActive()) pc->Update();
	}

    // some of the loads may alter the prescribed dofs, so we update the mesh again
	if (m_imp->m_meshUpdate)
	{
		mesh.Update(tp);
		m_imp->m_meshUpdate = false;
	}
    
	// do the callback
	DoCallback(CB_MODEL_UPDATE);
}

//-----------------------------------------------------------------------------
//! See if the BC's are setup correctly.
bool FEModel::InitBCs()
{
	// check the IC's
	int NIC = InitialConditions();
	for (int i = 0; i<NIC; ++i)
	{
		FEInitialCondition* pic = InitialCondition(i);
		if (pic->Init() == false)
		{
			std::string s = pic->GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Initial condition %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
	}

	// check the BC's
	int NBC = BoundaryConditions();
	for (int i=0; i<NBC; ++i)
	{
		FEBoundaryCondition* pbc = BoundaryCondition(i);
		if (pbc->Init() == false)
		{
			std::string s = pbc->GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Boundary condition %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
	}

    return true;
}

//-----------------------------------------------------------------------------
void FEModel::AddMaterial(FEMaterial* pm) 
{ 
	m_imp->m_MAT.push_back(pm);
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
FEMaterial* FEModel::FindMaterial(const std::string& matName)
{
	for (int i = 0; i<Materials(); ++i)
	{
		FEMaterial* mat = GetMaterial(i);
		if (mat->GetName() == matName) return mat;
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
			feLogError("Failed initializing material %d (name=\"%s\")", i+1, pmat->GetName().c_str());
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
			feLogError("Failed validating material %d (name=\"%s\")", i+1, pmat->GetName().c_str());
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Add a loadcurve to the model
void FEModel::AddLoadController(FELoadController* plc) 
{ 
	m_imp->m_LC.push_back(plc) ;
}

//-----------------------------------------------------------------------------
void FEModel::ReplaceLoadController(int n, FELoadController* plc)
{
	assert((n >= 0) && (n < LoadControllers()));
	delete m_imp->m_LC[n];
	m_imp->m_LC[n] = plc;
}

//-----------------------------------------------------------------------------
//! get a loadcurve
FELoadController* FEModel::GetLoadController(int i)
{ 
	return m_imp->m_LC[i];
}

//-----------------------------------------------------------------------------
//! get the number of loadcurves
int FEModel::LoadControllers() const 
{ 
	return (int)m_imp->m_LC.size();
}

//-----------------------------------------------------------------------------
//! Attach a load controller to a parameter
void FEModel::AttachLoadController(FEParam* param, int lc)
{
	Implementation::LoadParam lp;
	lp.param = param;
	lp.lc = lc;

	switch (param->type())
	{
	case FE_PARAM_DOUBLE: lp.m_scl  = param->value<double>(); break;
	case FE_PARAM_VEC3D : lp.m_vscl = param->value<vec3d >(); break;
	}

	m_imp->m_Param.push_back(lp);
}

//-----------------------------------------------------------------------------
void FEModel::AttachLoadController(FEParam* p, FELoadController* plc)
{
	AttachLoadController(p, plc->GetID());
}

//-----------------------------------------------------------------------------
//! Detach a load controller from a parameter
bool FEModel::DetachLoadController(FEParam* p)
{
	for (int i = 0; i < (int)m_imp->m_Param.size(); ++i)
	{
		Implementation::LoadParam& pi = m_imp->m_Param[i];
		if (pi.param == p)
		{
			m_imp->m_Param.erase(m_imp->m_Param.begin() + i);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Get a load controller for a parameter (returns null if the param is not under load control)
FELoadController* FEModel::GetLoadController(FEParam* p)
{
	for (int i = 0; i < (int)m_imp->m_Param.size(); ++i)
	{
		Implementation::LoadParam& pi = m_imp->m_Param[i];
		if (pi.param == p)
		{
			return (pi.lc >= 0 ? GetLoadController(pi.lc) : nullptr);
		}
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
//! Add a mesh data generator to the model
void FEModel::AddMeshDataGenerator(FEMeshDataGenerator* pmd)
{
	m_imp->m_MD.push_back(pmd);
}

//-----------------------------------------------------------------------------
FEMeshDataGenerator* FEModel::GetMeshDataGenerator(int i)
{
	return m_imp->m_MD[i];
}

//-----------------------------------------------------------------------------
//! get the number of mesh data generators
int FEModel::MeshDataGenerators() const
{
	return (int)m_imp->m_MD.size();
}

//-----------------------------------------------------------------------------
//! Initialize rigid force data
bool FEModel::InitModelLoads()
{
	// call the Init() function of all rigid forces
	for (int i=0; i<ModelLoads(); ++i)
	{
		FEModelLoad& FC = *ModelLoad(i);
		if (FC.Init() == false)
		{
			std::string s = FC.GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Load %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Does one-time initialization of the Mesh data. Call FEMesh::Reset for resetting 
//! the mesh data.
bool FEModel::InitMesh()
{
	FEMesh& mesh = GetMesh();

	// find and remove isolated vertices
	int ni = mesh.RemoveIsolatedVertices();
	if (ni != 0)
	{
		if (ni == 1)
			feLogWarning("%d isolated vertex removed.", ni);
		else
			feLogWarning("%d isolated vertices removed.", ni);
	}

	// Initialize shell data
	// This has to be done before the domains are initialized below
	InitShells();

	// reset data
	// TODO: Not sure why this is here
	mesh.Reset();

	// initialize all domains
	// Initialize shell domains first (in order to establish SSI)
	// TODO: I'd like to move the initialization of the SSI to InitShells, but I can't 
	//       do that because FESSIShellDomain::FindSSI depends on the FEDomain::m_Node array which is
	//       initialized in FEDomain::Init.
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		if (dom.Class() == FE_DOMAIN_SHELL)
			if (dom.Init() == false) return false;
	}
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		if (dom.Class() != FE_DOMAIN_SHELL)
			if (dom.Init() == false) return false;
	}

	// initialize surfaces
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		if (mesh.Surface(i).Init() == false) return false;
	}

	// All done
	return true;
}

//-----------------------------------------------------------------------------
void FEModel::InitShells()
{
	FEMesh& mesh = GetMesh();

	// calculate initial directors for shell nodes
	int NN = mesh.Nodes();
	vector<vec3d> D(NN, vec3d(0, 0, 0));
	vector<int> ND(NN, 0);

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		// Calculate the shell directors as the local node normals
		if (mesh.Domain(nd).Class() == FE_DOMAIN_SHELL)
		{
			FEShellDomain& sd = static_cast<FEShellDomain&>(mesh.Domain(nd));
			vec3d r0[FEElement::MAX_NODES];
			for (int i = 0; i<sd.Elements(); ++i)
			{
				FEShellElement& el = sd.Element(i);

				int n = el.Nodes();
				int* en = &el.m_node[0];

				// get the nodes
				for (int j = 0; j<n; ++j) r0[j] = mesh.Node(en[j]).m_r0;
				for (int j = 0; j<n; ++j)
				{
					int m0 = j;
					int m1 = (j + 1) % n;
					int m2 = (j == 0 ? n - 1 : j - 1);

					vec3d a = r0[m0];
					vec3d b = r0[m1];
					vec3d c = r0[m2];
					vec3d d = (b - a) ^ (c - a); d.unit();

					D[en[m0]] += d*el.m_h0[j];
					++ND[en[m0]];
				}
			}
		}
	}

	// assign initial directors to shell nodes
	// make sure we average the directors
	for (int i = 0; i<NN; ++i)
		if (ND[i] > 0) mesh.Node(i).m_d0 = D[i] / ND[i];

	// do any other shell initialization 
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		if (dom.Class() == FE_DOMAIN_SHELL)
		{
			FEShellDomain& shellDom = static_cast<FEShellDomain&>(dom);
			shellDom.InitShells();
		}
	}
}

//-----------------------------------------------------------------------------
//! Initializes contact data
bool FEModel::InitContact()
{
	// loop over all contact interfaces
    for (int i=0; i<SurfacePairConstraints(); ++i)
	{
		// get the contact interface
        FESurfacePairConstraint& ci = *SurfacePairConstraint(i);

		// initializes contact interface data
		if (ci.Init() == false)
		{
			std::string s = ci.GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Contact %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
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
		if (plc->Init() == false)
		{
			std::string s = plc->GetName();
			const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
			feLogError("Nonlinear constraint %d (%s) failed to initialize", i + 1, sz);
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function solves the FE problem by calling the solve method for each step.
bool FEModel::Solve()
{
	TRACK_TIME(Timer_ModelSolve);

	// error flag
	bool bok = true;

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
			bok = false;
			break;
		}

		// do callback
		DoCallback(CB_STEP_ACTIVE);

		// solve the analaysis step
		bok = m_imp->m_pStep->Solve();

		if (nstep + 1 == Steps())
		{
			// set the solved flag
			m_imp->m_bsolved = bok;
		}

		// do callbacks
		DoCallback(CB_STEP_SOLVED);

		// wrap it up
		m_imp->m_pStep->Deactivate();

		// break if the step has failed
		if (bok == false) break;
	}

	// do the callbacks
	DoCallback(CB_SOLVED);

	return bok;
}

//-----------------------------------------------------------------------------
bool FEModel::RCI_Rewind()
{
	return m_imp->PopState();
}

//-----------------------------------------------------------------------------
bool FEModel::RCI_ClearRewindStack()
{
	if (m_imp->m_dmp.size() == 0) return false;
	m_imp->m_dmp.clear();
	return true;
}

//-----------------------------------------------------------------------------
bool FEModel::RCI_Init()
{
	// start the timer
	GetTimer(Timer_ModelSolve)->start();

	// reset solver status flag
	m_imp->m_bsolved = false;

	// loop over all analysis steps
	int nstep = m_imp->m_nStep;
	m_imp->m_pStep = m_imp->m_Step[(int)nstep];

	FEAnalysis* step = m_imp->m_pStep;

	// intitialize step data
	if (step->Activate() == false)
	{
		return false;
	}

	// do callback
	DoCallback(CB_STEP_ACTIVE);

	// initialize the step's solver
	if (step->InitSolver() == false)
	{
		return false;
	}

	return true;
}

bool FEModel::RCI_Restart()
{
	FEAnalysis* step = GetCurrentStep();
	if (step == nullptr) return false;

	return step->InitSolver();
}

bool FEModel::RCI_Advance()
{
	// get the current step
	FEAnalysis* step = m_imp->m_pStep;
	if (step == nullptr) return false;

	// first see if the step has finished
	const double eps = step->m_tend * 1e-7;
	double currentTime = GetCurrentTime();
	if (step->m_tend - currentTime <= eps)
	{
		// TODO: not sure why this is needed.
		SetStartTime(GetCurrentTime());

		// wrap it up
		DoCallback(CB_STEP_SOLVED);
		step->Deactivate();

		// go to the next step
		int nstep = ++m_imp->m_nStep;
		if (nstep >= m_imp->m_Step.size())
		{
			// we're done
			m_imp->m_bsolved = true;
			return true;
		}
		else
		{
			// go to the next step
			step = m_imp->m_pStep = m_imp->m_Step[nstep];
			if (step->Activate() == false) return false;
			DoCallback(CB_STEP_ACTIVE);
			if (step->InitSolver() == false) return false;
		}
	}

	// store current state in case we need to rewind
	m_imp->PushState();

	// Inform that the time is about to change. (Plugins can use 
	// this callback to modify time step)
	DoCallback(CB_UPDATE_TIME);

	// update time
	FETimeInfo& tp = GetTime();
	double newTime = tp.currentTime + step->m_dt;
	tp.currentTime = newTime;
	tp.timeIncrement = step->m_dt;
	feLog("\n===== beginning time step %d : %lg =====\n", step->m_ntimesteps + 1, newTime);

	// initialize the solver step
	// (This basically evaluates all the parameter lists, but let's the solver
	//  customize this process to the specific needs of the solver)
	if (step->GetFESolver()->InitStep(newTime) == false) return false;

	// Solve the time step
	int ierr = step->SolveTimeStep();
	if (ierr != 0) return false;

	// update counters
	FESolver* psolver = step->GetFESolver();
	step->m_ntotref += psolver->m_ntotref;
	step->m_ntotiter += psolver->m_niter;
	step->m_ntotrhs += psolver->m_nrhs;

	// Yes! We have converged!
	feLog("\n------- converged at time : %lg\n\n", GetCurrentTime());

	// update nr of completed timesteps
	step->m_ntimesteps++;

	// call callback function
	if (DoCallback(CB_MAJOR_ITERS) == false)
	{
		feLogWarning("Early termination on user's request");
		return false;
	}

	return true;
}

bool FEModel::RCI_Finish()
{
	// stop the timer
	GetTimer(Timer_ModelSolve)->stop();

	// do the callbacks
	DoCallback(CB_SOLVED);
	return true;
}

//-----------------------------------------------------------------------------
// Model activation.
// BC's that are not assigned to a step will not have their Activate() member called
// so we do it here. This function is called in Init() and Reset()
void FEModel::Activate()
{
	// initial conditions
	// Must be activated before prescribed BC's
	// since relative prescribed BC's use the initial values
	for (int i=0; i<InitialConditions(); ++i)
	{
		FEInitialCondition& ic = *InitialCondition(i);
		if (ic.IsActive()) ic.Activate();
	}

	// Boundary conditions
	for (int i = 0; i<BoundaryConditions(); ++i)
	{
		FEBoundaryCondition& bc = *BoundaryCondition(i);
		if (bc.IsActive()) bc.Activate();
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

	// initialize material points before evaluating contact autopenalty
    m_imp->m_mesh.InitMaterialPoints();
    
	// contact interfaces
    for (int i=0; i<SurfacePairConstraints(); ++i)
	{
        FESurfacePairConstraint& ci = *SurfacePairConstraint(i);
		if (ci.IsActive()) ci.Activate();
	}
}

//-----------------------------------------------------------------------------
// TODO: temporary construction. Need to see if I can just use Activate(). 
//       This is called after remeshed
void FEModel::Reactivate()
{
	// reactivate BCs
	for (int i = 0; i < BoundaryConditions(); ++i)
	{
		FEBoundaryCondition& bc = *BoundaryCondition(i);
		if (bc.IsActive()) bc.Activate();
	}

	// reactivate model loads
	for (int i = 0; i < ModelLoads(); ++i)
	{
		FEModelLoad& ml = *ModelLoad(i);
		if (ml.IsActive()) ml.Activate();
	}

	// update surface interactions
	for (int i = 0; i < SurfacePairConstraints(); ++i)
	{
		FESurfacePairConstraint& ci = *SurfacePairConstraint(i);
		if (ci.IsActive()) ci.Activate();
	}

	// reactivate the linear constraints
	GetLinearConstraintManager().Activate();
}

//-----------------------------------------------------------------------------
//! \todo Do I really need this function. I think calling FEModel::Init achieves the
//! same effect.
bool FEModel::Reset()
{
	// reset all timers
	ResetAllTimers();

	// reset solved flag
	m_imp->m_bsolved = false;

	// initialize materials
	for (int i=0; i<Materials(); ++i)
	{
		FEMaterial* pmat = GetMaterial(i);
		pmat->Init();
	}

	// reset mesh data
	m_imp->m_mesh.Reset();

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
	m_imp->m_timeInfo.currentTime = 0;
	m_imp->m_timeInfo.timeIncrement = 0;
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

	// Reevaluate load parameters
	EvaluateLoadParameters();

	DoCallback(CB_RESET);

	return true;
}

//-----------------------------------------------------------------------------
//! Get the current time information.
FETimeInfo& FEModel::GetTime()
{
	return m_imp->m_timeInfo;
}

//-----------------------------------------------------------------------------
double FEModel::GetStartTime() const { return m_imp->m_ftime0; }

//-----------------------------------------------------------------------------
void FEModel::SetStartTime(double t) { m_imp->m_ftime0 = t; }

//-----------------------------------------------------------------------------
double FEModel::GetCurrentTime() const { return m_imp->m_timeInfo.currentTime; }

//-----------------------------------------------------------------------------
void FEModel::SetCurrentTime(double t) { m_imp->m_timeInfo.currentTime = t; }

//-----------------------------------------------------------------------------
void FEModel::SetCurrentTimeStep(double dt)
{ 
	FEAnalysis* step = GetCurrentStep(); assert(step);
	if (step) step->m_dt = dt;
}

//=============================================================================
//    P A R A M E T E R   F U N C T I O N S
//=============================================================================

//-----------------------------------------------------------------------------
FEParam* FEModel::FindParameter(const ParamString& s)
{
	// make sure it starts with the name of this model
	if (s != GetName()) return 0;
	return FECoreBase::FindParameter(s.next());
}

//-----------------------------------------------------------------------------
FEParamValue GetComponent(FEParamValue& p, const ParamString& c)
{
	switch (p.type())
	{
	case FE_PARAM_MAT3DS:
		if (c == "xx") return p.component(0);
		if (c == "xy") return p.component(1);
		if (c == "yy") return p.component(2);
		if (c == "xz") return p.component(3);
		if (c == "yz") return p.component(4);
		if (c == "zz") return p.component(5);
		break;
	}
	return FEParamValue();
}

//-----------------------------------------------------------------------------
FEParamValue GetComponent(vec3d& r, const ParamString& c)
{
	if (c == "x") return FEParamValue(0, &r.x, FE_PARAM_DOUBLE);
	if (c == "y") return FEParamValue(0, &r.y, FE_PARAM_DOUBLE);
	if (c == "z") return FEParamValue(0, &r.z, FE_PARAM_DOUBLE);
	return FEParamValue();
}

//-----------------------------------------------------------------------------
// helper function for evaluating mesh data
FEParamValue FEModel::GetMeshParameter(const ParamString& paramString)
{
	FEMesh& mesh = GetMesh();

	ParamString next = paramString.next();
	if (next == "node")
	{
		FENode* node = 0;
		int nid = next.Index();
		if (nid < 0)
		{
			ParamString fnc = next.next();
			if (fnc == "fromId")
			{
				nid = fnc.Index();
				node = mesh.FindNodeFromID(nid);
				next = next.next();
			}
		}
		else if ((nid >= 0) && (nid < mesh.Nodes()))
		{
			node = &mesh.Node(nid);
		}

		if (node)
		{
			ParamString paramString = next.next();
			if      (paramString == "position"      ) return GetComponent(node->m_rt, paramString.next());
			// TODO: the m_Fr is not a vec3d any more, so not sure what to do here.
			//       In any case, this should probably be handled by the FEMechModel
//			else if (paramString == "reaction_force") return GetComponent(node->m_Fr, paramString.next());
			else
			{
				// see if it corresponds to a solution variable
				int n = GetDOFIndex(paramString.c_str());
				if ((n >= 0) && (n < node->dofs()))
				{
					return FEParamValue(0, &node->get(n), FE_PARAM_DOUBLE);
				}
			}
		}
	}
	else if (next == "elem")
	{
		FEElement* elem = 0;
		int nid = next.Index();
		if (nid < 0)
		{
			ParamString fnc = next.next();
			if (fnc == "fromId")
			{
				nid = fnc.Index();
				elem = mesh.FindElementFromID(nid);
				next = next.next();
			}
		}
		else if ((nid >= 0) && (nid < mesh.Elements()))
		{
			elem = mesh.Element(nid);
		}

		if (elem)
		{
			ParamString paramString = next.next();
			FEDomain* dom = dynamic_cast<FEDomain*>(elem->GetMeshPartition());
			if (dom == nullptr) return FEParamValue();

			if (paramString == "var")
			{
				std::string varName = paramString.IDString();
				FEMaterial* mat = dom->GetMaterial();
				if (mat)
				{
					FEDomainParameter* var = mat->FindDomainParameter(varName);
					if (var)
					{
						FEParamValue v = var->value(*elem->GetMaterialPoint(0));
						if (v.isValid())
						{
							ParamString comp = paramString.next();
							return GetComponent(v, comp);
						}
					}
				}
			}
		}
	}
	else if (next == "domain")
	{
		int nid = next.Index();
		if ((nid >= 0) && (nid < mesh.Domains()))
		{
			FEDomain& dom = mesh.Domain(nid);
			ParamString paramName = next.next();
			FEParam* param = dom.FindParameter(paramName);
			if (param)
			{
				if (param->type() == FE_PARAM_DOUBLE_MAPPED)
				{
					FEParamDouble& v = param->value<FEParamDouble>();
					if (v.isConst()) return FEParamValue(param, &v.constValue(), FE_PARAM_DOUBLE);
				}
				return FEParamValue(param, param->data_ptr(), param->type());
			}
		}
	}

	// if we get here, we did not find it
	return FEParamValue();
}

//-----------------------------------------------------------------------------
//! Return a pointer to the named variable
//! This function returns a pointer to a named variable.

FEParamValue FEModel::GetParameterValue(const ParamString& paramString)
{
	// make sure it starts with the name of this model
	if (paramString != GetName()) return FEParamValue();

	ParamString paramComp = paramString.last();
	FEParam* param = FindParameter(paramString);
	if (param)
	{
		if ((strcmp(param->name(), paramComp.c_str()) != 0) || (paramComp.Index() != -1))
			return GetParameterComponent(paramComp, param);
		else
		{
			if (param->type() == FE_PARAM_DOUBLE_MAPPED)
			{
				FEParamDouble& v = param->value<FEParamDouble>();
				if (v.isConst()) return FEParamValue(param, &v.constValue(), FE_PARAM_DOUBLE);
			}
			return FEParamValue(param, param->data_ptr(), param->type());
		}
	}

	// see what the next reference is
	ParamString next = paramString.next();

	// if we get here, handle some special cases
	if (next == "mesh") return GetMeshParameter(next);

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
void FEModel::EvaluateLoadControllers(double time)
{
	const int NLC = LoadControllers();
	for (int i=0; i<NLC; ++i) GetLoadController(i)->Evaluate(time);
}

//-----------------------------------------------------------------------------
//! Evaluates all load curves at the specified time
void FEModel::EvaluateDataGenerators(double time)
{
	for (int i = 0; i < MeshDataGenerators(); ++i) GetMeshDataGenerator(i)->Evaluate(time);
}

//-----------------------------------------------------------------------------
//! Set the print parameters flag
void FEModel::SetPrintParametersFlag(bool b)
{
	m_imp->m_printParams = b;
}

//-----------------------------------------------------------------------------
//! Get the print parameter flag
bool FEModel::GetPrintParametersFlag() const
{
	return m_imp->m_printParams;
}

//-----------------------------------------------------------------------------
bool FEModel::EvaluateLoadParameters()
{
	feLog("\n");
	int NLC = LoadControllers();
	for (int i = 0; i<(int)m_imp->m_Param.size(); ++i)
	{
		Implementation::LoadParam& pi = m_imp->m_Param[i];
		int nlc = pi.lc;
		if ((nlc >= 0) && (nlc < NLC))
		{
			double s = GetLoadController(nlc)->Value();
			FEParam* p = pi.param;
			FECoreBase* parent = dynamic_cast<FECoreBase*>(p->parent());
			if (m_imp->m_printParams)
			{
				if (parent && (parent->GetName().empty() == false))
				{
					const char* pname = parent->GetName().c_str();
					feLog("Setting parameter \"%s.%s\" to : ", pname, p->name());
				}
				else
					feLog("Setting parameter \"%s\" to : ", p->name());
			};
			assert(p->IsVolatile());
			switch (p->type())
			{
			case FE_PARAM_INT: {
				p->value<int>() = (int)s;
				if (m_imp->m_printParams) feLog("%d\n", p->value<int>());
			}
			break;
			case FE_PARAM_DOUBLE: {
				p->value<double>() = pi.m_scl*s;
				if (m_imp->m_printParams) feLog("%lg\n", p->value<double>());
			}
			break;
			case FE_PARAM_BOOL: {
				p->value<bool>() = (s > 0 ? true : false); 
				if (m_imp->m_printParams) feLog("%s\n", (p->value<bool>() ? "true" : "false"));
			}
			break;
			case FE_PARAM_VEC3D: {
				vec3d& v = p->value<vec3d>();
				p->value<vec3d>() = pi.m_vscl*s;
				if (m_imp->m_printParams) feLog("%lg, %lg, %lg\n", v.x, v.y, v.z);
			}
			break;
			case FE_PARAM_DOUBLE_MAPPED: 
			{
				FEParamDouble& v = p->value<FEParamDouble>();
				double c = 1.0;
				if (v.isConst()) c = v.constValue();
				v.SetScaleFactor(s * pi.m_scl);
				if (m_imp->m_printParams) feLog("%lg\n", c*p->value<FEParamDouble>().GetScaleFactor());
			}
			break;
			case FE_PARAM_VEC3D_MAPPED :
			{
				FEParamVec3& v = p->value<FEParamVec3>();
				v.SetScaleFactor(s * pi.m_scl);
				if (m_imp->m_printParams) feLog("%lg\n", v.GetScaleFactor());
			}
			break;
			default:
				feLog("\n");
				assert(false);
			}
		}
		else
		{
			feLogError("Invalid load curve ID");
			return false;
		}
	}
	feLog("\n");

	return true;
}

//-----------------------------------------------------------------------------
DOFS& FEModel::GetDOFS()
{
	return m_imp->m_dofs;
}

//-----------------------------------------------------------------------------
int FEModel::GetDOFIndex(const char* sz) const
{
	return m_imp->m_dofs.GetDOF(sz);
}

//-----------------------------------------------------------------------------
int FEModel::GetDOFIndex(const char* szvar, int n) const
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
	catch (ForceConversion)
	{
		throw;
	}
	catch (IterationFailure)
	{
		throw;
	}
	catch (DoRunningRestart)
	{
		throw;
	}
	catch (std::exception e)
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
void FEModel::Log(int ntag, const char* msg)
{

}

//-----------------------------------------------------------------------------
void FEModel::Logf(int ntag, const char* msg, ...)
{
	if (m_imp->m_block_log) return;

	// get a pointer to the argument list
	va_list	args;

	// make the message
	char sztxt[2048] = { 0 };
	va_start(args, msg);
	vsprintf(sztxt, msg, args);
	va_end(args);

	Log(ntag, sztxt);
}

//-----------------------------------------------------------------------------
void FEModel::BlockLog()
{
	m_imp->m_block_log = true;
}

//-----------------------------------------------------------------------------
void FEModel::UnBlockLog()
{
	m_imp->m_block_log = false;
}

//-----------------------------------------------------------------------------
bool FEModel::LogBlocked() const
{
	return m_imp->m_block_log;
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
int FEModel::GlobalVariables() const
{
	return (int)m_imp->m_Var.size();
}

//-----------------------------------------------------------------------------
void FEModel::AddGlobalVariable(const string& s, double v)
{
	FEGlobalVariable* var = new FEGlobalVariable;
	var->v = v;
	var->name = s;
	AddParameter(var->v, var->name.c_str());
	m_imp->m_Var.push_back(var);
}

const FEGlobalVariable& FEModel::GetGlobalVariable(int n)
{
	return *m_imp->m_Var[n];
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
FEGlobalData* FEModel::FindGlobalData(const char* szname)
{
	for (int i = 0; i < m_imp->m_GD.size(); ++i)
	{
		if (m_imp->m_GD[i]->GetName() == szname) return m_imp->m_GD[i];
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
int FEModel::FindGlobalDataIndex(const char* szname)
{
	for (int i = 0; i < m_imp->m_GD.size(); ++i)
	{
		if (m_imp->m_GD[i]->GetName() == szname) return i;
	}
	return -1;
}

//-----------------------------------------------------------------------------
int FEModel::GlobalDataItems()
{
	return (int)m_imp->m_GD.size();
}

//-----------------------------------------------------------------------------
FECoreBase* CopyFEBioClass(FECoreBase* pc, FEModel* fem)
{
	if ((pc == nullptr) || (fem == nullptr))
	{
		assert(false);
		return nullptr;
	}

	const char* sztype = pc->GetTypeStr();

	// create a new material
	FECoreBase* pcnew = fecore_new<FECoreBase>(pc->GetSuperClassID(), sztype, fem);
	assert(pcnew);

	pcnew->SetID(pc->GetID());

	// copy parameters
	pcnew->GetParameterList() = pc->GetParameterList();

	// copy properties
	for (int i = 0; i < pc->PropertyClasses(); ++i)
	{
		FEProperty* prop = pc->PropertyClass(i);
		if (prop->size() > 0)
		{
			for (int j = 0; j < prop->size(); ++j)
			{
				FECoreBase* pci = prop->get(j);
				if (pc)
				{
					FECoreBase* pci_new = CopyFEBioClass(pci, fem); assert(pci_new);
					bool b = pcnew->SetProperty(i, pci_new); assert(b);
				}
			}
		}
	}

	return pcnew;
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
	m_imp->m_nStep = fem.m_imp->m_nStep;
	m_imp->m_timeInfo = fem.m_imp->m_timeInfo;
	m_imp->m_ftime0 = fem.m_imp->m_ftime0;
	m_imp->m_pStep = 0;

	// copy model variables
	// we only copy the user created parameters, which presumably don't exist yet
	// in this model.
	int NS = fem.GlobalVariables();
	if (NS > 0)
	{
		assert(GlobalVariables() == 0);
		FEParameterList& PL = GetParameterList();
		for (int i = 0; i < NS; ++i)
		{
			const FEGlobalVariable& var = fem.GetGlobalVariable(i);
			AddGlobalVariable(var.name, var.v);
		}
	}

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
		pnew->CopyFrom(ps);

		// copy the solver
		FESolver* psolver = ps->GetFESolver();
		const char* sztype = psolver->GetTypeStr();

		// create a new solver
		FESolver* pnew_solver = fecore_new<FESolver>(sztype, this);
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

		// copy the material
		FEMaterial* pnew = dynamic_cast<FEMaterial*>(CopyFEBioClass(pmat, this));
		assert(pnew);

		// copy the name
		pnew->SetName(pmat->GetName());

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
		FEDomain* pd = nullptr;
		switch (dom.Class())
		{
		case FE_DOMAIN_SOLID   : pd = fecore_new<FESolidDomain   >(sz, this); break;
		case FE_DOMAIN_SHELL   : pd = fecore_new<FEShellDomain   >(sz, this); break;
		case FE_DOMAIN_BEAM    : pd = fecore_new<FEBeamDomain    >(sz, this); break;
		case FE_DOMAIN_2D      : pd = fecore_new<FEDomain2D      >(sz, this); break;
		case FE_DOMAIN_DISCRETE: pd = fecore_new<FEDiscreteDomain>(sz, this); break;
		}
		assert(pd);
		pd->SetMaterial(GetMaterial(LUT[i]));

		// copy domain data
		pd->CopyFrom(&dom);

		// add it to the mesh
		mesh.AddDomain(pd);
	}

	// --- boundary conditions ---

	int NDC = fem.BoundaryConditions();
	for (int i=0; i<NDC; ++i)
	{
		FEBoundaryCondition* pbc = fem.BoundaryCondition(i);
		const char* sz = pbc->GetTypeStr();

		FEBoundaryCondition* pnew = fecore_new<FEBoundaryCondition>(sz, this);
		assert(pnew);

		pnew->CopyFrom(pbc);

		// add to model
		AddBoundaryCondition(pnew);
	}

	// --- contact interfaces ---
    int NCI = fem.SurfacePairConstraints();
	for (int i=0; i<NCI; ++i)
	{
		// get the next interaction
        FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
		const char* sztype = pci->GetTypeStr();

		// create a new contact interface
        FESurfacePairConstraint* pnew = fecore_new<FESurfacePairConstraint>(sztype, this);
		assert(pnew);

		// create a copy
		pnew->CopyFrom(pci);

		// add the new interface
        AddSurfacePairConstraint(pnew);

		// add the surfaces to the surface list
		mesh.AddSurface(pnew->GetSecondarySurface());
		mesh.AddSurface(pnew->GetPrimarySurface ());
	}

	// --- nonlinear constraints ---
	int NLC = fem.NonlinearConstraints();
	for (int i=0; i<NLC; ++i)
	{
		// get the next constraint
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		const char* sztype = plc->GetTypeStr();

		// create a new nonlinear constraint
		FENLConstraint* plc_new = fecore_new<FENLConstraint>(sztype, this);
		assert(plc_new);

		// create a copy
		plc_new->CopyFrom(plc);

		// add the nonlinear constraint
		AddNonlinearConstraint(plc_new);

		// add the surface to the mesh (if any)
        FESurfaceConstraint* psc = dynamic_cast<FESurfaceConstraint*>(plc_new);
        if (psc)
        {
            FESurface* ps = psc->GetSurface();
            if (ps) mesh.AddSurface(ps);
        }
	}

	// --- Load curves ---
	// copy load curves
	int NLD = fem.LoadControllers();
	for (int i = 0; i<NLD; ++i)
	{
		FELoadController* lc = fem.GetLoadController(i);
		FELoadController* newlc = fecore_new<FELoadController>(lc->GetTypeStr(), this); assert(newlc);
		AddLoadController(newlc);
	}

	// copy linear constraints
	if (fem.m_imp->m_LCM)
	{
		if (m_imp->m_LCM) delete m_imp->m_LCM;
		m_imp->m_LCM = new FELinearConstraintManager(this);
		m_imp->m_LCM->CopyFrom(*fem.m_imp->m_LCM);
	}

	// copy output data
	m_imp->m_plotData = fem.m_imp->m_plotData;

	// TODO: copy all the properties
//	assert(false);
}

//-----------------------------------------------------------------------------
// This function serializes data to a stream.
// This is used for running and cold restarts.
void FEModel::Implementation::Serialize(DumpStream& ar)
{
	if (ar.IsShallow())
	{
		// stream model data
		ar & m_timeInfo;

		// stream mesh
		m_fem->SerializeGeometry(ar);

		// serialize contact
		ar & m_CI;

		// serialize nonlinear constraints
		ar & m_NLC;

		// serialize step and solver data
		ar & m_Step;
	}
	else
	{
		if (ar.IsLoading()) m_fem->Clear();

		ar & m_moduleName;

		if (ar.IsLoading())
		{
			FECoreKernel::GetInstance().SetActiveModule(m_moduleName.c_str());
		}

		ar & m_timeInfo;
		ar & m_dofs;
		ar & m_Const;
		ar & m_GD;
		ar & m_ftime0;
		ar & m_bsolved;

		// we have to stream materials before the mesh
		ar & m_MAT;

		// we have to stream the mesh before any boundary conditions
		m_fem->SerializeGeometry(ar);

		// stream all boundary conditions
		ar & m_BC;
		ar & m_ML;
		ar & m_IC;
		ar & m_CI;
		ar & m_NLC;

		// stream step data next
		ar & m_nStep;
		ar & m_Step;
		ar & m_pStep; // This must be streamed after m_Step

		// serialize linear constraints
		if (m_LCM) m_LCM->Serialize(ar);

		// serialize data generators
		ar & m_MD;

		// load controllers and load parameters are streamed last
		// since they can depend on other model parameters.
		ar & m_LC;
		ar & m_Param;
	}
}

//-----------------------------------------------------------------------------
//! This is called to serialize geometry.
//! Derived classes can override this
void FEModel::SerializeGeometry(DumpStream& ar)
{
	ar & m_imp->m_mesh;
}

//-----------------------------------------------------------------------------
// This function serializes data to a stream.
// This is used for running and cold restarts.
void FEModel::Serialize(DumpStream& ar)
{
	TRACK_TIME(TimerID::Timer_Update);

	m_imp->Serialize(ar);
	DoCallback(ar.IsSaving() ? CB_SERIALIZE_SAVE : CB_SERIALIZE_LOAD);
}

//-----------------------------------------------------------------------------
void FEModel::BuildMatrixProfile(FEGlobalMatrix& G, bool breset)
{
	FEAnalysis* pstep = GetCurrentStep();
	FESolver* solver = pstep->GetFESolver();
	solver->BuildMatrixProfile(G, breset);
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
// reset all the timers
void FEModel::ResetAllTimers()
{
	for (size_t i = 0; i<m_imp->m_timers.size(); ++i)
	{
		Timer& ti = m_imp->m_timers[i];
		ti.reset();
	}
}

//-----------------------------------------------------------------------------
int FEModel::Timers()
{
	return (int)m_imp->m_timers.size();
}

//-----------------------------------------------------------------------------
Timer* FEModel::GetTimer(int i)
{
	return &(m_imp->m_timers[i]);
}

//-----------------------------------------------------------------------------
//! return number of mesh adaptors
int FEModel::MeshAdaptors()
{
	return (int)m_imp->m_MA.size();
}

//-----------------------------------------------------------------------------
//! retrieve a mesh adaptors
FEMeshAdaptor* FEModel::MeshAdaptor(int i)
{
	return m_imp->m_MA[i];
}

//-----------------------------------------------------------------------------
//! add a mesh adaptor
void FEModel::AddMeshAdaptor(FEMeshAdaptor* meshAdaptor)
{
	m_imp->m_MA.push_back(meshAdaptor);
}

//-----------------------------------------------------------------------------
void FEModel::SetUnits(const char* szunits)
{
	if (szunits)
		m_imp->m_units = szunits;
	else
		m_imp->m_units.clear();
}

//-----------------------------------------------------------------------------
const char* FEModel::GetUnits() const
{
	if (m_imp->m_units.empty()) return nullptr;
	else return m_imp->m_units.c_str();
}
