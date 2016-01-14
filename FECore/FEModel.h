#pragma once
#include "DOFS.h"
#include "FERigidSystem.h"
#include "FEMaterial.h"
#include "FEMesh.h"
#include "LoadCurve.h"
#include "BC.h"
#include "FEInitialCondition.h"
#include "FEModelLoad.h"
#include "FEBodyLoad.h"
#include "FESurfacePairInteraction.h"
#include "FEAnalysis.h"
#include "FESurfaceLoad.h"
#include "FENLConstraint.h"
#include "FELinearConstraint.h"
#include "FEObject.h"
#include "DumpStream.h"
#include "FEGlobalData.h"
#include "FETypes.h"
#include <string>
#include <vector>
#include <map>

//-----------------------------------------------------------------------------
// callback structure
#define CB_ALWAYS		0x00011111		//!< Call for all reasons
#define CB_INIT			0x00000001		//!< Call after model initialization (i.e. FEModel::Init())
#define CB_MAJOR_ITERS	0x00000010		//!< Call at the end of each major converged iteration
#define CB_MINOR_ITERS	0x00000100		//!< Call for each minor iteration
#define CB_SOLVED		0x00001000		//!< Call at the end of FEModel::Solve
#define CB_UPDATE_TIME	0x00010000		//!< Call when time is updated and right before time step is solved (in FEAnalysis::Solve)

typedef unsigned int FECORE_CB_WHEN;
typedef void (*FECORE_CB_FNC)(FEModel*,unsigned int,void*);
struct FECORE_CALLBACK {
	FECORE_CB_FNC	m_pcb;		// pointer to callback function
	void*			m_pd;		// pointer to user data
	FECORE_CB_WHEN	m_nwhen;	// when to call function
};

//-----------------------------------------------------------------------------
//! These output flags are used when calling FEModel::Write to report the state
//! of the solution. This can be used to decide whether to output the state data or not.
enum FE_OUTPUT_HINT {
	FE_UNKNOWN,				// the model is in an unknown state
	FE_INITIALIZED,			// the model has been initialized (i.e. FEModel::Init() was called)
	FE_UNCONVERGED,			// the model is in a non-converged state
	FE_CONVERGED,			// the time step has converged
	FE_AUGMENT,				// the model is entering an augmented (called before Augment())
	FE_SOLVED,				// the entire model is solved
	FE_STEP_INITIALIZED,	// the current step is initialized (i.e. FEAnalysis::Init() was called)
	FE_STEP_SOLVED,			// the step was solved
};

//-----------------------------------------------------------------------------
//! The FEModel class stores all the data for the finite element model, including
//! geometry, analysis steps, boundary and loading conditions, contact interfaces
//! and so on.
//!
class FEModel
{
public:
	enum {MAX_STRING = 256};

public:
	FEModel(void);
	virtual ~FEModel(void);

	// Initialization
	virtual bool Init();

	//! Resets data structures
	virtual bool Reset();

	// solve the model
	virtual bool Solve();

	// copy the model data
	virtual void CopyFrom(FEModel& fem);

	// clear all model data
	void Clear();

	// model activation
	void Activate();

public:
	// get the FE mesh
	FEMesh& GetMesh() { return m_mesh; }

	// get the rigid system
	FERigidSystem* GetRigidSystem() { return m_prs; }

	//! Initialize mesh data
	bool InitMesh();

	//! Validate BC's
	bool InitBCs();

public:
	//! set the problem title
	void SetTitle(const char* sz);

	//! get the problem title
	const char* GetTitle();

public:	// --- Load curve functions ----

	//! Add a loadcurve to the model
	void AddLoadCurve(FELoadCurve* plc) { m_LC.push_back(plc); }

	//! get a loadcurve
	FELoadCurve* GetLoadCurve(int i) { return m_LC[i]; }

	//! get the number of loadcurves
	int LoadCurves() { return m_LC.size(); }

public: // --- Material functions ---

	//! Add a material to the model
	void AddMaterial(FEMaterial* pm) { m_MAT.push_back(pm); }

	//! get the number of materials
	int Materials() { return m_MAT.size(); }

	//! return a pointer to a material
	FEMaterial* GetMaterial(int i) { return m_MAT[i]; }

	//! find a material based on its index
	FEMaterial* FindMaterial(int nid);

	//! material initialization
	bool InitMaterials();

	//! material validation
	bool ValidateMaterials();

public: // --- Boundary Conditions functions ---
	// fixed BC
	int FixedBCs() { return (int) m_BC.size(); }
	FEFixedBC* FixedBC(int i) { return m_BC[i]; }
	void AddFixedBC(FEFixedBC* pbc) { m_BC.push_back(pbc); }
	void AddFixedBC(int node, int bc);

	// prescribed BC's
	int PrescribedBCs() { return (int) m_DC.size(); }
	FEPrescribedBC* PrescribedBC(int i) { return m_DC[i]; }
	void AddPrescribedBC(FEPrescribedBC* pbc) { m_DC.push_back(pbc); }
	void ClearBCs();

	// initial conditions
	int InitialConditions() { return (int) m_IC.size(); }
	FEInitialCondition* InitialCondition(int i) { return m_IC[i]; }
	void AddInitialCondition(FEInitialCondition* pbc) { m_IC.push_back(pbc); }

	// nodal loads
	int NodalLoads() { return (int) m_FC.size(); }
	FENodalLoad* NodalLoad(int i) { return m_FC[i]; }
	void AddNodalLoad(FENodalLoad* pfc) { m_FC.push_back(pfc); }

	// surface loads
	int SurfaceLoads() { return (int) m_SL.size(); }
	FESurfaceLoad* SurfaceLoad(int i) { return m_SL[i]; }
	void AddSurfaceLoad(FESurfaceLoad* psl) { m_SL.push_back(psl); }

public: // --- Body load functions --- 

	//! Add a body load to the model
	void AddBodyLoad(FEBodyLoad* pf) { m_BL.push_back(pf); }

	//! get the number of body loads
	int BodyLoads() { return (int) m_BL.size(); }

	//! return a pointer to a body load
	FEBodyLoad* GetBodyLoad(int i) { return m_BL[i]; }

	//! see if there are any body loads
	bool HasBodyLoads() { return !m_BL.empty();}

	//! Init body loads
	bool InitBodyLoads();

public: // --- Analysis steps functions ---

	//! retrieve the number of steps
	int Steps() { return (int) m_Step.size(); }

	//! clear the steps
	void ClearSteps() { m_Step.clear(); }

	//! Add an analysis step
	void AddStep(FEAnalysis* pstep) { m_Step.push_back(pstep); }

	//! Get a particular step
	FEAnalysis* GetStep(int i) { return m_Step[i]; }

	//! Get the current step
	FEAnalysis* GetCurrentStep() { return m_pStep; }

	//! Set the current step
	void SetCurrentStep(FEAnalysis* pstep) { m_pStep = pstep; }

	//! Get the current time
	FETimePoint GetTime();

public: // --- Contact interface functions ---

	//! return number of surface pair interactions
	int SurfacePairInteractions() { return m_CI.size(); } 

	//! retrive a surface pair interaction
	FESurfacePairInteraction* SurfacePairInteraction(int i) { return m_CI[i]; }

	//! Add a surface pair interaction
	void AddSurfacePairInteraction(FESurfacePairInteraction* pci) { m_CI.push_back(pci); }

	//! Initializes contact data
	bool InitContact();

public: // --- Nonlinear constraints functions ---

	//! return number of nonlinear constraints
	int NonlinearConstraints() { return (int) m_NLC.size(); }

	//! retrieve a nonlinear constraint
	FENLConstraint* NonlinearConstraint(int i) { return m_NLC[i]; }

	//! add a nonlinear constraint
	void AddNonlinearConstraint(FENLConstraint* pnlc) { m_NLC.push_back(pnlc); }

	//! Initialize constraint data
	bool InitConstraints();

public:	// --- Model Loads ----
	//! return the number of model loads
	int ModelLoads() { return (int) m_ML.size(); }

	//! retrieve a model load
	FEModelLoad* ModelLoad(int i) { return m_ML[i]; }

	//! Add a model load
	void AddModelLoad(FEModelLoad* pml) { m_ML.push_back(pml); }

	//! initialize model loads
	bool InitModelLoads();

public: // --- parameter functions ---

	//! evaluate all load curves at some time
	void EvaluateLoadCurves(double time);

	//! evaluate all parameter lists
	bool EvaluateAllParameterLists();

	//! Evaluate parameter list
	bool EvaluateParameterList(FEParameterList& pl);

	//! Evaluate parameter list
	bool EvaluateParameterList(FECoreBase* pc);

	//! return a pointer to the named variable
	double* FindParameter(const char* szname);

public:	// --- Miscellaneous routines ---

	//! read/write the model state to a dump stream (used for running restarts)
	virtual void ShallowCopy(DumpStream& dmp, bool bsave);

	//! find a model componnet from its class ID
	FEModelComponent* FindModelComponent(int nid);

	//! set callback function
	void AddCallback(FECORE_CB_FNC pcb, unsigned int nwhen, void* pd);

	//! call the callback function
	void DoCallback(unsigned int nevent);

	//! I'd like to place the list of DOFS inside the model.
	//! As a first step, all classes that have access to the model
	//! should get the DOFS from this function
	DOFS& GetDOFS();

	//! Get the index of a DOF
	int GetDOFIndex(const char* sz);
	int GetDOFIndex(const char* szvar, int n);

public: // --- I/O functions

	//! write to plot file
	virtual void Write(FE_OUTPUT_HINT hint) {}

	//! serialize data
	virtual bool Serialize(DumpFile& ar) { return true; }

public: // Global data
	void AddGlobalData(FEGlobalData* psd);
	FEGlobalData* GetGlobalData(int i);
	int GlobalDataItems();

	// get/set global data
	void SetGlobalConstant(const string& s, double v);
	double GetGlobalConstant(const string& s);

public: // TODO: Find a better place for these parameters
	// I want to make this parameter part of the FEAnalysis, since 
	// it could be different analysis steps (in multi-step problems) may
	// require different solvers.
	int		m_nsolver;			//!< type of (linear) solver selected

	int		m_bwopt;			//!< bandwidth optimization flag
	int		m_nStep;			//!< current analysis step
	double	m_ftime;			//!< current time value
	double	m_ftime0;			//!< start time of current step
	int		m_nplane_strain;	//!< run analysis in plain strain mode \todo Move to the analysis class?
	double	m_ut4_alpha;		//!< UT4 integration alpha value
	bool	m_ut4_bdev;			//!< UT4 integration deviatoric formulation flag
	double	m_udghex_hg;		//!< hourglass parameter for UDGhex integration

protected:
	std::vector<FELoadCurve*>				m_LC;	//!< load curve data
	std::vector<FEMaterial*>				m_MAT;	//!< array of materials
	std::vector<FEFixedBC*>					m_BC;	//!< fixed constraints
	std::vector<FEPrescribedBC*>			m_DC;	//!< prescribed constraints
	std::vector<FEInitialCondition*>		m_IC;	//!< initial conditions
	std::vector<FENodalLoad*>				m_FC;	//!< concentrated nodal loads
	std::vector<FESurfaceLoad*>				m_SL;	//!< surface loads
	std::vector<FESurfacePairInteraction*>	m_CI;	//!< contact interface array
	std::vector<FEBodyLoad*>				m_BL;	//!< body load data
	std::vector<FENLConstraint*>			m_NLC;	//!< nonlinear constraints
	std::vector<FEModelLoad*>				m_ML;	//!< model loads

protected:
	std::vector<FEAnalysis*>	m_Step;		//!< array of analysis steps
	FEAnalysis*					m_pStep;	//!< pointer to current analysis step

protected:
	// DOFS data
	DOFS	m_dofs;				//!< list of degree of freedoms in this model

	// Geometry data
	FEMesh		m_mesh;			//!< the one and only FE mesh

	// the rigid body system
	FERigidSystem*		m_prs;	//!< the rigid body system manages rigid bodies

public:
	// linear constraint data
	list<FELinearConstraint>	m_LinC;		//!< linear constraints data
	vector<int>					m_LCT;		//!< linear constraint table
	vector<FELinearConstraint*>	m_LCA;		//!< linear constraint array (temporary solution!)

protected:
	char	m_sztitle[MAX_STRING];	//!< problem title
	list<FECORE_CALLBACK>	m_pcb;	//!< pointer to callback function

protected: // Global Data
	std::map<string, double> m_Const;	//!< Global model constants
	vector<FEGlobalData*>	m_GD;		//!< global data structures

public:
	static void SetDefaultSolver(int nsolver) { m_ndefault_solver = nsolver; } 
	static int m_ndefault_solver;
};
