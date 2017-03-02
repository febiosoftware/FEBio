#pragma once
#include "DOFS.h"
#include "FEMesh.h"
#include "FETypes.h"
#include "FEModelComponent.h"
#include "Callback.h"
#include <string>

//-----------------------------------------------------------------------------
// forward declarations
class FELoadCurve;
class FEMaterial;
class FEModelLoad;
class FENodalLoad;
class FEFixedBC;
class FEPrescribedBC;
class FEInitialCondition;
class FESurfaceLoad;
class FEEdgeLoad;
class FEBodyLoad;
class FENLConstraint;
class FESurfacePairInteraction;
class FERigidSystem;
class FEAnalysis;
class FEGlobalData;
class FEGlobalMatrix;
class FELinearConstraintManager;

//-----------------------------------------------------------------------------
//! The FEModel class stores all the data for the finite element model, including
//! geometry, analysis steps, boundary and loading conditions, contact interfaces
//! and so on.
//!
class FEModel : public FECoreBase, public CallbackHandler
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
	FEMesh& GetMesh();

	// get the rigid system
	FERigidSystem* GetRigidSystem();

	// get the linear constraint manager
	FELinearConstraintManager& GetLinearConstraintManager();

	//! Validate BC's
	bool InitBCs();

	//! Build the matrix profile for this model
	void BuildMatrixProfile(FEGlobalMatrix& G, bool breset);

public:
	//! set the problem title
	void SetTitle(const char* sz);

	//! get the problem title
	const char* GetTitle();

public:	// --- Load curve functions ----

	//! Add a loadcurve to the model
	void AddLoadCurve(FELoadCurve* plc);

	//! get a loadcurve
	FELoadCurve* GetLoadCurve(int i);

	//! get the number of loadcurves
	int LoadCurves() const;

public: // --- Material functions ---

	//! Add a material to the model
	void AddMaterial(FEMaterial* pm);

	//! get the number of materials
	int Materials();

	//! return a pointer to a material
	FEMaterial* GetMaterial(int i);

	//! find a material based on its index
	FEMaterial* FindMaterial(int nid);

	//! find a material based on its name
	FEMaterial* FindMaterial(const char* sz);

	//! material initialization
	bool InitMaterials();

	//! material validation
	bool ValidateMaterials();

public: // --- Boundary Conditions functions ---
	// fixed BC
	int FixedBCs();
	FEFixedBC* FixedBC(int i);
	void AddFixedBC(FEFixedBC* pbc);
	void AddFixedBC(int node, int bc);

	// prescribed BC's
	int PrescribedBCs();
	FEPrescribedBC* PrescribedBC(int i);
	void AddPrescribedBC(FEPrescribedBC* pbc);
	void ClearBCs();

	// initial conditions
	int InitialConditions();
	FEInitialCondition* InitialCondition(int i);
	void AddInitialCondition(FEInitialCondition* pbc);

	// nodal loads
	int NodalLoads();
	FENodalLoad* NodalLoad(int i);
	void AddNodalLoad(FENodalLoad* pfc);

	// surface loads
	int SurfaceLoads();
	FESurfaceLoad* SurfaceLoad(int i);
	void AddSurfaceLoad(FESurfaceLoad* psl);

	// edge loads
	int EdgeLoads();
	FEEdgeLoad* EdgeLoad(int i);
	void AddEdgeLoad(FEEdgeLoad* psl);

public: // --- Body load functions --- 

	//! Add a body load to the model
	void AddBodyLoad(FEBodyLoad* pf);

	//! get the number of body loads
	int BodyLoads();

	//! return a pointer to a body load
	FEBodyLoad* GetBodyLoad(int i);

	//! Init body loads
	bool InitBodyLoads();

public: // --- Analysis steps functions ---

	//! retrieve the number of steps
	int Steps();

	//! clear the steps
	void ClearSteps();

	//! Add an analysis step
	void AddStep(FEAnalysis* pstep);

	//! Get a particular step
	FEAnalysis* GetStep(int i);

	//! Get the current step
	FEAnalysis* GetCurrentStep();

	//! Set the current step index
	int GetCurrentStepIndex() const;

	//! Set the current step
	void SetCurrentStep(FEAnalysis* pstep);

	//! Set the current step index
	void SetCurrentStepIndex(int n);

	//! Get the current time
	FETimeInfo GetTime();

	//! Get the start time
	double GetStartTime() const;

	//! Set the start time
	void SetStartTime(double t);

	//! Get the current time
	double GetCurrentTime() const;

	//! Set the current time
	void SetCurrentTime(double t);

public: // --- Contact interface functions ---

	//! return number of surface pair interactions
	int SurfacePairInteractions();

	//! retrive a surface pair interaction
	FESurfacePairInteraction* SurfacePairInteraction(int i);

	//! Add a surface pair interaction
	void AddSurfacePairInteraction(FESurfacePairInteraction* pci);

	//! Initializes contact data
	bool InitContact();

public: // --- Nonlinear constraints functions ---

	//! return number of nonlinear constraints
	int NonlinearConstraints();

	//! retrieve a nonlinear constraint
	FENLConstraint* NonlinearConstraint(int i);

	//! add a nonlinear constraint
	void AddNonlinearConstraint(FENLConstraint* pnlc);

	//! Initialize constraint data
	bool InitConstraints();

public:	// --- Model Loads ----
	//! return the number of model loads
	int ModelLoads();

	//! retrieve a model load
	FEModelLoad* ModelLoad(int i);

	//! Add a model load
	void AddModelLoad(FEModelLoad* pml);

	//! initialize model loads
	bool InitModelLoads();

	//! find a surface load based on the name
	FESurfaceLoad* FindSurfaceLoad(const char* sz);

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
	double* FindParameter(const ParamString& param);

public:	// --- Miscellaneous routines ---

	//! find a model componnet from its class ID
	FEModelComponent* FindModelComponent(int nid);

	//! call the callback function
	//! This function returns fals if the run is to be aborted
	bool DoCallback(unsigned int nevent);

	//! I'd like to place the list of DOFS inside the model.
	//! As a first step, all classes that have access to the model
	//! should get the DOFS from this function
	DOFS& GetDOFS();

	//! Get the index of a DOF
	int GetDOFIndex(const char* sz);
	int GetDOFIndex(const char* szvar, int n);

	//! serialize data for restarts
	void Serialize(DumpStream& ar);

	//! Get the linear solver type
	int GetLinearSolverType() const;

	//! set the linear solver type
	void SetLinearSolverType(int ntype);

	//! see if we need to optimize bandwidth of linear system
	bool OptimizeBandwidth() const;

	//! Set the optimize band width flag
	void SetOptimizeBandwidth(bool b);

protected:
	// helper functions for serialization
	void SerializeLoadData    (DumpStream& ar);
	void SerializeGlobals     (DumpStream& ar);
	void SerializeMaterials   (DumpStream& ar);
	void SerializeGeometry    (DumpStream& ar);
	void SerializeContactData (DumpStream& ar);
	void SerializeBoundaryData(DumpStream& ar);
	void SerializeAnalysisData(DumpStream& ar);

public: // Global data
	void AddGlobalData(FEGlobalData* psd);
	FEGlobalData* GetGlobalData(int i);
	int GlobalDataItems();

	// get/set global data
	void SetGlobalConstant(const string& s, double v);
	double GetGlobalConstant(const string& s);

public:
	// TODO: put this somewhere else
	double	m_ut4_alpha;		//!< UT4 integration alpha value
	bool	m_ut4_bdev;			//!< UT4 integration deviatoric formulation flag
	double	m_udghex_hg;		//!< hourglass parameter for UDGhex integration

private:
	class Implementation;
	Implementation*	m_imp;
};
