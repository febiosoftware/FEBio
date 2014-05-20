#pragma once
#include "FEMaterial.h"
#include "FEMesh.h"
#include "LoadCurve.h"
#include "BC.h"
#include "FEBodyLoad.h"
#include "FESurfacePairInteraction.h"
#include "FEAnalysis.h"
#include "FESurfaceLoad.h"
#include "FENLConstraint.h"
#include "FELinearConstraint.h"
#include "FEObject.h"
#include "DumpStream.h"
#include "DataStore.h"
#include "FEGlobalData.h"
#include "FETypes.h"
#include <string>
#include <vector>
#include <map>

//-----------------------------------------------------------------------------
// Forward declaration of the FEModel class.
class FEModel;
class DataRecord;
class PlotFile;

//-----------------------------------------------------------------------------
// callback structure
#define CB_ALWAYS		0x00000011
#define CB_MAJOR_ITERS	0x00000001
#define CB_MINOR_ITERS	0x00000010

typedef unsigned int FECORE_CB_WHEN;
typedef void (*FECORE_CB_FNC)(FEModel*,void*);
struct FECORE_CALLBACK {
	FECORE_CB_FNC	m_pcb;		// pointer to callback function
	void*			m_pd;		// pointer to user data
	FECORE_CB_WHEN	m_nwhen;	// when to call function
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

public:
	// get the FE mesh
	FEMesh& GetMesh() { return m_mesh; }

	// get the rigid object
	FEObject* Object(int i) { return m_Obj[i]; }

	// return number of rigid objects
	int Objects() { return (int) m_Obj.size(); }

	//! Initialize mesh data
	bool InitMesh();

	//! Initialize the objects
	bool InitObjects();

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

	//! material initialization
	bool InitMaterials();

public: // --- Boundary Conditions functions ---

	// prescribed BC's
	int PrescribedBCs() { return (int) m_DC.size(); }
	FEPrescribedBC* PrescribedBC(int i) { return m_DC[i]; }
	void AddPrescribedBC(FEPrescribedBC* pbc) { m_DC.push_back(pbc); }
	void ClearBCs();

	// nodal loads
	int NodalLoads() { return (int) m_FC.size(); }
	FENodalForce* NodalLoad(int i) { return m_FC[i]; }
	void AddNodalLoad(FENodalForce* pfc) { m_FC.push_back(pfc); }

	// surface loads
	int SurfaceLoads() { return (int) m_SL.size(); }
	FESurfaceLoad* SurfaceLoad(int i) { return m_SL[i]; }
	void AddSurfaceLoad(FESurfaceLoad* psl) { m_SL.push_back(psl); }

	// rigid nodes
	int RigidNodes() { return (int) m_RN.size(); }
	FERigidNode* RigidNode(int i) { return m_RN[i]; }
	void AddRigidNode(FERigidNode* prn) { m_RN.push_back(prn); }

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
	FETimePoint GetCurrentTime();

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

public: // --- parameter functions ---

	//! evaluate all parameter lists
	void EvaluateAllParameterLists();

	//! Evaluate parameter list
	void EvaluateParameterList(FEParameterList& pl);

	//! Evaluate material parameters
	void EvaluateMaterialParameters(FEMaterial* pm);

	//! return a pointer to the named variable
	double* FindParameter(const char* szname);

public:	// --- Miscellaneous routines ---

	//! read/write the model state to a dump stream (used for running restarts)
	virtual void ShallowCopy(DumpStream& dmp, bool bsave);

	//! find a boundary condition from the ID
	FEBoundaryCondition* FindBC(int nid);

	//! set callback function
	void AddCallback(FECORE_CB_FNC pcb, unsigned int nwhen, void* pd);

	//! call the callback function
	void DoCallback(unsigned int nevent);

public: // --- I/O functions

	//! write to plot file
	virtual void Write() {}

	//! write data to log file
	virtual void WriteData() {}

	//! write data to dump file
	virtual void DumpData() {}

	//! serialize data
	virtual bool Serialize(DumpFile& ar) { return true; }

	//! Add data record
	virtual void AddDataRecord(DataRecord* pd) {}

public: // Global data
	void AddGlobalData(FEGlobalData* psd);
	FEGlobalData* GetGlobalData(int i);
	int GlobalDataItems();

	// get/set global data
	void SetGlobalConstant(const string& s, double v);
	double GetGlobalConstant(const string& s);

public:
	//! set the debug level
	void SetDebugFlag(bool b) { m_debug = b; }

	//! get the debug level
	bool GetDebugFlag() { return m_debug; }

public: // TODO: Find a better place for these parameters
	int		m_nsolver;			//!< type of solver selected
	int		m_bwopt;			//!< bandwidth optimization flag
	int		m_nStep;			//!< current analysis step
	double	m_ftime;			//!< current time value
	double	m_ftime0;			//!< start time of current step
	int		m_nplane_strain;	//!< run analysis in plain strain mode \todo Move to the analysis class?
	bool	m_debug;			//!< debug flag
	double	m_ut4_alpha;		//!< UT4 integration alpha value
	bool	m_ut4_bdev;			//!< UT4 integration deviatoric formulation flag
	double	m_udghex_hg;		//!< hourglass parameter for UDGhex integration

protected:
	std::vector<FELoadCurve*>				m_LC;	//!< load curve data
	std::vector<FEMaterial*>				m_MAT;	//!< array of materials
	std::vector<FEBodyLoad*>				m_BL;	//!< body load data
	std::vector<FESurfacePairInteraction*>	m_CI;	//!< contact interface array
	std::vector<FENLConstraint*>			m_NLC;	//!< nonlinear constraints

protected:
	std::vector<FEAnalysis*>	m_Step;		//!< array of analysis steps
	FEAnalysis*					m_pStep;	//!< pointer to current analysis step

protected:
	// Geometry data
	FEMesh		m_mesh;					//!< the one and only FE mesh
	std::vector<FEObject*>		m_Obj;	//!< FE Object array (NOTE: only used for rigid bodies)

protected:
	// Boundary Conditions
	std::vector<FEPrescribedBC*>	m_DC;	//!< prescribed constraints
	std::vector<FENodalForce*>		m_FC;	//!< concentrated nodal loads
	std::vector<FESurfaceLoad*>		m_SL;	//!< surface loads
	std::vector<FERigidNode*>		m_RN;	//!< rigid nodes

public:
	// Boundary conditions for rigid bodies
	// TODO: I'd like to do something different with this. Perhaps place them in the BC or in some constraint section.
	vector<FERigidBodyFixedBC*>			m_RBC;	//!< rigid body fixed
	vector<FERigidBodyDisplacement*>	m_RDC;	//!< rigid body displacements
	vector<FERigidBodyForce*>			m_RFC;	//!< rigid body forces

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
};
