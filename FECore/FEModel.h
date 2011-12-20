#pragma once
#include "FEMaterial.h"
#include "FEMesh.h"
#include "LoadCurve.h"
#include "BC.h"
#include "FEBodyForce.h"
#include "FEElasticMaterial.h"		// TODO: I want to delete this
#include "FEAnalysis.h"
#include <vector>
#include <map>

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// FEBIO callback structure
typedef void (*FEBIO_CB_FNC)(FEModel*,void*);
struct FEBIO_CALLBACK {
	FEBIO_CB_FNC	m_pcb;
	void*	m_pd;
};

//-----------------------------------------------------------------------------
class FEModel
{
public:
	FEModel(void);
	virtual ~FEModel(void);

	// input from file
	virtual bool Input(const char* szfile) = 0;

	// Initialization
	virtual bool Init() = 0;

	// solve the model
	virtual bool Solve() = 0;

	// get the FE mesh
	FEMesh& GetMesh() { return m_mesh; }

	// facilities for (re)storing the model state data
	virtual void PushState() = 0;
	virtual void PopState () = 0;

public:	// Load curve functions

	//! Add a loadcurve to the model
	void AddLoadCurve(FELoadCurve* plc) { m_LC.push_back(plc); }

	//! get a loadcurve
	FELoadCurve* GetLoadCurve(int i) { return m_LC[i]; }

	//! get the number of loadcurves
	int LoadCurves() { return m_LC.size(); }

public: // material functions

	//! Add a material to the model
	void AddMaterial(FEMaterial* pm) { m_MAT.push_back(pm); }

	//! get the number of materials
	int Materials() { return m_MAT.size(); }

	//! return a pointer to a material
	FEMaterial* GetMaterial(int id) { return m_MAT[id]; }

	//! return the elastic material
	// TODO: this is only a temp solution
	virtual FEElasticMaterial* GetElasticMaterial(int id) = 0;
	virtual FEElasticMaterial* GetElasticMaterial(FEMaterial* pm) = 0;

public: // body force functions

	//! Add a body force to the model
	void AddBodyForce(FEBodyForce* pf) { m_BF.push_back(pf); }

	//! get the number of body forces
	int BodyForces() { return (int) m_BF.size(); }

	//! return a point to a body force
	FEBodyForce* GetBodyForce(int i) { return m_BF[i]; }

	//! see if there are any body forces
	bool HasBodyForces() { return !m_BF.empty();}

public: // analysis step functions
	//! retrieve the number of steps
	int Steps() { return (int) m_Step.size(); }

	//! clear the steps
	void ClearSteps() { m_Step.clear(); }

	//! Add an analysis step
	void AddStep(FEAnalysis* pstep) { m_Step.push_back(pstep); }

	//! Get a particular step
	FEAnalysis* GetStep(int i) { return m_Step[i]; }

public:	// Miscellaneous routines

	//! set callback function
	void AddCallback(FEBIO_CB_FNC pcb, void* pd);

	//! call the callback function
	void DoCallback();

	// get/set global data
	static void SetGlobalConstant(const string& s, double v);
	static double GetGlobalConstant(const string& s);

protected:
	std::vector<FELoadCurve*>	m_LC;	//!< load curve data
	std::vector<FEMaterial*>	m_MAT;	//!< array of materials
	std::vector<FEBodyForce*>	m_BF;	//!< body force data
	std::vector<FEAnalysis*>	m_Step;	//!< array of analysis steps

public:
	// Geometry data
	FEMesh		m_mesh;		//!< the FE mesh

	// Boundary Conditions
	std::vector<FEPrescribedBC*>	m_DC;	//!< prescribed constraints
	std::vector<FENodalForce*>		m_FC;	//!< concentrated nodal loads
	std::vector<FERigidNode*>		m_RN;	//!< rigid nodes

	list<FEBIO_CALLBACK>	m_pcb;	//!< pointer to callback function

protected:
	static std::map<string, double> m_Const;
};
