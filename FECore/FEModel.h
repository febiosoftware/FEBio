#pragma once
#include "FEMaterial.h"
#include "FEMesh.h"
#include "LoadCurve.h"
#include "BC.h"
#include "FEBodyForce.h"
#include <vector>

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

public: // body force functions

	//! Add a body force to the model
	void AddBodyForce(FEBodyForce* pf) { m_BF.push_back(pf); }

	//! get the number of body forces
	int BodyForces() { return (int) m_BF.size(); }

	//! return a point to a body force
	FEBodyForce* GetBodyForce(int i) { return m_BF[i]; }

	//! see if there are any body forces
	bool HasBodyForces() { return !m_BF.empty();}

protected:
	std::vector<FELoadCurve*>	m_LC;	//!< load curve data
	std::vector<FEMaterial*>	m_MAT;	//!< array of materials
	std::vector<FEBodyForce*>	m_BF;	//!< body force data

public:
	// Geometry data
	FEMesh		m_mesh;		//!< the FE mesh

	// Boundary Conditions
	std::vector<FEPrescribedBC*>	m_DC;	//!< prescribed constraints
	std::vector<FENodalForce*>		m_FC;	//!< concentrated nodal loads
	std::vector<FERigidNode*>		m_RN;	//!< rigid nodes
};
