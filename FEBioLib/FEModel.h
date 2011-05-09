#pragma once
#include "LoadCurve.h"
#include "BC.h"
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

protected:
	std::vector<FELoadCurve*>	m_LC;	//!< load curve data

public:
	// Boundary Conditions
	std::vector<FEPrescribedBC*>	m_DC;	//!< prescribed constraints
	std::vector<FENodalForce*>		m_FC;	//!< concentrated nodal loads
	std::vector<FERigidNode*>		m_RN;	//!< rigid nodes
};
