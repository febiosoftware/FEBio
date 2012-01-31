#pragma once
#include "FECore/FEShellDomain.h"
#include "FEElasticDomain.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomain : public FEShellDomain, public FEElasticDomain
{
public:
	FEElasticShellDomain(FEMesh* pm, FEMaterial* pmat) : FEShellDomain(FE_SHELL_DOMAIN, pm, pmat) {}

	//! TODO: do I really need this?
	FEElasticShellDomain& operator = (FEElasticShellDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! Clone this domain
	FEDomain* Clone();

	//! Reset element data
	void Reset();

	//! Initialize elements
	void InitElements();

	//! Initialize domain
	bool Initialize(FEModel& fem);

	//! Serialize domain data to archive
	void Serialize(DumpFile& ar);

	//! Unpack shell element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public: // overrides from FEElasticDomain

	//! calculates the residual
//	void Residual(FENLSolver* psolver, vector<double>& R);

	//! internal stress forces
	void InternalForces(FENLSolver* psolver, vector<double>& R);

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FENLSolver* psolver, vector<double>& R, vector<double>& F) { assert(false); }

	//! calculate body force
	void BodyForce(FENLSolver* psolver, FEBodyForce& bf, vector<double>& R);

	// update stresses
	void UpdateStresses(FEModel& fem);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver);

	// inertial stiffness (TODO: implement this)
	void InertialStiffness   (FENLSolver* psolver) { assert(false); }

	// body force stiffness (TODO: implement this)
	void BodyForceStiffness  (FENLSolver* psolver) { assert(false); }

public:

	// --- S T I F F N E S S --- 

	//! calculates the shell element stiffness matrix
	void ElementStiffness(int iel, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for shell elements
	void ElementInternalForce(FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe);
};
