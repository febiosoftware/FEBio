#pragma once
#include "FECore/FETrussDomain.h"
#include "FEElasticDomain.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D truss elements
class FEElasticTrussDomain : public FETrussDomain, public FEElasticDomain
{
public:
	//! Constructor
	FEElasticTrussDomain(FEMesh* pm, FEMaterial* pmat) : FETrussDomain(FE_TRUSS_DOMAIN, pm, pmat) {}

	//! Clone this domain
	FEDomain* Clone();

	//! copy operator
	FEElasticTrussDomain& operator = (FEElasticTrussDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! Reset data
	void Reset();

	//! Initialize elements
	void InitElements();

	//! Unpack truss element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public: // overloads from FEElasticDomain

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver);

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FENLSolver* psolver, vector<double>& R, vector<double>& F) { assert(false); }

	//! calculates the residual
//	void Residual(FENLSolver* psolver, vector<double>& R);

	//! internal stress forces
	void InternalForces(FENLSolver* psolver, vector<double>& R);

	//! calculate body force (TODO: implement this)
	void BodyForce(FENLSolver* psolver, FEBodyForce& bf, vector<double>& R) { assert(false); }

	//! update the truss stresses
	void UpdateStresses(FEModel& fem);

protected:
	//! calculates the truss element stiffness matrix
	void ElementStiffness(int iel, matrix& ke);

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForces(FETrussElement& el, vector<double>& fe);
};
