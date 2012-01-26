#pragma once
#include "FECore/FETrussDomain.h"
#include "FEBioLib/FEElasticDomain.h"

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

	//! calculates the residual
	void Residual(FENLSolver* psolver, vector<double>& R);

protected:
	//! calculates the truss element stiffness matrix
	void ElementStiffness(int iel, matrix& ke);

	//! Calculates the internal stress vector for solid elements
	void InternalForces(FETrussElement& el, vector<double>& fe);

	//! update the truss stresses
	void UpdateStresses(FEModel& fem);
};
