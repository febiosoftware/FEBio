#pragma once

#include "FECore/FESolidDomain.h"
#include "FECore/FENLSolver.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
class FELinearElasticDomain
{
public:
	virtual ~FELinearElasticDomain(){}
	virtual void StiffnessMatrix(FENLSolver* psolver) = 0;
	virtual void RHS(FENLSolver* psolver, vector<double>& R) = 0;
	virtual void UpdateStresses(FEModel& fem) = 0;
};

//-----------------------------------------------------------------------------
//! Class describing a linear elastic solid domain
class FELinearSolidDomain : public FESolidDomain, public FELinearElasticDomain
{
public:
	//! constructor
	FELinearSolidDomain(FEMesh* pm, FEMaterial* pmat) : FESolidDomain(FE_LINEAR_SOLID_DOMAIN, pm, pmat) {}

	//! Clone the data
	FEDomain* Clone();

	//! Initialization
	bool Initialize(FEModel& fem);

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! reset element data
	void Reset();

	//! initialize elements
	void InitElements();

public: // overrides from FELinearElasticDomain

	//! Build the stiffness matrix
	void StiffnessMatrix(FENLSolver* psolver);

	// Calculate the RHS vector
	void RHS(FENLSolver* psolver, vector<double>& R);

	//! Update the element stresses
	void UpdateStresses(FEModel& fem);

protected:
	void InitialStress(FESolidElement& el, vector<double>& fe);
	void InternalForce(FESolidElement& el, vector<double>& fe);

	void ElementStiffness(FESolidElement& el, matrix& ke);
};
