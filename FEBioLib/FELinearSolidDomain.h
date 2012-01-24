#pragma once

#include "FECore/FESolidDomain.h"
#include "FECore/FENLSolver.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
//! Class describing a linear elastic solid domain
class FELinearSolidDomain : public FESolidDomain
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

	//! Build the stiffness matrix
	void StiffnessMatrix(FENLSolver* psolver);

	//! Update the element stresses
	void UpdateStresses(FEModel& fem);

	// Calculate the RHS vector
	void RHS(FENLSolver* psolver, vector<double>& R);

	//! reset element data
	void Reset();

	//! initialize elements
	void InitElements();

protected:
	void InitialStress(FESolidElement& el, vector<double>& fe);
	void InternalForce(FESolidElement& el, vector<double>& fe);

protected:
	void ElementStiffness(FESolidElement& el, matrix& ke);
};
