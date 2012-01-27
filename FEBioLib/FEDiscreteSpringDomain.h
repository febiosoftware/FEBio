#pragma once
#include "FECore/FEDiscreteDomain.h"
#include "FEElasticDomain.h"

//-----------------------------------------------------------------------------
//! domain for discrete elements
class FEDiscreteSpringDomain : public FEDiscreteDomain, public FEElasticDomain
{
public:
	//! constructor
	FEDiscreteSpringDomain(FEMesh* pm, FEMaterial* pmat) : FEDiscreteDomain(FE_DISCRETE_DOMAIN, pm, pmat) {}

	//! Clone this domain
	FEDomain* Clone();

	//! Unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! Serialize data to archive
	void Serialize(DumpFile& ar);

public: // overridden from FEElasticDomain

	//! calculate stiffness matrix
	void StiffnessMatrix(FENLSolver* psolver);

	//! calculate residual
	void Residual(FENLSolver* psolver, vector<double>& R);

protected:
	void UpdateStresses(FEModel& fem){}	// not used for discrete springs
};
