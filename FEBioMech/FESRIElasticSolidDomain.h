#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! Class implementing selective reduced integration for the evaluation of the internal force vector
class FESRIElasticSolidDomain : public FEElasticSolidDomain
{
public:
	FESRIElasticSolidDomain(FEMesh* pmesh, FEMaterial* pmat);

	//! create a copy
	FEDomain* Copy();

public:
	//! internal stress forces
	virtual void ElementInternalForce(FESolidElement& el, vector<double>& fe);
};
