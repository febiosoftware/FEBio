#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! Class implementing selective reduced integration for the evaluation of the internal force vector
class FESRIElasticSolidDomain : public FEElasticSolidDomain
{
public:
	FESRIElasticSolidDomain(FEModel* pfem);

public:
	//! internal stress forces
	virtual void ElementInternalForce(FESolidElement& el, vector<double>& fe);
};
