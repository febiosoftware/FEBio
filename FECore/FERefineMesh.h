#pragma once
#include "fecore_api.h"
class FEModel;

class FECORE_API FERefineMesh
{
public:
	FERefineMesh();

	bool Apply(FEModel& fem);

protected:
	bool DoMeshRefinement(FEModel& fem);
};
