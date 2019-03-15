#pragma once
#include "FERefineMesh.h"

class FETetRefine : public FERefineMesh
{
public:
	FETetRefine(FEModel* fem);

	bool Apply(int iteration);

protected:
	bool DoTetRefinement(FEModel& fem);
};
