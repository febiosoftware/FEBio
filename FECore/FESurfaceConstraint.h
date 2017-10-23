#pragma once
#include "FENLConstraint.h"

class FESurface;

// Base class for nonlinear constraints that are defined using a surface.
class FESurfaceConstraint : public FENLConstraint
{
public:
	FESurfaceConstraint(FEModel* fem);

	// return the surface
	virtual FESurface* GetSurface() { return 0; }
};
