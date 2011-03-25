#pragma once
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
//! This is the base class for all loads that are applied to surfaces
class FESurfaceLoad : public FEBoundaryCondition
{
public:
	FESurfaceLoad(void);
	~FESurfaceLoad(void);
};
