#pragma once
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
//! This class is the base class of all boundary conditions

//! Specific boundary conditions can be defined be inheriting from this class.
class FECORE_API FEBoundaryCondition : public FEModelComponent
{
	FECORE_SUPER_CLASS

public:
	//! constructor
	FEBoundaryCondition(FEModel* pfem);

	//! desctructor
	virtual ~FEBoundaryCondition(){}
};
