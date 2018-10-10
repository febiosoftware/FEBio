#pragma once
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
//! This class is the base class of all boundary conditions

//! Specific boundary conditions can be defined be inheriting from this class.
class FECORE_API FEBoundaryCondition : public FEModelComponent
{
	DECLARE_SUPER_CLASS(FEBC_ID);

public:
	//! constructor
	FEBoundaryCondition(SUPER_CLASS_ID sid, FEModel* pfem);

	//! desctructor
	virtual ~FEBoundaryCondition(){}
};
