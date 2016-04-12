#pragma once
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
//! This class is the base class of all boundary conditions

//! Specific boundary conditions can be defined be inheriting from this class.
class FEBoundaryCondition : public FEModelComponent
{
public:
	//! constructor
	FEBoundaryCondition(SUPER_CLASS_ID sid, FEModel* pfem);

	//! desctructor
	virtual ~FEBoundaryCondition(){}
};
