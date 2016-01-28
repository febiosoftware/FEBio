#pragma once
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;

//-----------------------------------------------------------------------------
//! Base class for body-loads
//! \todo This is a work in progress
class FEBodyLoad : public FEModelComponent
{
public:
	FEBodyLoad(FEModel* pfem);
	virtual ~FEBodyLoad();

public:
	//! initialization
	virtual bool Init() { return true; }

	//! update
	virtual void Update(){}
};
