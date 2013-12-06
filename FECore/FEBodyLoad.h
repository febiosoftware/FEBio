#pragma once
#include "DumpFile.h"
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

	//! serialize data to archive
	virtual void Serialize(DumpFile& ar){}

	//! update
	virtual void Update(){}
};
