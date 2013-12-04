#pragma once
#include "DumpFile.h"
#include "FEParameterList.h"

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;

//-----------------------------------------------------------------------------
//! Base class for body-loads
//! \todo This is a work in progress
class FEBodyLoad : public FEParamContainer
{
public:
	FEBodyLoad(FEModel* pfem);
	virtual ~FEBodyLoad();

	//! Get the FE model
	FEModel* GetFEModel() { return m_pfem; }

public:
	//! initialization
	virtual bool Init() { return true; }

	//! serialize data to archive
	virtual void Serialize(DumpFile& ar){}

	//! update
	virtual void Update(){}

private:
	FEModel* m_pfem;
};
