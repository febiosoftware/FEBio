#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// The FECoreTask class is the base class for all tasks.
// A task is simply the highest level module which defines what the code will do.
class FECORE_EXPORT FECoreTask : public FECoreBase
{
public:
	FECoreTask(FEModel* pfem);
	virtual ~FECoreTask(void);

	//! initialize the task
	//! make sure to call FEModel::Init at some point
	virtual bool Init(const char* szfile) = 0;

	//! Run the task.
	virtual bool Run() = 0;

	//! return the FEModel
	FEModel* GetFEModel() { return m_pfem; }

protected:
	FEModel*	m_pfem;
};
