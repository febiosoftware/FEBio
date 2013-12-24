#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// The FECoreTask class is the base class for all tasks.
// A task is simply the highest level module which defines what the code will do.
class FECoreTask : public FECoreBase
{
public:
	FECoreTask(FEModel* pfem);
	virtual ~FECoreTask(void);

	//! Run the task.
	virtual bool Run (const char* szfile) = 0;

protected:
	FEModel*	m_pfem;
};
