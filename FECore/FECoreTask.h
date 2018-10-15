#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// The FECoreTask class is the base class for all tasks.
// A task is simply the highest level module which defines what the code will do.
class FECORE_API FECoreTask : public FECoreBase
{
	DECLARE_SUPER_CLASS(FETASK_ID);

public:
	FECoreTask(FEModel* fem);
	virtual ~FECoreTask(void);

	//! initialize the task
	//! make sure to call FEModel::Init at some point
	virtual bool Init(const char* szfile) = 0;

	//! Run the task.
	virtual bool Run() = 0;
};
