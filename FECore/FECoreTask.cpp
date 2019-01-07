#include "stdafx.h"
#include "FECoreTask.h"

REGISTER_SUPER_CLASS(FECoreTask, FETASK_ID);

//-----------------------------------------------------------------------------
FECoreTask::FECoreTask(FEModel* fem) : FECoreBase(fem, FETASK_ID) 
{
}

//-----------------------------------------------------------------------------
FECoreTask::~FECoreTask(void)
{
}
