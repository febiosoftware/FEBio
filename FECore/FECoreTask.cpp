#include "stdafx.h"
#include "FECoreTask.h"

//-----------------------------------------------------------------------------
FECoreTask::FECoreTask(FEModel* pfem) : FECoreBase(FETASK_ID) { m_pfem = pfem; }

//-----------------------------------------------------------------------------
FECoreTask::~FECoreTask(void){}
