#pragma once
#include "FEBioTask.h"
#include <vector>

//-----------------------------------------------------------------------------
// This is the FEBio kernel class which manages the interactions between the 
// different executable modules. In particular, it manages the factory classes
// which are responsible for the life of different classes
class FEBioKernel
{
	struct TASK_DESCR
	{
		const char*			sztag;
		FEBioTaskFactory*	ptf;
	};

public:
	void RegisterTask(FEBioTaskFactory* ptf, const char* sztag);
	FEBioTask* CreateTask(const char* sztag);

protected:
	std::vector<TASK_DESCR>		m_Task;
};
