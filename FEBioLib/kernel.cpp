#include "stdafx.h"
#include "febio.h"
using namespace std;

void FEBioKernel::RegisterTask(FEBioTaskFactory *ptf, const char *sztag)
{
	TASK_DESCR td;
	td.sztag = sztag;
	td.ptf = ptf;
	m_Task.push_back(td);
}

FEBioTask* FEBioKernel::CreateTask(const char *sztag)
{
	for (int i=0; i<(int) m_Task.size(); ++i)
	{
		TASK_DESCR& td = m_Task[i];
		if (strcmp(td.sztag, sztag) == 0) return td.ptf->CreateTask();
	}
	return 0;
}
