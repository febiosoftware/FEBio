#include "stdafx.h"
#include "febio.h"
using namespace std;

void FEBioKernel::RegisterTask(FEBioTaskFactory *ptf, const char *sztag)
{
	TASK_DESCR td;
	td.sztag = sztag;
	td.pfac = ptf;
	m_Task.push_back(td);
}

FEBioTask* FEBioKernel::CreateTask(const char *sztag, FEModel* pfem)
{
	for (int i=0; i<(int) m_Task.size(); ++i)
	{
		TASK_DESCR& td = m_Task[i];
		if (strcmp(td.sztag, sztag) == 0) return td.pfac->Create(pfem);
	}
	return 0;
}

void FEBioKernel::RegisterBodyForce(FEBodyForceFactory *pfac, const char *sztag)
{
	BODY_FORCE_DESCR cd;
	cd.sztag = sztag;
	cd.pfac = pfac;
	m_BF.push_back(cd);
}

FEBodyForce* FEBioKernel::CreateBodyForce(const char *sztag, FEModel* pfem)
{
	for (int i=0; i<(int) m_Task.size(); ++i)
	{
		BODY_FORCE_DESCR& cd = m_BF[i];
		if (strcmp(cd.sztag, sztag) == 0) return cd.pfac->Create(pfem);
	}
	return 0;
}
