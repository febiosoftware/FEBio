#include "stdafx.h"
#include "febio.h"
using namespace std;

FEBioKernel* FEBioKernel::m_pKernel = 0;

FEBioKernel& FEBioKernel::GetInstance()
{
	if (m_pKernel == 0) m_pKernel = new FEBioKernel;
	return *m_pKernel;
}

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

void FEBioKernel::RegisterMaterial(FEMaterialFactory* pmf, const char* sztag)
{
	MATERIAL_DESCR cd;
	cd.sztag = sztag;
	cd.pfac = pmf;
	m_Mat.push_back(cd);
}

FEMaterial* FEBioKernel::CreateMaterial(const char* sztag, FEModel* pfem)
{
	for (int i=0; i<(int) m_Mat.size(); ++i)
	{
		MATERIAL_DESCR& cd = m_Mat[i];
		if (strcmp(cd.sztag, sztag) == 0) return cd.pfac->Create(pfem);
	}
	return 0;
}
