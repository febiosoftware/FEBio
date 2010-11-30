#include "StdAfx.h"
#include "FEPlotDataFactory.h"

FEPlotDataFactory* FEPlotDataFactory::m_pThis = 0;

FEPlotDataFactory* FEPlotDataFactory::GetInstance()
{
	if (m_pThis == 0) m_pThis = new FEPlotDataFactory;
	return m_pThis;
}

FEPlotData* FEPlotDataFactory::Create(const char* sz)
{
	list<ClassDescriptor*> cdl = GetInstance()->m_list;

	list<ClassDescriptor*>::iterator pi;
	for (pi=cdl.begin(); pi != cdl.end(); ++pi)
	{
		if (strcmp((*pi)->GetName(), sz) == 0) return (*pi)->Create();
	}
	return 0;
}
