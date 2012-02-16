// Document.cpp: implementation of the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Document.h"

//-----------------------------------------------------------------------------
void CTask::SetFileName(const char* szfile)
{
	m_szfile[0] = 0;
	int l = strlen(szfile)+1;
	assert((l>1) && (l<MAX_FILE));
	if ((l > 1) && (l<MAX_FILE)) strncpy(m_szfile, szfile, l);
}

//-----------------------------------------------------------------------------
CDocument::CDocument()
{
}

//-----------------------------------------------------------------------------
CDocument::~CDocument()
{
}

//-----------------------------------------------------------------------------
CTask* CDocument::AddTask(const char* szfile)
{
	CTask* pt = new CTask;
	pt->SetFileName(szfile);
	m_Task.push_back(pt);
	return pt;
}
