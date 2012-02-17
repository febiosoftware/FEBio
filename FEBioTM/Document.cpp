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
	// create a new task
	CTask* pt = new CTask;
	pt->SetFileName(szfile);

	// create a text buffer
	Fl_Text_Buffer* pb = new Fl_Text_Buffer;
	int nret = pb->appendfile(szfile);
	pt->SetTextBuffer(pb);

	m_Task.push_back(pt);
	return pt;
}


//-----------------------------------------------------------------------------
bool CDocument::RunTask(int i)
{
	// get the task
	CTask* pt = GetTask(i);

	// create the FEM object
	FEM fem;

	// load the data from file
	if (fem.Input(pt->GetFileName()) == false) return false;

	// initialize FE data
	if (fem.Init() == false) return false;

	// solve the problem
	return fem.Solve();
}
