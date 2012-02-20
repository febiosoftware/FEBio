// Document.cpp: implementation of the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Document.h"
#include "Wnd.h"
#include "MainApp.h"

extern void InitFEBioLibrary();

//-----------------------------------------------------------------------------
void CTask::SetFileName(const char* szfile)
{
	m_szfile[0] = 0;
	int l = strlen(szfile)+1;
	assert((l>1) && (l<MAX_FILE));
	if ((l > 1) && (l<MAX_FILE)) strncpy(m_szfile, szfile, l);
}

//-----------------------------------------------------------------------------
void LogBuffer::print(const char* sz)
{
	m_plog->insert(sz);
	m_plog->show_insert_position();
//	m_plog->redraw();
	Fl::flush();
}

//-----------------------------------------------------------------------------
CDocument::CDocument()
{
	// initialize FEBio library
	InitFEBioLibrary();
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

	// create a log buffer
	pb = new Fl_Text_Buffer;
	pt->SetLogBuffer(pb);

	m_Task.push_back(pt);
	return pt;
}


//-----------------------------------------------------------------------------
bool CDocument::RunTask(int i)
{
	// get the task
	CTask* pt = GetTask(i);

	CWnd* pwnd = FLXGetMainWnd();

	// create a log buffer
	LogBuffer* plog = new LogBuffer(pwnd->GetLogWnd());
	clog.SetLogStream(plog);

	// create the FEM object
	FEM fem;

	// load the data from file
	if (fem.Input(pt->GetFileName()) == false) return false;

	// initialize FE data
	if (fem.Init() == false) return false;

	// progress tracker
	FETMProgress prg(fem);

	// solve the problem
	bool bret = fem.Solve(prg);

	// don't forget to clean up
	delete plog;

	// all done!
	return bret;
}
