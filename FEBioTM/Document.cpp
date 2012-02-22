// Document.cpp: implementation of the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Document.h"
#include "Wnd.h"
#include "MainApp.h"

extern void InitFEBioLibrary();
extern const char* wnd_title;

//-----------------------------------------------------------------------------
void CTask::SetFileName(const char* szfile)
{
	m_szfile[0] = 0;
	int l = strlen(szfile)+1;
	assert((l>1) && (l<MAX_FILE));
	if ((l > 1) && (l<MAX_FILE)) strncpy(m_szfile, szfile, l);
}

//-----------------------------------------------------------------------------
const char* CTask::GetFileTitle()
{
	char* c1 = strrchr(m_szfile, '\\');
	char* c2 = strrchr(m_szfile, '/');
	if ((c1 == 0) && (c2 == 0)) return m_szfile;
	if (c1 == 0) return c2+1;
	if (c2 == 0) return c1+1;
	if (c1 < c2) return c2+1; else return c1+1;
}

//-----------------------------------------------------------------------------
void CTask::GetFilePath(char* szpath)
{
	char* c1 = strrchr(m_szfile, '\\');
	char* c2 = strrchr(m_szfile, '/');
	if ((c1 == 0) && (c2 == 0)) strcpy(szpath, m_szfile);
	if ((c1 == 0) || (c2 > c1)) strncpy(szpath, m_szfile, c2 - m_szfile);
	if ((c2 == 0) || (c1 > c2)) strncpy(szpath, m_szfile, c1 - m_szfile);
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
void FETMProgress::SetProgress(double f)
{
	static char sz[1024] = {0};
	int n = (int) f;
	sprintf(sz, "(%d%%) %s - %s", n, m_pTask->GetFileTitle(), wnd_title);
	m_pWnd->label(sz);
	m_pw->value((float) f); 
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

CTask* CDocument::GetTask(int i)
{
	if ((i>=0) && (i<(int)m_Task.size())) return m_Task[i];
	else return 0;
}

//-----------------------------------------------------------------------------
CTask* CDocument::AddTask(const char* szfile)
{
	// create a new task
	CTask* pt = new CTask;
	pt->SetFileName(szfile);

	// create a text buffer
	Fl_Text_Buffer* pb = new Fl_Text_Buffer;
	// TODO: the buffer size has to be a multiple of 100 due to a bug in the Fl_Text_Buffer. 
	//       I need to test whether this bug has been fixed.
	int nret = pb->appendfile(szfile, 100*1024);
	pt->SetTextBuffer(pb);

	// create a log buffer
	pb = new Fl_Text_Buffer;
	pt->SetLogBuffer(pb);

	m_Task.push_back(pt);
	return pt;
}

//-----------------------------------------------------------------------------
void CDocument::RemoveTask(int n)
{
	CTask* pt = GetTask(n);
	assert(pt);
	if (pt == 0) return;
	vector<CTask*>::iterator it = m_Task.begin() + n;
	m_Task.erase(it);
	delete pt;
}

//-----------------------------------------------------------------------------
bool CDocument::RunTask(int i)
{
	// get the task
	CTask* pt = GetTask(i);

	// save the file
	if (pt->GetStatus() == CTask::MODIFIED) pt->Save();

	CWnd* pwnd = FLXGetMainWnd();

	CTaskBrowser* ptb = pwnd->GetTaskBrowser();

	// create a log buffer
	LogBuffer* plog = new LogBuffer(pwnd->GetLogWnd());
	clog.SetLogStream(plog);

	// clear the log
	pt->Clearlog();

	// create the FEM object
	FEM fem;

	// set the default output file names
	char szbase[1024] = {0}, szfile[1024] = {0};
	strcpy(szbase, pt->GetFileName());
	char* ch = strrchr(szbase, '.'); assert(ch);
	if (ch) *ch = 0;
	sprintf(szfile, "%s.log", szbase); fem.SetLogFilename(szfile);
	sprintf(szfile, "%s.plt", szbase); fem.SetPlotFilename(szfile);
	sprintf(szfile, "%s.dmp", szbase); fem.SetDumpFilename(szfile);
	fem.SetInputFilename(pt->GetFileName());

	// load the data from file
	if (fem.Input(pt->GetFileName()) == false)
	{
		pt->SetStatus(CTask::FAILED);
		return false;
	}

	// initialize FE data
	if (fem.Init() == false) 
	{
		pt->SetStatus(CTask::FAILED);
		return false;
	}

	// progress tracker
	FETMProgress prg(pwnd, pt, ptb->TrackSelectedTask());

	pt->SetStatus(CTask::RUNNING);

	// solve the problem
	bool bret = fem.Solve(prg);

	ptb->DoneTracking();

	pt->SetStatus(bret?CTask::COMPLETED:CTask::FAILED);

	// don't forget to clean up
	delete plog;

	// all done!
	return bret;
}
