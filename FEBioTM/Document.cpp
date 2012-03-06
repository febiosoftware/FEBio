// Document.cpp: implementation of the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Document.h"
#include "Wnd.h"
#include "MainApp.h"
#include "XMLWriter.h"
#include "FEBioXML/XMLReader.h"
#include <FL/threads.h>

extern void InitFEBioLibrary();

//-----------------------------------------------------------------------------
void LogBuffer::print(const char* sz)
{
	// obtain a lock before we change the display
	Fl::lock();

	// update the UI
	m_plog->insert(sz);
	m_plog->show_insert_position();

	// release the lock
	Fl::unlock();

	// notify the main thread to update the display
	Fl::awake((void*) 0);
}

//-----------------------------------------------------------------------------
FETMProgress::FETMProgress(FEM* pfem, CWnd* pwnd, CTask* pt, Fl_Progress* pw) : m_pfem(pfem), m_pWnd(pwnd), m_pTask(pt), m_pw(pw)
{
	pw->maximum(100.f); 
	pw->minimum(0.f); 
	pw->value(0.f); 
}

//-----------------------------------------------------------------------------
void FETMProgress::SetProgress(double f)
{
	// obtain a lock before we change the progress bar
	Fl::lock();

	m_pw->value((float) f);

	// releas the lock
	Fl::unlock();

	// notify the main thread to update the display
	Fl::awake((void*) 0);
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
void CDocument::NewSession()
{
	while (m_Task.empty() == false) RemoveTask(0);
}

//-----------------------------------------------------------------------------
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
	int nret = pb->appendfile(szfile);
	pt->SetTextBuffer(pb);

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
void* febio_func(void* pd)
{
	FETMProgress* prg = (FETMProgress*)(pd);
	FEM* pfem = prg->GetFEM();

	CTask* pt = prg->GetTask();
	Fl::lock();
	CTask::m_prun = pt;
	pt->SetStatus(CTask::RUNNING);
	Fl::unlock();

	bool bret = pfem->Solve(*prg);

	Fl::lock();
	pt->SetStatus(bret?CTask::COMPLETED:CTask::FAILED);
	CWnd* pwnd = prg->GetWnd();
	CTaskBrowser* ptb = pwnd->GetTaskBrowser();
	ptb->DoneTracking();
	Fl::unlock();

	Fl::awake((void*)0);

	// clean up
	delete pfem;
	delete prg;

	// done
	return 0;
}

//-----------------------------------------------------------------------------
void CDocument::RunTask(CTask* pt)
{
	// make sure no task is running
	if (CTask::m_prun != 0) return;

	// save the file
	if (pt->GetStatus() == CTask::MODIFIED) pt->Save();

	CWnd* pwnd = FLXGetMainWnd();

	CTaskBrowser* ptb = pwnd->GetTaskBrowser();

	// create a log buffer
	static LogBuffer* plog = new LogBuffer(pwnd->GetLogWnd());
	clog.SetLogStream(plog);

	// clear the output window
	pwnd->ClearOutputWnd();

	// create the FEM object
	FEM* pfem = new FEM(pt);

	// set the default output file names
	char szbase[1024] = {0}, szfile[1024] = {0};
	strcpy(szbase, pt->GetFileName());
	char* ch = strrchr(szbase, '.'); assert(ch);
	if (ch) *ch = 0;
	sprintf(szfile, "%s.log", szbase); pfem->SetLogFilename(szfile);
	sprintf(szfile, "%s.plt", szbase); pfem->SetPlotFilename(szfile);
	sprintf(szfile, "%s.dmp", szbase); pfem->SetDumpFilename(szfile);
	pfem->SetInputFilename(pt->GetFileName());

	// load the data from file
	if (pfem->Input(pt->GetFileName()) == false)
	{
		pt->SetStatus(CTask::FAILED);
		delete pfem;
		return;
	}

	// initialize FE data
	if (pfem->Init() == false) 
	{
		pt->SetStatus(CTask::FAILED);
		delete pfem;
		return;
	}

	// progress tracker
	FETMProgress* prg = new FETMProgress(pfem, pwnd, pt, ptb->TrackSelectedTask());

	// solve the problem
	Fl_Thread thread_id;
	fl_create_thread(thread_id, febio_func, prg);
}

//-----------------------------------------------------------------------------
bool CDocument::SaveSession(const char* szfile)
{
	XMLWriter xml;
	xml.open(szfile);
	xml.add_branch("febio_tm_session");
	for (int i=0; i<Tasks(); ++i)
	{
		CTask* pt = GetTask(i);
		XMLElement e;
		e.name("file");
		e.add_attribute("name", pt->GetFileName());
		xml.add_empty(e);
	}
	xml.close_branch();
	xml.close();
	return true;
}

//-----------------------------------------------------------------------------
bool CDocument::OpenSession(const char* szfile)
{
	XMLReader xml;
	if (xml.Open(szfile) == false) return false;

	XMLTag tag;
	if (xml.FindTag("febio_tm_session", tag) == false) return false;

	NewSession();

	++tag;
	do
	{
		if (tag == "file")
		{
			const char* szfile = tag.AttributeValue("name");
			AddTask(szfile);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	xml.Close();
	return true;
}
