// Document.cpp: implementation of the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Document.h"
#include "Wnd.h"
#include "MainApp.h"
#include "XMLWriter.h"
#include "FEBioXML/XMLReader.h"
#include <FEBioHeat/FEBioHeat.h>
#include <FEBioMix/FEBioMix.h>
#include "flx_threads.h"
#include <time.h>

extern void InitFEBioLibrary();

static Fl_Thread thread_id = -1;

//-----------------------------------------------------------------------------
void LogBuffer::print(const char* sz)
{
	// obtain a lock before we change the display
	Fl::lock();

	// update the UI
	int N = m_plog->buffer()->length();
	m_plog->insert_position(N);
	m_plog->insert(sz);
	m_plog->show_insert_position();

	// release the lock
	Fl::unlock();

	// notify the main thread to update the display
	Fl::awake((void*) 0);
}

//-----------------------------------------------------------------------------
Progress::Progress(Fl_Progress* pw) : m_pw(pw)
{
	pw->maximum(100.f); 
	pw->minimum(0.f); 
	pw->value(0.f); 
}

//-----------------------------------------------------------------------------
void Progress::SetProgress(double f)
{
	static char sz[256] = {0};

	// obtain a lock before we change the progress bar
	Fl::lock();

	m_pw->value((float) f);
	sprintf(sz, "%d%%", (int) f);
	m_pw->label(sz);

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
	FEBioHeat::InitModule();
	FEBioMix::InitModule();
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
	// get the GUI components
	CWnd* pwnd = (CWnd*) pd;
	CDocument* pdoc = pwnd->GetDocument();
	CTaskBrowser* ptb = pwnd->GetTaskBrowser();

	LogBuffer log(pwnd->GetOutputWnd());
	clog.SetLogStream(&log);

	// continue looping until queue is empty
	int itask = -1;
	do
	{
		// find a queued task
		Fl::lock();
		itask = -1;
		for (int i=0; i<pdoc->Tasks(); ++i)
		{
			if (pdoc->GetTask(i)->GetStatus() == CTask::QUEUED)
			{
				itask = i;
				break;
			}
		}
		Fl::unlock();

		// if task found, run the task
		if (itask != -1)
		{
			time_t time0, time1;

			// get the task
			CTask* pt = pdoc->GetTask(itask);

			// lock the display
			Fl_Progress* pw = 0;
			Fl::lock();
			{
				// add entry to log
				pwnd->AddLogEntry("Running: %s\n", pt->GetFileName());

				time(&time0);
				pwnd->AddLogEntry("\tstart time: %s", ctime(&time0));

				// clear the output wnd
				pwnd->ClearOutputWnd();

				// tell the task browser to track this task
				pw = ptb->TrackTask(itask);
			}
			Fl::unlock();
			Fl::awake((void*)0);

			// set-up the progress tracker
			Progress prg(pw);

			// run the task
			pt->Run(prg);

			// done tracking
			Fl::lock();
			{
				ptb->DoneTracking();
				time(&time1);
				pwnd->AddLogEntry("\tend time  : %s", ctime(&time1));

				double sec = difftime(time1, time0);
				int ih = (int) (sec / 3600.0); sec -= ih*3600;
				int im = (int) (sec / 60.0); sec -= im*60.0;
				int is = (int) sec;
				pwnd->AddLogEntry("\trun time  : %d:%02d:%02d\n\n", ih, im, is);
			}
			Fl::unlock();
			Fl::awake((void*)0);
		}
	}
	while (itask != -1);

	// done
	Fl::lock();
	thread_id = -1;
	Fl::unlock();
	return 0;
}

//-----------------------------------------------------------------------------
void CDocument::RunQueue()
{
	if (thread_id == -1)
	{
		fl_create_thread(thread_id, febio_func, (void*) FLXGetMainWnd());
	}
}

//-----------------------------------------------------------------------------
void CDocument::RunTask(CTask* pt)
{
	// make sure this task is not running
	if (pt->GetStatus() == CTask::RUNNING) return;

	// save the file if necessary
	if (pt->GetStatus() == CTask::MODIFIED) pt->Save();

	// add the file to the queue
	pt->SetStatus(CTask::QUEUED);

	// run the queue
	RunQueue();
}

//-----------------------------------------------------------------------------
void CDocument::RunSession()
{
	// add all tasks to the queue
	for (int i=0; i<Tasks(); ++i)
	{
		// get the next task
		CTask* pt = GetTask(i);

		// save if necessary
		if (pt->GetStatus() == CTask::MODIFIED) pt->Save();

		// only add the task if it is not running
		if (pt->GetStatus() != CTask::RUNNING) pt->SetStatus(CTask::QUEUED);
	}

	// run the queue
	RunQueue();
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
