// Document.cpp: implementation of the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Document.h"
#include "Wnd.h"
#include "MainApp.h"
#include "XMLWriter.h"
#include "FEBioXML/XMLReader.h"
#include "flx_threads.h"
#include <time.h>
#include <assert.h>

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
	m_ptest = 0;

	// initialize FEBio library
	InitFEBioLibrary();
}

//-----------------------------------------------------------------------------
CDocument::~CDocument()
{
	if (m_ptest) delete m_ptest;
}

//-----------------------------------------------------------------------------
void CDocument::AddTest(CTest* ptest)
{ 
	if (m_ptest) delete m_ptest;
	m_ptest = ptest; 
}

//-----------------------------------------------------------------------------
void CDocument::RunTest()
{
	CWnd* pwnd = FLXGetMainWnd();
	if (m_ptest)
	{
		const char* sz[] = {"Failed", "Passed", "No Data"};
		const int M = 12;
		char szret[M], sztime[M], sziter[M], szrhs[M], szref[M];

		m_ptest->Run();

		pwnd->SetTestFormat(0);
		pwnd->AddTestEntry(" Name                           | Return code| Timesteps  | Iterations | RHS evals  |  Reforms   | Status\n");
		pwnd->AddTestEntry("-------------------------------------------------------------------------------------------------------------------------\n");
		int N = m_session.Tasks();
		for (int i=0; i<N; ++i)
		{
			CTask& task = *m_session.GetTask(i);
			CTest::STATS data =  m_ptest->GetStats(task);
			CTask::STATS stat = task.m_stats;
			if (stat.nreturn == data.nreturn) sprintf(szret , "%d", stat.nreturn); else sprintf(szret , "%d(%d)", stat.nreturn, data.nreturn);
			if (stat.ntime   == data.ntime  ) sprintf(sztime, "%d", stat.ntime  ); else sprintf(sztime, "%d(%d)", stat.ntime  , data.ntime  );
			if (stat.niters  == data.niters ) sprintf(sziter, "%d", stat.niters ); else sprintf(sziter, "%d(%d)", stat.niters , data.niters );
			if (stat.nrhs    == data.nrhs   ) sprintf(szrhs , "%d", stat.nrhs   ); else sprintf(szrhs , "%d(%d)", stat.nrhs   , data.nrhs   );
			if (stat.nreform == data.nreform) sprintf(szref , "%d", stat.nreform); else sprintf(szref , "%d(%d)", stat.nreform, data.nreform);

			int nresult = m_ptest->GetResult(i);
			if (nresult == 0) pwnd->SetTestFormat(1); else pwnd->SetTestFormat(0);
			pwnd->AddTestEntry("%-32s|%12s|%12s|%12s|%12s|%12s| %s\n", task.GetFileTitle(), szret, sztime, sziter, szrhs, szref, sz[nresult]);
		}
	}
}

//-----------------------------------------------------------------------------
void CDocument::NewSession()
{
	if (m_ptest) delete m_ptest;
	m_session.Clear();
}

//-----------------------------------------------------------------------------
void* febio_func(void* pd)
{
	// get the GUI components
	CWnd* pwnd = (CWnd*) pd;
	CDocument* pdoc = pwnd->GetDocument();
	CTaskBrowser* ptb = pwnd->GetTaskBrowser();

	TMSession& session = pdoc->GetSession();

	LogBuffer log(pwnd->GetOutputWnd());
	felog.SetLogStream(&log);

	// continue looping until queue is empty
	int itask = -1;
	do
	{
		// find a queued task
		Fl::lock();
		itask = -1;
		for (int i=0; i<session.Tasks(); ++i)
		{
			if (session.GetTask(i)->GetStatus() == CTask::QUEUED)
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
			CTask* pt = session.GetTask(itask);

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
	for (int i=0; i<m_session.Tasks(); ++i)
	{
		// get the next task
		CTask* pt = m_session.GetTask(i);

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
	{
		xml.add_branch("Files");
		for (int i=0; i<m_session.Tasks(); ++i)
		{
			CTask* pt = m_session.GetTask(i);
			XMLElement e;
			e.name("file");
			e.add_attribute("name", pt->GetFileName());
			xml.add_empty(e);
		}
		xml.close_branch();
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
		if (tag == "Files")
		{
			++tag;
			do
			{
				if (tag == "file")
				{
					const char* szfile = tag.AttributeValue("name");
					m_session.AddTask(szfile);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	xml.Close();
	return true;
}
