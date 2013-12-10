// Document.h: interface for the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
#define AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/log.h"
#include <FL/Fl_Progress.H>
#include <FL/Fl.H>
#include "Session.h"
#include "Test.h"

class CWnd;
class FEBioModel;

//-----------------------------------------------------------------------------
class Progress
{
public:
	Progress(Fl_Progress* pw);
	void SetProgress(double f);

	Fl_Progress* GetWidget() { return m_pw; }

	void SetTask(CTask* pt);
	CTask* GetTask() { return m_pt; }

public:
	Fl_Progress*	m_pw;
	CTask*			m_pt;
	FEBioModel*		m_pfem;
};

//-----------------------------------------------------------------------------
// class that will direct log output to the Log window
class LogBuffer : public LogStream
{
public:
	LogBuffer(Fl_Text_Display* pb) : m_plog(pb) {}
	void print(const char* sz);
private:
	Fl_Text_Display*	m_plog;
};

//-----------------------------------------------------------------------------
// Document class - stores all data
//
class CDocument  
{
public:
	// constructor/destructor
	CDocument();
	virtual ~CDocument();

	// get the session
	TMSession& GetSession() { return m_session; }

	// run a task
	void RunTask(CTask* pt);

	// run all queued tasks
	void RunSession();

	// sessions
	void NewSession();
	bool SaveSession(const char* szfile);
	bool OpenSession(const char* szfile);

	// Add a test
	void AddTest(CTest* ptest);

	// run test
	void RunTest();

	// get the active test
	CTest* GetActiveTest() { return m_ptest; }

public:
	void OnTimer();

protected:
	void RunQueue();

protected:
	TMSession	m_session;
	CTest*		m_ptest;
};

#endif // !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
