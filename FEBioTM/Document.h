// Document.h: interface for the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
#define AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEM.h"
#include "FEBioLib/log.h"
#include <FL/Fl_Progress.H>
#include <FL/Fl.H>
#include "Task.h"

class CWnd;

//-----------------------------------------------------------------------------
class FETMProgress : public Progress
{
public:
	FETMProgress(Fl_Progress* pw);
	void SetProgress(double f);

protected:
	Fl_Progress*	m_pw;
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

	// add a taks to the document
	CTask* AddTask(const char* szfile);

	// remove a task from the queue
	void RemoveTask(int n);

	// get the number of tasks
	int Tasks() { return (int) m_Task.size(); }

	// get a task
	CTask* GetTask(int i);

	// run a task
	void RunTask(CTask* pt);

	// run all queued tasks
	void RunSession();

	// sessions
	void NewSession();
	bool SaveSession(const char* szfile);
	bool OpenSession(const char* szfile);

protected:
	void RunQueue();

protected:
	vector<CTask*>	m_Task;
};

#endif // !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
