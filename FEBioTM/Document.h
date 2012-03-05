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
	FETMProgress(FEM* pfem, CWnd* pwnd, CTask* pt, Fl_Progress* pw) : m_pfem(pfem), m_pWnd(pwnd), m_pTask(pt), m_pw(pw) { m_bstatus = false; pw->maximum(100.f); pw->minimum(0.f); pw->value(0.f); }
	void SetProgress(double f);

	FEM* GetFEM() { return m_pfem; }
	void SetStatus(bool b) { m_bstatus = b; }
	bool GetStatus() { return m_bstatus; }

protected:
	bool			m_bstatus;
	FEM*			m_pfem;
	Fl_Progress*	m_pw;
	CWnd*			m_pWnd;
	CTask*			m_pTask;
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

	// sessions
	void NewSession();
	bool SaveSession(const char* szfile);
	bool OpenSession(const char* szfile);

protected:
	vector<CTask*>	m_Task;
};

#endif // !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
