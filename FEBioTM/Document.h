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
#include <FL/Fl_Text_Display.H>

//-----------------------------------------------------------------------------
class CTask
{
	enum {MAX_FILE = 512};

public:
	enum { QUEUED, RUNNING, COMPLETED };

public:
	CTask() { m_szfile[0] = 0; m_pfile = 0; m_plog = 0; m_nstatus = QUEUED; }
	~CTask() { delete m_pfile; delete m_plog; }

	void SetFileName(const char* szfile);
	const char* GetFileName() { return m_szfile; }

	void SetTextBuffer(Fl_Text_Buffer* pb) { m_pfile = pb; }
	Fl_Text_Buffer* GetTextBuffer() { return m_pfile; }

	void SetLogBuffer(Fl_Text_Buffer* pb) { m_plog = pb; }
	Fl_Text_Buffer* GetLogBuffer() { return m_plog; }

	void Clearlog()
	{
		m_plog->select(0, m_plog->length());
		m_plog->remove_selection();	
	}

	void SetStatus(int n) { m_nstatus = n; }
	int GetStatus() { return m_nstatus; }

	void Save() { m_pfile->savefile(m_szfile); }

protected:
	char			m_szfile[MAX_FILE];		//!< file name
	Fl_Text_Buffer*	m_pfile;				//!< text buffer for editing
	Fl_Text_Buffer*	m_plog;					//!< log buffer
	int				m_nstatus;				//!< status
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

	// get the number of tasks
	int Tasks() { return (int) m_Task.size(); }

	// get a task
	CTask* GetTask(int i) { return m_Task[i]; }

	// run a task
	bool RunTask(int i);

protected:
	vector<CTask*>	m_Task;
};

#endif // !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
