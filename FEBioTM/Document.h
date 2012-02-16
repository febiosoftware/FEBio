// Document.h: interface for the CDocument class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
#define AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEM.h"

//-----------------------------------------------------------------------------
class CTask
{
	enum {MAX_FILE = 512};

public:
	void SetFileName(const char* szfile);
	const char* GetFileName() { return m_szfile; }

protected:
	char	m_szfile[MAX_FILE];
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
	CTask* Task(int i) { return m_Task[i]; }

protected:
	vector<CTask*>	m_Task;
	FEM		m_fem;
};

#endif // !defined(AFX_DOCUMENT_H__E699CE33_76BC_46FB_8CFC_4FA83D106B4C__INCLUDED_)
