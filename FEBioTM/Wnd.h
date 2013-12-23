// Wnd.h: interface for the CWnd class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_WND_H__793D79A3_EBE4_4660_8EE3_0016B7467520__INCLUDED_)
#define AFX_WND_H__793D79A3_EBE4_4660_8EE3_0016B7467520__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <Flx_Wnd.h>
#include "Document.h"
#include "Menu.h"
#include "FileBrowser.h"
#include "TaskBrowser.h"
#include "SettingsView.h"
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Text_Editor.H>
#include <FL/Fl_Tabs.H>

//-----------------------------------------------------------------------------
class CWnd : public Flx_Wnd
{
public:
	// constructor
	CWnd(int w, int h, const char* sztitle, CDocument* pdoc);

	// destructor
	virtual ~CWnd();

	// retrieve the document
	CDocument* GetDocument() { return m_pDoc; }

	// open a file
	bool OpenFile(const char* szfile);

	// update the GUI
	void Update();

	Fl_Text_Display* GetOutputWnd() { return m_pOut; }
	Fl_Text_Display* GetLogWnd() { return m_pLog; }

	CTaskBrowser* GetTaskBrowser() { return m_pTask; }

	CTask* GetSelectedTask();

	CTask* GetFirstSelectedTask();
	CTask* GetNextSelectedTask(CTask* pt);
	bool IsTaskSelected(int n) { return m_pTask->IsTaskSelected(n); }

	void ClearOutputWnd();
	void ClearLogWnd();
	void ClearTestWnd();

	void AddLogEntry(const char* sz, ...);
	void AddTestEntry(const char* sz, ...);
	void SetTestFormat(int i) { m_nTestFmt = i; }

public:	// --- M E N U   H A N D L E R S ---
	void OnFileOpen       (Fl_Widget* pw, void* pd);
	void OnFileSave       (Fl_Widget* pw, void* pd);
	void OnFileSaveAs     (Fl_Widget* pw, void* pd);
	void OnFileRevert     (Fl_Widget* pw, void* pd);
	void OnFileClose      (Fl_Widget* pw, void* pd);
	void OnFileCloseAll   (Fl_Widget* pw, void* pd);
	void OnFileOpenSession(Fl_Widget* pw, void* pd);
	void OnFileSaveSession(Fl_Widget* pw, void* pd);
	void OnFileExit       (Fl_Widget* pw, void* pd);

	void OnEditSelectAll  (Fl_Widget* pw, void* pd);
	void OnEditFind       (Fl_Widget* pw, void* pd);
	void OnEditFindAgain  (Fl_Widget* pw, void* pd);
	void OnEditGoToLine   (Fl_Widget* pw, void* pd);
	void OnEditFilter     (Fl_Widget* pw, void* pd);
	void OnEditClearFilter(Fl_Widget* pw, void* pd);

	void OnRunSelected      (Fl_Widget* pw, void* pd);
	void OnRunSession       (Fl_Widget* pw, void* pd);
	void OnRunCancelSelected(Fl_Widget* pw, void* pd);
	void OnRunCancelAll     (Fl_Widget* pw, void* pd);

	void OnToolsCreateTest(Fl_Widget* pw, void* pd);
	void OnToolsRunTest   (Fl_Widget* pw, void* pd);
	void OnToolsSaveData  (Fl_Widget* pw, void* pd);
	void OnToolsLoadData  (Fl_Widget* pw, void* pd);

	// other envent handlers
	void OnSelectFile(Fl_Widget* pw, void* pd);
	void OnSelectTab (Fl_Widget* pw, void* pd);
	void OnChangeText(Fl_Widget* pw, void* pd);

	void SelectFile();

protected:
	// handle events
	int handle(int nevent);	

protected:
	CDocument*	m_pDoc;	//!< pointer to document
	char		m_szfind[1024];	//!< last find string
	bool		m_bcase;		//!< case sensitive search

	CMenu*				m_pMenu;	//!< the menu

	Fl_Tabs*			m_pTabs;
	CFileBrowser*		m_pFile;	//!< the file browser
	CTaskBrowser*		m_pTask;	//!< the task browser
	Fl_Text_Editor*		m_pText;	//!< input file display
	Fl_Text_Display*	m_pOut;		//!< output text display
	Fl_Text_Display*	m_pLog;		//!< run-log display
	Fl_Text_Display*	m_pTest;	//!< Test View
	CSettingsView*		m_pOps;		//!< file settings view

	int					m_nTestFmt;	//!< format for test output
	Fl_Text_Buffer*		m_pTestAtt;	//!< test attributes

	Fl_Text_Display*	m_pSel;		//!< Target for editing
};

#endif // !defined(AFX_WND_H__793D79A3_EBE4_4660_8EE3_0016B7467520__INCLUDED_)
