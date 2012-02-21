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

	Fl_Text_Display* GetLogWnd() { return m_pLog; }

	CTaskBrowser* GetTaskBrowser() { return m_pTask; }

public:	// --- M E N U   H A N D L E R S ---
	void OnFileOpen(Fl_Widget* pw, void* pd);
	void OnFileExit(Fl_Widget* pw, void* pd);

	void OnRunSelected(Fl_Widget* pw, void* pd);

	// other envent handlers
	void OnSelectFile(Fl_Widget* pw, void* pd);
	void OnSelectTab (Fl_Widget* pw, void* pd);

protected:
	// handle events
	int handle(int nevent);	

protected:
	CDocument*	m_pDoc;	//!< pointer to document

	CMenu*				m_pMenu;	//!< the menu

	Fl_Tabs*			m_pTabs;
	CFileBrowser*		m_pFile;	//!< the file browser
	CTaskBrowser*		m_pTask;	//!< the task browser
	Fl_Text_Editor*		m_pText;	//!< Text display
	Fl_Text_Display*	m_pLog;		//!< log text display
};

#endif // !defined(AFX_WND_H__793D79A3_EBE4_4660_8EE3_0016B7467520__INCLUDED_)
