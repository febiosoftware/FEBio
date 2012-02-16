// MainApp.h: interface for the CMainApp class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MAINAPP_H__5432BE41_2F50_46E1_80E3_E83BAC5BF8D8__INCLUDED_)
#define AFX_MAINAPP_H__5432BE41_2F50_46E1_80E3_E83BAC5BF8D8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <Flx_App.h>

#include "Wnd.h"
#include "Document.h"

class CMainApp : public Flx_App 
{
public:
	CMainApp();
	virtual ~CMainApp();

	bool Init(); // initializes the application
	bool Load(const char* szfilename);
	int Run();	// runs the application

	CDocument* GetDocument() { return m_pDoc; }
	CWnd* GetMainWnd() { return m_pMainWnd; }


protected:
	CWnd*		m_pMainWnd;	// main window
	CDocument*	m_pDoc;		// document class
};

CMainApp* FLX_GetMainApp();
CWnd*	FLX_GetMainWnd();

#endif // !defined(AFX_MAINAPP_H__5432BE41_2F50_46E1_80E3_E83BAC5BF8D8__INCLUDED_)
