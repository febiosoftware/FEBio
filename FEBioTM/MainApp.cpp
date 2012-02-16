// MainApp.cpp: implementation of the CMainApp class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MainApp.h"
#include <FL/fl_ask.H>
#include <FL/Fl_File_Icon.H>
#include <FL/Fl_Shared_Image.H>

//////////////////////////////////////////////////////////////////////
// the one and only application object
CMainApp theapp;
//////////////////////////////////////////////////////////////////////

CMainApp* FLX_GetMainApp() { return &theapp; }
CWnd* FLX_GetMainWnd() { return theapp.GetMainWnd(); }

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMainApp::CMainApp()
{
	m_pMainWnd = 0;
	m_pDoc     = 0;
}

CMainApp::~CMainApp()
{
	delete m_pDoc;
}

bool CMainApp::Init()
{
	// setup FLTK
	fl_message_font(FL_HELVETICA, FL_NORMAL_SIZE);
	fl_register_images();
	Fl_File_Icon::load_system_icons();
	
	// Create the document
	m_pDoc = new CDocument;
	if (m_pDoc == 0) return false;

	// Create the main window and views
	m_pMainWnd = new CWnd(Fl::w() - 128, Fl::h() - 128, "FEBio - Task Manager", m_pDoc);
	if (m_pMainWnd == 0) return false;

	Fl::visual(FL_DOUBLE|FL_INDEX);

	m_pMainWnd->show();

	return true;
}

bool CMainApp::Load(const char* szfilename)
{
	// load the initial project
	if (strlen(szfilename)>0)
	{
//		m_pMainWnd->OpenFile(szfilename);
	}

	return true;
}

int CMainApp::Run() 
{
	// message loop is handled by FLTK
	return Fl::run();
}
