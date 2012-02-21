// MainApp.cpp: implementation of the CMainApp class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MainApp.h"
#include <FL/fl_ask.H>
#include <FL/Fl_File_Icon.H>
#include <FL/Fl_Shared_Image.H>

  static short	plain[] = {	// Plain file icon
		Fl_File_Icon::COLOR, 0, FL_WHITE, Fl_File_Icon::POLYGON,
		Fl_File_Icon::VERTEX, 2000, 1000, Fl_File_Icon::VERTEX, 2000, 9000,
		Fl_File_Icon::VERTEX, 6000, 9000, Fl_File_Icon::VERTEX, 8000, 7000,
		Fl_File_Icon::VERTEX, 8000, 1000, Fl_File_Icon::END,
		Fl_File_Icon::COLOR, 0, FL_DARK3, Fl_File_Icon::LINE, 
		Fl_File_Icon::VERTEX, 2000, 1000, Fl_File_Icon::VERTEX, 2000, 9000, Fl_File_Icon::END,
		Fl_File_Icon::COLOR, 0, FL_DARK2, Fl_File_Icon::LINE, 
		Fl_File_Icon::VERTEX, 2000, 9000, Fl_File_Icon::VERTEX, 6000, 9000, Fl_File_Icon::END,
		Fl_File_Icon::COLOR, 0, FL_BLACK, Fl_File_Icon::LINE, 
		Fl_File_Icon::VERTEX, 6000, 9000, Fl_File_Icon::VERTEX, 8000, 7000,
		Fl_File_Icon::VERTEX, 8000, 1000, Fl_File_Icon::VERTEX, 2000, 1000, Fl_File_Icon::END,
		Fl_File_Icon::COLOR, 0, FL_BLACK, Fl_File_Icon::POLYGON,
		Fl_File_Icon::VERTEX, 4000, 4000, Fl_File_Icon::VERTEX, 4000, 6000,
		Fl_File_Icon::VERTEX, 5000, 7000, Fl_File_Icon::VERTEX, 6000, 6000,
		Fl_File_Icon::VERTEX, 6000, 4000, Fl_File_Icon::END,
		Fl_File_Icon::COLOR, 0, FL_RED, Fl_File_Icon::POLYGON,
		Fl_File_Icon::VERTEX, 4000, 6000, Fl_File_Icon::VERTEX, 5000, 7000,
		Fl_File_Icon::VERTEX, 6000, 6000, Fl_File_Icon::VERTEX, 5000, 5000, Fl_File_Icon::END,
		Fl_File_Icon::END
		};


static Fl_Color c = fl_rgb_color(255, 235, 133);
static short c1 = (c>>16)&(0xFFFF);
static short c2 = c&(0xFFFF);

 static short	dir[] = {	// Directory icon
		  Fl_File_Icon::COLOR, c1, c2, Fl_File_Icon::POLYGON, Fl_File_Icon::VERTEX, 1000, 1000,
		  Fl_File_Icon::VERTEX, 1000, 7500,  Fl_File_Icon::VERTEX, 9000, 7500,
		  Fl_File_Icon::VERTEX, 9000, 1000, Fl_File_Icon::END,
		  Fl_File_Icon::POLYGON, Fl_File_Icon::VERTEX, 1000, 7500, Fl_File_Icon::VERTEX, 2500, 9000,
		  Fl_File_Icon::VERTEX, 5000, 9000, Fl_File_Icon::VERTEX, 6500, 7500, Fl_File_Icon::END,
		  Fl_File_Icon::COLOR, 0, FL_WHITE, Fl_File_Icon::LINE, Fl_File_Icon::VERTEX, 1500, 1500,
		  Fl_File_Icon::VERTEX, 1500, 7000, Fl_File_Icon::VERTEX, 9000, 7000, Fl_File_Icon::END,
		  Fl_File_Icon::COLOR, 0, FL_BLACK, Fl_File_Icon::LINE, Fl_File_Icon::VERTEX, 9000, 7500,
		  Fl_File_Icon::VERTEX, 9000, 1000, Fl_File_Icon::VERTEX, 1000, 1000, Fl_File_Icon::END,
		  Fl_File_Icon::COLOR, 0, FL_GRAY, Fl_File_Icon::LINE, Fl_File_Icon::VERTEX, 1000, 1000,
		  Fl_File_Icon::VERTEX, 1000, 7500, Fl_File_Icon::VERTEX, 2500, 9000,
		  Fl_File_Icon::VERTEX, 5000, 9000, Fl_File_Icon::VERTEX, 6500, 7500,
		  Fl_File_Icon::VERTEX, 9000, 7500, Fl_File_Icon::END,
		  Fl_File_Icon::END
		};


//////////////////////////////////////////////////////////////////////
// the one and only application object
CMainApp theapp;
//////////////////////////////////////////////////////////////////////

CMainApp* FLXGetMainApp() { return &theapp; }
CWnd* FLXGetMainWnd() { return theapp.GetMainWnd(); }

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMainApp::CMainApp() : m_prefs(Fl_Preferences::USER, "FEBio2", "FEBioTM")
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

	new Fl_File_Icon("*.feb", Fl_File_Icon::PLAIN, sizeof(plain) / sizeof(plain[0]), plain);
	new Fl_File_Icon("*", Fl_File_Icon::DIRECTORY, sizeof(dir) / sizeof(dir[0]), dir);

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
		m_pMainWnd->OpenFile(szfilename);
	}

	return true;
}

int CMainApp::Run() 
{
	// message loop is handled by FLTK
	return Fl::run();
}
