// Wnd.cpp: implementation of the CWnd class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Wnd.h"
#include "Document.h"
#include "TaskBar.h"
#include <FL/Fl.H>
#include <FL/Fl_Tile.H>
#include <Flex.h>
#include <Flx_Dialog.h>
#include <FL/Fl_Tabs.H>

//-----------------------------------------------------------------------------
CWnd::CWnd(int w, int h, const char* sztitle, CDocument* pdoc) : Flx_Wnd(w, h, "FEBio Task Manager")
{
	int hm = 27;	// menu height
	int wf = 400;	// file browser width

	// store pointer to document
	m_pDoc = pdoc;

	// set the custom user interface settings
	// normal theme
	Fl::background(236, 233, 216);
	Fl::set_color((Fl_Color)1, 255, 200, 200);

	Fl_Group* pg;
	begin();
	{
		// add the menu
		m_pMenu = new CMenu(w, hm, this);

		Fl_Tile* pt = new Fl_Tile(0, hm, w, h-hm);
		{
			m_pFile = new CFileBrowser(0, hm, wf, h-hm, this);

			m_pTask = new CTaskBrowser(wf, hm, w-wf, 400, this);

			Fl_Tabs* ptabs = new Fl_Tabs(wf, hm+400, w-wf, h-hm-400);
			{
				m_pText = new Fl_Text_Display(wf, hm+400+25, w-wf, h-hm-400-25, "    Input    ");
				m_pText->textfont(FL_COURIER);
				m_pText->buffer(new Fl_Text_Buffer);
				m_pText->box(FL_DOWN_BOX);

				pg = new Fl_Group(wf, hm+400+25, w-wf, h-hm-400-25, "    Log    ");
				{
				}
				pg->end();
			}
			ptabs->end();
//			ptabs->color(FL_DARK3);
		}
		pt->end();
		pt->resizable(m_pFile);

		resizable(pt);
	}
	end();

//	box(FL_NO_BOX); // no background filling
	color(FL_DARK3);

	size_range(400, 300);

}

//-----------------------------------------------------------------------------
CWnd::~CWnd()
{

}

//-----------------------------------------------------------------------------
void CWnd::Update()
{
	m_pTask->Update();
}

//-----------------------------------------------------------------------------
int CWnd::handle(int nevent)
{
	switch (nevent)
	{
	case FL_KEYDOWN:
		switch (Fl::event_key())
		{
		case FL_Escape:
			{
				// capture the escape key, since FLTK's default
				// behaviour of closing the window is not desired
				// for the main window.
				return 1;
			}
			break;
		}
		break;
	};
	
	return Flx_Wnd::handle(nevent);
}

//-----------------------------------------------------------------------------
bool CWnd::OpenFile(const char* szfile)
{
	CTask* pt = m_pDoc->AddTask(szfile);
	if (pt)
	{
		m_pTask->AddTask(pt);
		m_pText->buffer()->appendfile(szfile);
	}
	return (pt != 0);
}

//-----------------------------------------------------------------------------
void CWnd::OnFileOpen(Fl_Widget *pw, void *pd)
{
	char szfilename[512] = {0};
	char szfilter[] = "FEBio input files\t*.feb\n";
	if (flx_file_open(szfilename, szfilter) == FLX_OK)
	{
		OpenFile(szfilename);
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnFileExit(Fl_Widget *pw, void *pd)
{
	hide();
}
