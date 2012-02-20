// Wnd.cpp: implementation of the CWnd class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Wnd.h"
#include "Document.h"
#include <FL/Fl.H>
#include <FL/Fl_Tile.H>
#include <Flex.h>
#include <Flx_Dialog.h>
#include <flx_message.h>

//-----------------------------------------------------------------------------
CWnd::CWnd(int w, int h, const char* sztitle, CDocument* pdoc) : Flx_Wnd(w, h, "FEBio Task Manager"), m_pDoc(pdoc)
{
	int hm = 27;	// menu height
	int ht = 200;	// task browser height
	int wf = 400;	// file browser width

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

			m_pTask = new CTaskBrowser(wf, hm, w-wf, ht, this);

			pg = new Fl_Group(wf, hm+ht, w-wf, h-hm-ht);
			{
				m_pTabs = new Fl_Tabs(wf, hm+ht, w-wf, h-hm-ht);
				{
					Fl_Group* pg = new Fl_Group(wf, hm+ht+24, w-wf, h-hm-ht-24, "    Input    ");
					{
						m_pText = new Fl_Text_Display(wf, hm+ht+24, w-wf, h-hm-ht-24);
						m_pText->textfont(FL_COURIER);
						m_pText->box(FL_DOWN_BOX);
						pg->resizable(m_pText);
					}
					pg->end();
					m_pTabs->resizable(pg);
					pg->labelsize(11);

					pg = new Fl_Group(wf, hm+ht+24, w-wf, h-hm-ht-24, "    Log    ");
					{
						m_pLog = new Fl_Text_Display(wf, hm+ht+24, w-wf, h-hm-ht-24);
						m_pLog->textfont(FL_COURIER);
						m_pLog->box(FL_DOWN_BOX);
						m_pLog->color(FL_BLACK);
						m_pLog->textcolor(FL_WHITE);
						pg->resizable(m_pLog);
					}
					pg->end();
					pg->labelsize(11);
				}
				m_pTabs->end();
				pg->resizable(m_pTabs);
				AddCallback(m_pTabs, (FLX_CALLBACK) &CWnd::OnSelectTab);
			}
			pg->end();
			pg->box(FL_FLAT_BOX);
			pg->color(FL_DARK2);
		}
		pt->end();
		resizable(pt);
	}
	end();

	box(FL_FLAT_BOX); // no background filling
	color(FL_DARK3);
	size_range(400, 300);

	m_pTabs->do_callback();
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
		m_pText->buffer(pt->GetTextBuffer());
		m_pLog->buffer(pt->GetLogBuffer());
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

//-----------------------------------------------------------------------------
void CWnd::OnSelectFile(Fl_Widget* pw, void* pd)
{
	CTaskTable* pt = dynamic_cast<CTaskTable*>(pw);
	assert(pt);
	if (pt->callback_context() == Fl_Table::CONTEXT_CELL)
	{
		int n = m_pTask->SelectedTask();
		assert(n >= 0);
		CTask* pt = GetDocument()->GetTask(n);
		m_pText->buffer(pt->GetTextBuffer());
		m_pLog->buffer(pt->GetLogBuffer());
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnRunSelected(Fl_Widget *pw, void *pd)
{
	int n = m_pTask->SelectedTask();
	if ((n < 0) || (n >= m_pDoc->Tasks())) flx_error("No task selected");
	else 
	{
		m_pTabs->value(m_pTabs->child(1));
		m_pTabs->do_callback();
		Fl::flush();
		m_pDoc->RunTask(n);
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnSelectTab(Fl_Widget* pw, void* pd)
{
	Fl_Widget* ps = m_pTabs->value();
	int n = m_pTabs->children();
	for (int i=0; i<n; ++i) 
	{
		Fl_Widget* pc = m_pTabs->child(i);
		pc->labelfont(FL_HELVETICA);
		pc->selection_color(FL_DARK2);
	}
	ps->labelfont(FL_HELVETICA_BOLD);
	ps->selection_color(FL_GRAY);
}
