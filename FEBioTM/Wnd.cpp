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
#include <FL/Fl_Preferences.H>
#include "MainApp.h"
#include <flx_message.h>
#include "DlgEditFind.h"

#ifdef WIN32
#define IDI_ICON1	101
extern HINSTANCE fl_display;
#endif

//-----------------------------------------------------------------------------
const char* wnd_title = "FEBio Task Manager";

//-----------------------------------------------------------------------------
CWnd::CWnd(int w, int h, const char* sztitle, CDocument* pdoc) : Flx_Wnd(w, h, wnd_title), m_pDoc(pdoc)
{
	int hm = 27;	// menu height
	int ht = 200;	// task browser height
	int wf = 400;	// file browser width

	// set the custom user interface settings
	// normal theme
	Fl::background(236, 233, 216);
	Fl::set_color((Fl_Color)1, 255, 200, 200);

	m_szfind[0] = 0;
	m_bcase = false;

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
						m_pText = new Fl_Text_Editor(wf, hm+ht+24, w-wf, h-hm-ht-24);
						m_pText->textfont(FL_COURIER);
						m_pText->box(FL_DOWN_BOX);
						pg->resizable(m_pText);
						AddCallback(m_pText, (FLX_CALLBACK) &CWnd::OnChangeText);
						m_pText->when(FL_WHEN_CHANGED);
					}
					pg->end();
					m_pTabs->resizable(pg);
					pg->labelsize(11);

					pg = new Fl_Group(wf, hm+ht+24, w-wf, h-hm-ht-24, "    Output    ");
					{
						m_pOut = new Fl_Text_Display(wf, hm+ht+24, w-wf, h-hm-ht-24);
						m_pOut->textfont(FL_COURIER);
						m_pOut->box(FL_DOWN_BOX);
						m_pOut->color(FL_BLACK);
						m_pOut->textcolor(FL_WHITE);
						m_pOut->buffer(new Fl_Text_Buffer);
						pg->resizable(m_pOut);
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

	// set the windows callback
	AddCallback(this, (FLX_CALLBACK) &CWnd::OnFileExit);

#ifdef WIN32
	// set the icon of the window
	icon((char*) LoadIcon(fl_display, MAKEINTRESOURCE(IDI_ICON1)));
	show();
#endif

	m_pTabs->do_callback();
}

//-----------------------------------------------------------------------------
CWnd::~CWnd()
{

}

//-----------------------------------------------------------------------------
// clear the output window
void CWnd::ClearOutputWnd()
{
	Fl_Text_Buffer* plog = m_pOut->buffer();
	plog->select(0, plog->length());
	plog->remove_selection();	
}

//-----------------------------------------------------------------------------
void CWnd::Update()
{
	
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
void CWnd::OnFileSave(Fl_Widget* pw, void* pd)
{
	int n = m_pTask->SelectedTask();
	CTask* pt = m_pDoc->GetTask(n);
	if (pt == 0) flx_error("No task selected");
	else 
	{
		pt->Save();
		pt->SetStatus(CTask::READY);
		m_pTask->redraw();
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnFileSaveAs(Fl_Widget* pw, void* pd)
{
	int n = m_pTask->SelectedTask();
	CTask* pt = m_pDoc->GetTask(n);
	if (pt == 0) flx_error("No task selected");
	else 
	{
		char szfile[1024] = {0};
		if (flx_file_save(szfile, "FEBio files (*.feb)\t*.feb") == FLX_OK)
		{
			pt->Save(szfile);
			pt->SetStatus(CTask::READY);
			m_pTask->redraw();
		}
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnFileClose(Fl_Widget* pw, void* pd)
{
	int n = m_pTask->SelectedTask();
	if (n>=0)
	{
		m_pText->buffer(0);
		m_pDoc->RemoveTask(n);
		m_pTask->RemoveTask(n);
		SelectFile();
	}
	else flx_alert("Nothing to remove.");
}

//-----------------------------------------------------------------------------
void CWnd::OnFileCloseAll(Fl_Widget* pw, void* pd)
{
	m_pText->buffer(0);
	m_pOut->buffer(0);
	m_pDoc->NewSession();
	m_pTask->Update();
}

//-----------------------------------------------------------------------------
void CWnd::OnFileExit(Fl_Widget *pw, void *pd)
{
	// save current working directory to preferences
	Fl_Preferences& pref = FLXGetMainApp()->GetPreferences();
	pref.set("cwd", m_pFile->GetCWD());

	// close application
	hide();
}

//-----------------------------------------------------------------------------
void CWnd::OnSelectFile(Fl_Widget* pw, void* pd)
{
	CTaskTable* pt = dynamic_cast<CTaskTable*>(pw);
	assert(pt);
	if (pt->callback_context() == Fl_Table::CONTEXT_CELL) SelectFile();
}

//-----------------------------------------------------------------------------
void CWnd::SelectFile()
{
	CTask* pt = GetDocument()->GetTask(m_pTask->SelectedTask());
	if (pt)
	{
		m_pText->buffer(pt->GetTextBuffer());
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnEditFind(Fl_Widget* pw, void* pd)
{
	// TODO: find which text buffer is active
	//       for now, let's assume the input buffer
	CDlgEditFind dlg;
	strcpy(dlg.m_sztxt, m_szfind);
	dlg.m_bcase = m_bcase;
	if (dlg.DoModal() == FLX_OK)
	{
		strcpy(m_szfind, dlg.m_sztxt);
		m_bcase = dlg.m_bcase;
		if (m_szfind[0] != 0) OnEditFindAgain(pw, pd);
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnEditFindAgain(Fl_Widget* pw, void* pd)
{
	int npos = m_pText->insert_position();
	Fl_Text_Buffer* pbuf = m_pText->buffer();
	int found = pbuf->search_forward(npos, m_szfind, &npos, (m_bcase?1:0));
	if (found)
	{
		pbuf->select(npos, npos+strlen(m_szfind));
		m_pText->insert_position(npos + strlen(m_szfind));
		m_pText->show_insert_position();
	}
	else flx_alert("Could not find string:\n\n%s", m_szfind);
}

//-----------------------------------------------------------------------------
void CWnd::OnRunSelected(Fl_Widget *pw, void *pd)
{
	CTask* pt = m_pDoc->GetTask(m_pTask->SelectedTask());
	if (pt == 0) { flx_error("No task selected"); return; }

	m_pTabs->value(m_pTabs->child(1));
	m_pTabs->do_callback();
	Fl::check();
	m_pDoc->RunTask(pt);
	m_pTask->redraw();
	label(wnd_title);
}

//-----------------------------------------------------------------------------
void CWnd::OnRunSession(Fl_Widget* pw, void* pd)
{
	// Add all tasks to the queue
	for (int i=0; i<m_pDoc->Tasks(); ++i)
	{
		CTask* pt = m_pDoc->GetTask(i);
		if (pt->GetStatus() == CTask::MODIFIED) pt->Save();
		pt->SetStatus(CTask::QUEUED);
	}
	m_pTabs->value(m_pTabs->child(1));
	m_pTask->redraw();
	Fl::flush();

	// run all tasks
	for (int i=0; i<m_pDoc->Tasks(); ++i) 
	{
		CTask* pt = m_pDoc->GetTask(i);
		m_pTask->SelectTask(i);
		SelectFile();
		m_pDoc->RunTask(pt);
		Fl::flush();
	}

	label(wnd_title);
}

//-----------------------------------------------------------------------------
void CWnd::OnRunStop(Fl_Widget* pw, void* pd)
{
	int n = m_pTask->SelectedTask();
	CTask* pt = m_pDoc->GetTask(n);
	if (pt && (pt->GetStatus() == CTask::RUNNING))
	{
		pt->SetStatus(CTask::CANCELLED);
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

//-----------------------------------------------------------------------------
void CWnd::OnChangeText(Fl_Widget* pw, void* pd)
{
	int n = m_pTask->SelectedTask();
	CTask* pt = m_pDoc->GetTask(n);
	if (pt && (pt->GetStatus() != CTask::MODIFIED))
	{
		pt->SetStatus(CTask::MODIFIED);
		m_pTask->redraw();
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnFileOpenSession(Fl_Widget* pw, void* pd)
{
	char szfile[1024] = {0};
	if (flx_file_open(szfile, "Session files\t*.ftm") == FLX_OK)
	{
		if (m_pDoc->OpenSession(szfile) == false) flx_error("Failed opening session");
		m_pTask->Update();
		SelectFile();
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnFileSaveSession(Fl_Widget* pw, void* pd)
{
	char szfile[1024] = {0};
	if (flx_file_save(szfile, "Session files\t*.ftm") == FLX_OK)
	{
		if (m_pDoc->SaveSession(szfile) == false) flx_error("Failed saving session");
	}
}
