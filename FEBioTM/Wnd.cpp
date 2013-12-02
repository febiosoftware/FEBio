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
#include <assert.h>

#ifdef WIN32
#define IDI_ICON1	101
extern HINSTANCE fl_display;
#endif

//-----------------------------------------------------------------------------
// style table entries
static Fl_Text_Display::Style_Table_Entry stable[] = {
       // FONT COLOR      FONT FACE   FONT SIZE
       // --------------- ----------- --------------
       {  FL_BLACK     , FL_COURIER, 12 }, // A - default
       {  FL_BLUE      , FL_COURIER, 12 }, // B - keyword
       {  FL_DARK_RED  , FL_COURIER, 12 }, // C - text (attribute value)
	   {  FL_MAGENTA   , FL_COURIER, 12 }, // D - attribute
	   {  FL_DARK_GREEN, FL_COURIER, 12 }, // E - comment
   };

//-----------------------------------------------------------------------------
// style table entries for test window
static Fl_Text_Display::Style_Table_Entry stable_test[] = {
       // FONT COLOR      FONT FACE   FONT SIZE
       // --------------- ----------- --------------
       {  FL_BLACK     , FL_COURIER, 12 }, // A - default
       {  FL_RED       , FL_COURIER, 12 }, // B - keyword
   };


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

	m_pTabs = 0;
	m_pText = 0;
	m_pOut = 0;
	m_pLog = 0;
	m_pTest = 0;

	m_nTestFmt = 0;
	m_pTestAtt = 0;

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
					Fl_Group* pg = new Fl_Group(wf, hm+ht+24, w-wf, h-hm-ht-24, "    Input   ");
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

					pg = m_pOps = new CSettingsView(this, wf, hm+ht+24, w-wf, h-hm-ht-24, "  Settings  ");
					pg->labelsize(11);

					pg = new Fl_Group(wf, hm+ht+24, w-wf, h-hm-ht-24, "   Output   ");
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

					pg = new Fl_Group(wf, hm+ht+24, w-wf, h-hm-ht-24, "    Log     ");
					{
						m_pLog = new Fl_Text_Display(wf, hm+ht+24, w-wf, h-hm-ht-24);
						m_pLog->textfont(FL_COURIER);
						m_pLog->box(FL_DOWN_BOX);
						m_pLog->buffer(new Fl_Text_Buffer);
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

	// set the windows callback
	AddCallback(this, (FLX_CALLBACK) &CWnd::OnFileExit);

#ifdef WIN32
	// set the icon of the window
	icon((char*) LoadIcon(fl_display, MAKEINTRESOURCE(IDI_ICON1)));
	show();
#endif

	m_pSel = m_pText;
	m_pTabs->do_callback();

	ClearLogWnd();
}

//-----------------------------------------------------------------------------
CWnd::~CWnd()
{

}

//-----------------------------------------------------------------------------
// clear the output window
void CWnd::ClearOutputWnd()
{
	Fl_Text_Buffer* pb = m_pOut->buffer();
	pb->select(0, pb->length());
	pb->remove_selection();	
}

//-----------------------------------------------------------------------------
void CWnd::ClearLogWnd()
{
	Fl_Text_Buffer* pb = m_pLog->buffer();
	pb->select(0, pb->length());
	pb->remove_selection();	

	m_pLog->insert(" FEBio Task Manager\n");
	m_pLog->insert("----------------------------------------------\n");
}

//-----------------------------------------------------------------------------
void CWnd::ClearTestWnd()
{
	Fl_Text_Buffer* pb = m_pTest->buffer();
	pb->select(0, pb->length());
	pb->remove_selection();	

	m_pTestAtt->select(0, m_pTestAtt->length());
	m_pTestAtt->remove_selection();
}

//-----------------------------------------------------------------------------
void CWnd::AddLogEntry(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	static char sztxt[1024] = {0};
	va_start(args, sz);
	vsprintf(sztxt, sz, args);
	va_end(args);

	// add the message to the log
	m_pLog->insert_position(m_pLog->buffer()->length());
	m_pLog->insert(sztxt);
	m_pLog->show_insert_position();	
}

//-----------------------------------------------------------------------------
void CWnd::AddTestEntry(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	static char sztxt[1024] = {0};
	va_start(args, sz);
	vsprintf(sztxt, sz, args);
	va_end(args);

	// add the message to the log
	m_pTest->insert_position(m_pTest->buffer()->length());
	m_pTest->insert(sztxt);
	m_pTest->show_insert_position();

	static char szstyle[2] = {'A', 'B'};
	int l = strlen(sztxt);
	for (int i=0; i<l; ++i) sztxt[i] = szstyle[m_nTestFmt];
	m_pTestAtt->insert(m_pTestAtt->length(), sztxt);

	m_pTest->highlight_data(m_pTestAtt, stable_test, sizeof(stable_test)/sizeof(stable[0]), 'A', 0, 0);
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
	CTask* pt = m_pDoc->GetSession().AddTask(szfile);
	if (pt)
	{
		m_pTask->AddTask(pt);
		if (m_pText)
		{
			m_pText->buffer(pt->GetTextBuffer());
			m_pText->highlight_data(pt->GetStyleBuffer(), stable, sizeof(stable)/sizeof(stable[0]), 'A', 0, 0);
		}
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
	if (m_pSel == 0) return;
	if (m_pSel == m_pText)
	{
		int n = m_pTask->SelectedTask();
		CTask* pt = m_pDoc->GetSession().GetTask(n);
		if (pt == 0) flx_error("No task selected");
		else 
		{
			pt->Save();
			pt->SetStatus(CTask::READY);
			m_pTask->redraw();
		}
	}
	else OnFileSaveAs(pw, pd);
}

//-----------------------------------------------------------------------------
void CWnd::OnFileSaveAs(Fl_Widget* pw, void* pd)
{
	if (m_pSel == 0) return;

	char szfile[1024] = {0};
	if (m_pSel == m_pText)
	{
		int n = m_pTask->SelectedTask();
		CTask* pt = m_pDoc->GetSession().GetTask(n);
		if (pt == 0) flx_error("No task selected");
		else 
		{
			if (flx_file_save(szfile, "FEBio files (*.feb)\t*.feb") == FLX_OK)
			{
				pt->Save(szfile);
				pt->SetStatus(CTask::READY);
				m_pTask->redraw();
			}
		}
	}
	else
	{
		if (flx_file_save(szfile, "Text files (*.txt)\t*.txt") == FLX_OK)
		{
			Fl_Text_Buffer* pb = m_pSel->buffer();
			if (pb->savefile(szfile) != 0) flx_error("Failed saving text buffer to file.");
		}
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnFileRevert(Fl_Widget* pw, void* pd)
{
	int n = m_pTask->SelectedTask();
	CTask* pt = m_pDoc->GetSession().GetTask(n);
	if (pt == 0) flx_error("No task selected");
	else 
	{
		if (pt->GetStatus() == CTask::MODIFIED)
		{
			if (flx_choice("This file has been modified. Are you sure you want to revert it?", "yes", "no", 0) == 1) return;
		}

		pt->Revert();
		m_pTask->redraw();
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnFileClose(Fl_Widget* pw, void* pd)
{
	int n = m_pTask->SelectedTask();
	if (n>=0)
	{
		if (m_pText) m_pText->buffer(0);
		m_pDoc->GetSession().RemoveTask(n);
		m_pTask->RemoveTask(n);
		SelectFile();
	}
	else flx_alert("Nothing to remove.");
}

//-----------------------------------------------------------------------------
void CWnd::OnFileCloseAll(Fl_Widget* pw, void* pd)
{
	m_pText->buffer(0);
	m_pDoc->NewSession();
	m_pTask->Update();
}

//-----------------------------------------------------------------------------
void CWnd::OnFileOpenSession(Fl_Widget* pw, void* pd)
{
	char szfile[1024] = {0};
	if (flx_file_open(szfile, "Session files\t*.ftm") == FLX_OK)
	{
		// close the current session
		OnFileCloseAll(0,0);
		ClearLogWnd();

		// open the new session
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
	CTask* pt = GetDocument()->GetSession().GetTask(m_pTask->SelectedTask());
	if (pt)
	{
		if (m_pText) 
		{
			m_pText->buffer(pt->GetTextBuffer());
			m_pText->highlight_data(pt->GetStyleBuffer(), stable, sizeof(stable)/sizeof(stable[0]), 'A', 0, 0);
		}
		if (m_pOps ) m_pOps->Update();
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnEditSelectAll(Fl_Widget* pw, void* pd)
{
	m_pTask->SelectAll();
}

//-----------------------------------------------------------------------------
void CWnd::OnEditFind(Fl_Widget* pw, void* pd)
{
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
	int npos = m_pSel->insert_position();
	Fl_Text_Buffer* pbuf = m_pSel->buffer();
	int found = pbuf->search_forward(npos, m_szfind, &npos, (m_bcase?1:0));
	if (found)
	{
		// search from the current position
		pbuf->select(npos, npos+strlen(m_szfind));
		m_pSel->insert_position(npos + strlen(m_szfind));
		m_pSel->show_insert_position();
		m_pSel->take_focus();
	}
	else 
	{
		// wrap around to the beginning
		int found = pbuf->search_forward(0, m_szfind, &npos, (m_bcase?1:0));
		if (found)
		{
			pbuf->select(npos, npos+strlen(m_szfind));
			m_pSel->insert_position(npos + strlen(m_szfind));
			m_pSel->show_insert_position();
			m_pSel->take_focus();
		}
		else flx_alert("Could not find string:\n\n%s", m_szfind);
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnEditFilter(Fl_Widget* pw, void* pd)
{
	CDlgEditFilter dlg;
	if (dlg.DoModal() == FLX_OK)
	{
		m_pTask->SelectAll(0);
		TMSession& session = m_pDoc->GetSession();
		int N = session.Tasks();
		int ncnt = 0;
		for (int i=0; i<N; ++i)
		{
			CTask& task = *session.GetTask(i);
			Fl_Text_Buffer* pbuf = task.GetTextBuffer();
			int npos;
			if (pbuf->search_forward(0, dlg.m_sztxt, &npos, dlg.m_bcase) > 0) { m_pTask->SelectTask(i); ncnt++; }
		}
		flx_alert("%d files found.", ncnt);
	}
}

//-----------------------------------------------------------------------------
void CWnd::OnEditGoToLine(Fl_Widget* pw, void* pd)
{
	CDlgEditGoToLine dlg;
	if (dlg.DoModal() == FLX_OK)
	{
		int nline = dlg.m_nline - 1;
		if (nline < 0) nline = 0;

		Fl_Text_Buffer* pbuf = m_pSel->buffer();
		int npos = pbuf->skip_lines(0, nline);
		m_pSel->insert_position(npos);
		m_pSel->show_insert_position();
		m_pSel->take_focus();
	}
}

//-----------------------------------------------------------------------------
CTask* CWnd::GetSelectedTask()
{
	return m_pDoc->GetSession().GetTask(m_pTask->SelectedTask());
}

//-----------------------------------------------------------------------------
void CWnd::OnRunSelected(Fl_Widget *pw, void *pd)
{
	// get the selected task
	int N = m_pTask->SelectedTasks();
	if (N == 0) { flx_error("No task selected"); return; }

	// show the output window
	m_pTabs->value(m_pTabs->child(2));
	m_pTabs->do_callback();

	TMSession& session = m_pDoc->GetSession();
	N = session.Tasks();
	for (int i=0; i<N; ++i)
	{
		if (m_pTask->IsTaskSelected(i))
		{
			CTask* pt = session.GetTask(i);

			// Add the task to the queue
			m_pDoc->RunTask(pt);
		}
	}

	m_pTask->redraw();
}

//-----------------------------------------------------------------------------
void CWnd::OnRunSession(Fl_Widget* pw, void* pd)
{
	// show the output window
	if (m_pTabs) m_pTabs->value(m_pTabs->child(2));
	m_pTask->redraw();

	// run the session
	m_pDoc->RunSession();
}

//-----------------------------------------------------------------------------
// Stop the task that the user has selected
void CWnd::OnRunCancelSelected(Fl_Widget* pw, void* pd)
{
	TMSession& session = m_pDoc->GetSession();
	int N = session.Tasks();
	for (int i=0; i<N; ++i)
	{
		if (m_pTask->IsTaskSelected(i))
		{
			CTask* pt = session.GetTask(i);
			if (pt)
			{
				int n = pt->GetStatus();
				if ((n==CTask::RUNNING)||(n==CTask::QUEUED)) pt->SetStatus(CTask::CANCELLED);
			}
		}
	}
	m_pTask->redraw();
}

//-----------------------------------------------------------------------------
// Stop all queued tasks
void CWnd::OnRunCancelAll(Fl_Widget* pw, void* pd)
{
	for (int i=0; i<m_pDoc->GetSession().Tasks(); ++i)
	{
		CTask* pt = m_pDoc->GetSession().GetTask(i);
		int n = pt->GetStatus();
		if ((n==CTask::RUNNING)||(n==CTask::QUEUED)) pt->SetStatus(CTask::CANCELLED);
	}
}

//-----------------------------------------------------------------------------
// Create a new test
void CWnd::OnToolsCreateTest(Fl_Widget* pw, void* pd)
{
	int X = m_pTabs->x();
	int Y = m_pTabs->y();
	int W = m_pTabs->w();
	int H = m_pTabs->h();

	if (m_pTest == 0)
	{
		Fl_Group* pg;
		m_pTabs->add(pg = new Fl_Group(X, Y+24, W, H-24, "    Test     "));
		{
			m_pTest = new Fl_Text_Display(X, Y+24, W, H-24);
			m_pTest->textfont(FL_COURIER);
			m_pTest->box(FL_DOWN_BOX);
			m_pTest->buffer(new Fl_Text_Buffer);
			pg->resizable(m_pTest);
			m_pTestAtt = new Fl_Text_Buffer;
		}
		pg->end();
		pg->labelsize(11);
		OnSelectTab(0, 0);
		m_pTabs->redraw();
	}

	CTest* ptest = new CTest(m_pDoc->GetSession());
	m_pDoc->AddTest(ptest);
}

//-----------------------------------------------------------------------------
// run a test
void CWnd::OnToolsRunTest(Fl_Widget* pw, void* pd)
{
	ClearTestWnd();
	m_pDoc->RunTest();
}

//-----------------------------------------------------------------------------
void CWnd::OnToolsLoadData(Fl_Widget* pw, void* pd)
{
	CTest* ptest = m_pDoc->GetActiveTest();
	if (ptest)
	{
		char szfilename[512] = {0};
		char szfilter[] = "All files\t*\n";
		if (flx_file_open(szfilename, szfilter) == FLX_OK)
		{
			bool bret = ptest->LoadData(szfilename);
			if (bret == false) flx_error("Failed reading meta data.");
		}
	}
	else flx_error("No active test.");
}

//-----------------------------------------------------------------------------
void CWnd::OnToolsSaveData(Fl_Widget* pw, void* pd)
{
	CTest* ptest = m_pDoc->GetActiveTest();
	if (ptest)
	{
		char szfilename[512] = {0};
		char szfilter[] = "All files\t*\n";
		if (flx_file_save(szfilename, szfilter) == FLX_OK)
		{
			bool bret = ptest->SaveData(szfilename);
			if (bret == false) flx_error("Failed saving meta data.");
		}
	}
	else flx_error("No active test.");
}

//-----------------------------------------------------------------------------
void CWnd::OnSelectTab(Fl_Widget* pw, void* pd)
{
	Fl_Group* ps = dynamic_cast<Fl_Group*>(m_pTabs->value()); assert(ps);
	int n = m_pTabs->children();
	for (int i=0; i<n; ++i) 
	{
		Fl_Widget* pc = m_pTabs->child(i);
		pc->labelfont(FL_HELVETICA);
		pc->selection_color(FL_DARK2);
	}
	ps->labelfont(FL_HELVETICA_BOLD);
	ps->selection_color(FL_GRAY);
	Fl_Widget* pc = ps->child(0);
	m_pSel = 0;
	if (pc == m_pText) m_pSel = m_pText;
	if (pc == m_pOut ) m_pSel = m_pOut;
	if (pc == m_pLog ) m_pSel = m_pLog;
	if (pc == m_pTest) m_pSel = m_pTest;
}

//-----------------------------------------------------------------------------
void CWnd::OnChangeText(Fl_Widget* pw, void* pd)
{
	int n = m_pTask->SelectedTask();
	CTask* pt = m_pDoc->GetSession().GetTask(n);

	// we need to update the format of the modified line
	int npos = m_pText->insert_position();
	int nstart = m_pText->line_start(npos);
	int nend = m_pText->line_end(nstart, true);

//	pt->format_style(nstart, nend);

	if (pt && (pt->GetStatus() != CTask::MODIFIED))
	{
		pt->SetStatus(CTask::MODIFIED);
		m_pTask->redraw();
	}
}
