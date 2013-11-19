#include "stdafx.h"
#include "SettingsView.h"
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/fl_draw.H>
#include "Wnd.h"

Fl_Menu_Item log_menu[] = {
	{"Default"},
	{"Never"},
	{"Progress"},
	{"Major Iterations"},
	{"Minor Iterations"},
	{0}
};

//-----------------------------------------------------------------------------
CSettingsView::CSettingsView(CWnd* pwnd, int X, int Y, int W, int H, const char* sz) : Flx_Group(X, Y, W, H, sz)
{
	m_pWnd = pwnd;
	m_bdebug = false;
	m_nlog = 0;

	Fl_Check_Button* pc;
	Fl_Choice* pl;
	begin();
	{
		pc = new Fl_Check_Button(X+70, Y+10, 200, 20, "Debug mode"); DDX_BOOL(pc, &m_bdebug);
		AddCallback(pc, (FLX_CALLBACK) &CSettingsView::OnChange);

		pl = new Fl_Choice(X+70, Y+35, 200, 20, "Log level"); pl->copy(log_menu); DDX_CHOICE(pl, &m_nlog);
		AddCallback(pl, (FLX_CALLBACK) &CSettingsView::OnChange);
	}
	end();
	resizable(0);
	box(FL_DOWN_BOX);
	fl_color(FL_GRAY);
	UpdateData(false);
}

//-----------------------------------------------------------------------------
void CSettingsView::OnChange(Fl_Widget* pw, void* pd)
{
	UpdateData();

	TMSession& session = m_pWnd->GetDocument()->GetSession();

	int N = session.Tasks();
	for (int i=0; i<N; ++i)
	{
		if (m_pWnd->IsTaskSelected(i))
		{
			CTask* pt = session.GetTask(i);
			if (pt)
			{
				pt->m_bdebug = m_bdebug;
				pt->m_nlog   = m_nlog;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void CSettingsView::Update()
{
	CTask* pt = m_pWnd->GetSelectedTask();
	if (pt)
	{
		m_bdebug = pt->m_bdebug;
		m_nlog   = pt->m_nlog;
	}
	UpdateData(false);
}
