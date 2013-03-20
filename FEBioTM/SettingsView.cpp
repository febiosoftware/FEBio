#include "stdafx.h"
#include "SettingsView.h"
#include <FL/Fl_Check_Button.H>
#include <FL/fl_draw.H>
#include "Wnd.h"

//-----------------------------------------------------------------------------
CSettingsView::CSettingsView(CWnd* pwnd, int X, int Y, int W, int H, const char* sz) : Flx_Group(X, Y, W, H, sz)
{
	m_pWnd = pwnd;
	m_bdebug = false;

	Fl_Check_Button* pc;
	begin();
	{
		pc = new Fl_Check_Button(X+10, Y+10, 200, 20, "Debug mode"); DDX_BOOL(pc, &m_bdebug);
		AddCallback(pc, (FLX_CALLBACK) &CSettingsView::OnChange);
	}
	end();
	box(FL_DOWN_BOX);
	fl_color(FL_GRAY);
}

//-----------------------------------------------------------------------------
void CSettingsView::OnChange(Fl_Widget* pw, void* pd)
{
	UpdateData();

	CTask* pt = m_pWnd->GetSelectedTask();
	if (pt)
	{
		pt->m_bdebug = m_bdebug;
	}
}

//-----------------------------------------------------------------------------
void CSettingsView::Update()
{
	CTask* pt = m_pWnd->GetSelectedTask();
	if (pt)
	{
		m_bdebug = pt->m_bdebug;
	}
	UpdateData(false);
}
