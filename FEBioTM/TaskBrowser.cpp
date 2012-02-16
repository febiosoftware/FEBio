#include "stdafx.h"
#include "TaskBrowser.h"
#include "Wnd.h"
#include "Document.h"

CTaskBrowser::CTaskBrowser(int x, int y, int w, int h, CWnd* pwnd) : Flx_Group(x, y, w, h), m_pWnd(pwnd)
{
	begin();
	{
		Fl_Scroll* ps = new Fl_Scroll(x, y, w, h);
		{
			Fl_Pack* pk = m_pg = new Fl_Pack(x, y, w, h);
			{
			}
			pk->end();
			pk->type(Fl_Pack::VERTICAL);
			pk->box(FL_NO_BOX);
//			pk->color(FL_DARK3);
			ps->resizable(pk);
		}
		ps->end();
		ps->box(FL_NO_BOX);
		ps->type(Fl_Scroll::VERTICAL);
	}
	end();
	box(FL_DOWN_BOX);
}

void CTaskBrowser::Update()
{
}

void CTaskBrowser::AddTask(CTask *pt)
{
	m_pg->add(new CTaskBar(x(), y(), w(), pt));
	redraw();
}
