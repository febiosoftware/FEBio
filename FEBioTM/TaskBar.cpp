#include "stdafx.h"
#include "TaskBar.h"
#include "Document.h"
#include <FL/Fl_Box.H>

CTaskBar::CTaskBar(int x, int y, int w, CTask* pt) : Flx_Group(x, y, w, 30), m_pTask(pt)
{
	int H = h();
	Fl_Box* ps;
	begin();
	{
		ps = new Fl_Box(x+50, y, w-50, H, pt->GetFileName());
		ps->align(FL_ALIGN_INSIDE | FL_ALIGN_LEFT);
	}
	end();

	box(FL_ROUND_UP_BOX);
	color(FL_LIGHT2);
}
