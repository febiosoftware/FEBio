#pragma once
#include <Flx_Group.h>
#include "TaskBar.h"
#include <FL/Fl_Pack.H>
#include <FL/Fl_Scroll.H>

class CWnd;

class CTaskBrowser : public Flx_Group
{
public:
	CTaskBrowser(int x, int y, int w, int h, CWnd* pwnd);

	void Update();

	void AddTask(CTask* pt);

protected:
	CWnd*	m_pWnd;
	Fl_Pack*	m_pg;
};
