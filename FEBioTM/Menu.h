#pragma once
#include <FL/Fl_Menu_Bar.H>

class CWnd;

class CMenu : public Fl_Menu_Bar
{
public:
	CMenu(int x, int h, CWnd* pwnd);
	virtual ~CMenu() { m_pWnd = 0; }

protected:
	CWnd*	m_pWnd;
};
