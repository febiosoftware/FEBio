#include "stdafx.h"
#include "Menu.h"
#include "Wnd.h"

//-----------------------------------------------------------------------------
CMenu::CMenu(int w, int h, CWnd* pwnd) : Fl_Menu_Bar(0, 0, w, h), m_pWnd(pwnd)
{
	// construct the menu
	static Fl_Menu_Item m[] = {
		{"File", 0, 0, 0, FL_SUBMENU },
			{"&Open ...", FL_CTRL + 'o', FLX_MENU_HANDLER(pwnd, CWnd::OnFileOpen)},
			{"&Exit"    ,             0, FLX_MENU_HANDLER(pwnd, CWnd::OnFileExit)},
			{0},
		{0}
	};

	// set the textsize
	textsize(11);

	// set the menu items
	menu(m);
}
