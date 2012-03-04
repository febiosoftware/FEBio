#include "stdafx.h"
#include "Menu.h"
#include "Wnd.h"

//-----------------------------------------------------------------------------
CMenu::CMenu(int w, int h, CWnd* pwnd) : Fl_Menu_Bar(0, 0, w, h), m_pWnd(pwnd)
{
	// construct the menu
	static Fl_Menu_Item m[] = {
		{"File", 0, 0, 0, FL_SUBMENU },
			{"&Open ..."  , FL_CTRL + 'o', FLX_MENU_HANDLER(pwnd, CWnd::OnFileOpen)},
			{"Save"       , FL_CTRL + 's', FLX_MENU_HANDLER(pwnd, CWnd::OnFileSave)},
			{"Save as ...", FL_CTRL + 'a', FLX_MENU_HANDLER(pwnd, CWnd::OnFileSaveAs), FL_MENU_DIVIDER},
			{"Close"     ,       0, FLX_MENU_HANDLER(pwnd, CWnd::OnFileClose)},
			{"Close all" ,       0, FLX_MENU_HANDLER(pwnd, CWnd::OnFileCloseAll), FL_MENU_DIVIDER},
			{"Open Session ...", 0, FLX_MENU_HANDLER(pwnd, CWnd::OnFileOpenSession)},
			{"Save Session ...", 0, FLX_MENU_HANDLER(pwnd, CWnd::OnFileSaveSession), FL_MENU_DIVIDER},
			{"&Exit"      ,             0, FLX_MENU_HANDLER(pwnd, CWnd::OnFileExit)},
			{0},
		{"Edit", 0, 0, 0, FL_SUBMENU},
			{"Find ...", FL_CTRL + 'f', FLX_MENU_HANDLER(pwnd, CWnd::OnEditFind)},
			{0},
		{"Run", 0, 0, 0, FL_SUBMENU },
			{"Run Selected", FL_CTRL + FL_F + 5, FLX_MENU_HANDLER(pwnd, CWnd::OnRunSelected)},
			{"Run Session" ,           FL_F + 5, FLX_MENU_HANDLER(pwnd, CWnd::OnRunSession)},
			{"Stop"        ,                  0, FLX_MENU_HANDLER(pwnd, CWnd::OnRunStop)},
			{0},
		{0}
	};

	// set the textsize
	textsize(11);

	// set the menu items
	menu(m);
}
