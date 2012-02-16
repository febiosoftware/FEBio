// Wnd.cpp: implementation of the CWnd class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Wnd.h"
#include <FL/Fl.H>

//-----------------------------------------------------------------------------
CWnd::CWnd(int w, int h, const char* sztitle, CDocument* pdoc) : Flx_Wnd(w, h)
{
	m_pDoc = pdoc;

	// set the custom user interface settings
	// normal theme
	Fl::background(236, 233, 216);

	// Graphite theme
//	Fl::background(100, 100, 100);
//	Fl::foreground(255,255,255);
	Fl::set_color((Fl_Color)1, 255, 200, 200);

	begin();
	{
	}
	end();

//	box(FL_NO_BOX); // no background filling

	size_range(400, 300);

}

//-----------------------------------------------------------------------------
CWnd::~CWnd()
{

}

//-----------------------------------------------------------------------------
bool CWnd::OpenFile(const char* szfile)
{
	return m_pDoc->OpenFile(szfile);
}
