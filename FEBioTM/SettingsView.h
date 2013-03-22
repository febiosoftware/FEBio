#pragma once
#include <Flx_Group.h>

//-----------------------------------------------------------------------------
// forward declarations
class CWnd;

//-----------------------------------------------------------------------------
class CSettingsView : public Flx_Group
{
public:
	CSettingsView(CWnd* pwnd, int X, int Y, int W, int H, const char* sz);

	void Update();

protected:
	void OnChange(Fl_Widget* pw, void* pd);

public:
	bool	m_bdebug;	//!< debug mode
	int		m_nlog;		//!< log level

private:
	CWnd*	m_pWnd;
};
