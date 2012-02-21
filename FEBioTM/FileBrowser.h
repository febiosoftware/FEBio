#pragma once
#include <Flx_Group.h>
#include <FL/Fl_File_Browser.H>
#include <FL/Fl_Preferences.H>

class CWnd;

//-----------------------------------------------------------------------------
class CFileBrowser : public Flx_Group
{
public:
	CFileBrowser(int x, int y, int w, int h, CWnd* pwnd);
	virtual ~CFileBrowser(void);

protected:
	void OnSelectFile(Fl_Widget* pw, void* pd);

protected:
	bool is_dir(const char* szfile);
	void set_dir(const char* szfile);

protected:
	CWnd*		m_pWnd;
	Fl_File_Browser*	m_pfile;

	char	m_szdir[1024];	// current directory

private:
//	static Fl_Preferences	m_prefs;
};
