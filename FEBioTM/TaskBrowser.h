#pragma once
#include <Flx_Group.h>
#include <FL/Fl_Table_Row.H>
#include <FL/Fl_Progress.H>

class CWnd;
class CTask;

//-----------------------------------------------------------------------------
class CTaskTable : public Fl_Table_Row
{
	enum { WIDTH = 100 };

public:
	CTaskTable(int X, int Y, int W, int H, CWnd* pwnd);

public:
	void draw_cell(TableContext context, int ROW, int COL, int X, int Y, int W, int H);

	void resize(int X, int Y, int W, int H);

	void show_progress(int nrow);
	void hide_progress();

	Fl_Progress* GetProgressBar() { return m_pg; }

protected:
	CWnd*	m_pWnd;
	Fl_Progress*	m_pg;
	int				m_nrow; // row where progress bar is visible
};

//-----------------------------------------------------------------------------
class CTaskBrowser : public Flx_Group
{
public:
	CTaskBrowser(int x, int y, int w, int h, CWnd* pwnd);

	void Update();

	void AddTask(CTask *pt);

	void RemoveTask(int n);

	int SelectedTask();

	Fl_Progress* TrackSelectedTask();

	void DoneTracking();

protected:
	CWnd*			m_pWnd;
	CTaskTable*		m_pg;
};
