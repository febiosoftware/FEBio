#pragma once
#include <Flx_Group.h>
#include <FL/Fl_Table_Row.H>
#include <FL/Fl_Progress.H>

class CWnd;
class CTask;

//-----------------------------------------------------------------------------
class CTaskTable : public Fl_Table_Row
{
	enum { WIDTH = 100, HEIGHT = 10 };

public:
	CTaskTable(int X, int Y, int W, int H, CWnd* pwnd);

public:
	void draw_cell(TableContext context, int ROW, int COL, int X, int Y, int W, int H);

	void resize(int X, int Y, int W, int H);

	int handle(int nevent);

	void show_progress(int nrow);
	void hide_progress();

	Fl_Progress* GetProgressBar() { return m_pg; }

public:
	void SelectNext();
	void SelectPrev();

protected:
	CWnd*	m_pWnd;
	Fl_Progress*	m_pg;
	int				m_nrow; // row where progress bar is visible
	int				m_nfocus;	// row that has focus
};

//-----------------------------------------------------------------------------
class CTaskBrowser : public Flx_Group
{
public:
	CTaskBrowser(int x, int y, int w, int h, CWnd* pwnd);

	void Update();

	void AddTask(CTask *pt);

	void RemoveTask(int n);

	void SelectTask(int n);

	void SelectAll(int n = 1);

	int SelectedTask();

	//! return the number of selected tasks
	int SelectedTasks();

	//! see if a task is selected
	bool IsTaskSelected(int n);

	Fl_Progress* TrackTask(int nrow);

	void DoneTracking();

protected:
	CWnd*			m_pWnd;
	CTaskTable*		m_pg;
};
