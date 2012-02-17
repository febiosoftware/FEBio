#pragma once
#include <Flx_Group.h>
#include <FL/Fl_Table_Row.H>

class CWnd;
class CTask;

//-----------------------------------------------------------------------------
class CTaskTable : public Fl_Table_Row
{
public:
	CTaskTable(int X, int Y, int W, int H, CWnd* pwnd);

public:
	void draw_cell(TableContext context, int ROW, int COL, int X, int Y, int W, int H);

	void resize(int X, int Y, int W, int H);

protected:
	CWnd*	m_pWnd;
};

//-----------------------------------------------------------------------------
class CTaskBrowser : public Flx_Group
{
public:
	CTaskBrowser(int x, int y, int w, int h, CWnd* pwnd);

	void Update();

	void AddTask(CTask *pt);

	int SelectedTask();

protected:
	CWnd*			m_pWnd;
	CTaskTable*		m_pg;
};
