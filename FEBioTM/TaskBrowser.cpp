#include "stdafx.h"
#include "TaskBrowser.h"
#include "Wnd.h"
#include "Document.h"

//=============================================================================
CTaskTable::CTaskTable(int X, int Y, int W, int H, CWnd* pwnd) : Fl_Table_Row(X, Y, W, H), m_pWnd(pwnd)
{
    rows(0);				// how many rows
    row_header(0);          // enable row headers (along left)
    row_height_all(20);     // default height of rows
    row_resize(0);          // disable row resizing

    cols(2);				// how many columns
    col_header(1);          // enable column headers (along top)
    col_width(0,W-WIDTH);      // default width of columns
    col_width(1,WIDTH-4);      // default width of columns
    col_resize(0);          // enable column resizing

	// the one-and-only progress bar
	m_pg = new Fl_Progress(X+W/2, Y+H/2, 0, 0);
	m_pg->hide();
	m_pg->selection_color(FL_GREEN);
	m_nrow = -1;
	end();

	type(Fl_Table_Row::SELECT_SINGLE);
	selection_color(fl_rgb_color(216,216,255));

	color(FL_WHITE);
}

//-----------------------------------------------------------------------------
void CTaskTable::draw_cell(TableContext context, int ROW, int COL, int X, int Y, int W, int H)
{
	static char* szc[] = {"Path", "Status"};
	static char* szs[] = {"", "Modified", "Queued", "Running", "Completed", "Failed", "Cancelled"};

    switch ( context ) 
	{
     case CONTEXT_STARTPAGE:                   // before page is drawn..
        fl_font(FL_HELVETICA, 11);              // set the font for our drawing operations
        return; 
     case CONTEXT_COL_HEADER:                  // Draw column headers
		fl_push_clip(X,Y,W,H);
			fl_draw_box(FL_THIN_UP_BOX, X,Y,W,H, row_header_color());
			fl_color(FL_BLACK);
			fl_draw(szc[COL], X+5,Y,W-5,H, FL_ALIGN_LEFT);
		fl_pop_clip();
		return; 
     case CONTEXT_CELL:                        // Draw data in cells
		 {
			 bool b = (row_selected(ROW)==1);
			CDocument* pdoc = m_pWnd->GetDocument();
			CTask* pt = pdoc->GetTask(ROW);
			 fl_push_clip(X,Y,W,H);
			 // Draw cell bg
			 fl_color((b?selection_color():FL_WHITE)); fl_rectf(X,Y,W,H);
			 if ((ROW == m_nrow) && (COL == 1) && m_pg->visible()) return;
			 // Draw cell data
			 fl_color(FL_GRAY0); fl_draw((COL == 0? pt->GetFileName() : szs[pt->GetStatus()]), X+5,Y,W-5,H, (COL == 0 ? FL_ALIGN_LEFT : FL_ALIGN_CENTER));
			 // Draw box border
//			 fl_color(color()); fl_rect(X,Y,W,H);
			 fl_pop_clip();
		 }
		 return;
     default:
		 return;
    }
}

//-----------------------------------------------------------------------------
void CTaskTable::resize(int X, int Y, int W, int H)
{
	Fl_Table_Row::resize(X, Y, W, H);
    col_width(0,W-WIDTH);      // default width of columns
    col_width(1,WIDTH-4);      // default width of columns
}

//-----------------------------------------------------------------------------
void CTaskTable::show_progress(int nrow)
{
	int X, Y, W, H;
	m_nrow = nrow;
	find_cell(CONTEXT_CELL, nrow, 1, X, Y, W, H);
	m_pg->resize(X+1, Y+2, W-2, 18);
	m_pg->show();
}

//-----------------------------------------------------------------------------
void CTaskTable::hide_progress()
{
	if (m_pg->visible()) m_pg->hide();
}

//=============================================================================
CTaskBrowser::CTaskBrowser(int x, int y, int w, int h, CWnd* pwnd) : Flx_Group(x, y, w, h), m_pWnd(pwnd)
{
	begin();
	{
		m_pg = new CTaskTable(x, y, w, h, pwnd);
		pwnd->AddCallback(m_pg, (FLX_CALLBACK) &CWnd::OnSelectFile);
	}
	end();
	box(FL_NO_BOX);
	resizable(m_pg);
}

//-----------------------------------------------------------------------------
void CTaskBrowser::Update()
{
	CDocument* pdoc = m_pWnd->GetDocument();
	m_pg->rows(pdoc->Tasks());
	if (m_pg->rows() > 0) m_pg->select_row(0);
}

//-----------------------------------------------------------------------------
void CTaskBrowser::AddTask(CTask *pt)
{
	// add a row
	int N = m_pg->rows();
	m_pg->rows(N+1);
	m_pg->select_row(N);
}

//-----------------------------------------------------------------------------
void CTaskBrowser::RemoveTask(int n)
{
	int N = m_pg->rows();
	if (N > 0) m_pg->rows(N - 1);
}

//-----------------------------------------------------------------------------
void CTaskBrowser::SelectTask(int n)
{
	m_pg->select_row(n);
	m_pg->row_position(n);
}

//-----------------------------------------------------------------------------
int CTaskBrowser::SelectedTask()
{
	for (int i=0; i<m_pg->rows(); ++i)
		if (m_pg->row_selected(i) == 1) return i;
	return -1;
}

//-----------------------------------------------------------------------------
Fl_Progress* CTaskBrowser::TrackSelectedTask()
{
	int nrow = SelectedTask();
	m_pg->show_progress(nrow);
	return m_pg->GetProgressBar();
}

//-----------------------------------------------------------------------------
void CTaskBrowser::DoneTracking()
{
	m_pg->hide_progress();
}
