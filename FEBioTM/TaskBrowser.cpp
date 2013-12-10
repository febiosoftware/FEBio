#include "stdafx.h"
#include "TaskBrowser.h"
#include "Wnd.h"
#include "Document.h"

//=============================================================================
CTaskTable::CTaskTable(int X, int Y, int W, int H, CWnd* pwnd) : Fl_Table_Row(X, Y, W, H), m_pWnd(pwnd)
{
	m_pg = 0;

    rows(0);				// how many rows
    row_header(0);          // enable row headers (along left)
    row_resize(0);          // disable row resizing

    cols(4);				// how many columns
    col_header(1);          // enable column headers (along top)
    col_width(0,W-3*WIDTH);      // default width of columns
    col_width(1,WIDTH);      // default width of columns
    col_width(2,WIDTH);      // default width of columns
    col_width(3,WIDTH-4);      // default width of columns
    col_resize(0);          // enable column resizing

	// the one-and-only progress bar
	begin();
	{
		m_pg = new Fl_Progress(X+W/2, Y+H/2, 10, 10);
		m_pg->hide();
		m_pg->selection_color(FL_GREEN);
		m_nrow = -1;
		m_pg->box(FL_BORDER_BOX);
	}
	end();

	type(Fl_Table_Row::SELECT_MULTI);
	selection_color(fl_rgb_color(216,216,255));

	color(FL_WHITE);
}

//-----------------------------------------------------------------------------
void CTaskTable::draw_cell(TableContext context, int ROW, int COL, int X, int Y, int W, int H)
{
	static char* szc[] = {"Path", "File Size", "Status", "Run time"};
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
	 case CONTEXT_RC_RESIZE:
		if (m_pg && m_pg->visible()) 
		{
			if ((ROW == 0)&&(COL == 0))
			{
				find_cell(CONTEXT_TABLE, m_nrow, 2, X, Y, W, H);
				m_pg->resize(X, Y, W, H);
			}
		}
		break;
     case CONTEXT_CELL:                        // Draw data in cells
		 {
			 if ((ROW == m_nrow) && (COL == 2) && m_pg->visible()) return;

			 bool b = (row_selected(ROW)==1);
			CDocument* pdoc = m_pWnd->GetDocument();
			CTask* pt = pdoc->GetSession().GetTask(ROW);
			 fl_push_clip(X,Y,W,H);
			 // Draw cell bg
			 fl_color((b?selection_color():FL_WHITE)); fl_rectf(X,Y,W,H);
			 // Draw cell data
			 fl_color(FL_GRAY0); 
			 char szl[256] = {0};
			 float g;
			switch (COL)
			{
			case 0: fl_draw(pt->GetFileName(), X+5,Y,W-5,H, FL_ALIGN_LEFT); break;
			case 1: g = pt->FileSize() / 1024.f; sprintf(szl, "%.2f KB", g); fl_draw(szl, X+5,Y,W-10,H, FL_ALIGN_RIGHT); break;
			case 2: fl_draw(szs[pt->GetStatus()], X+5,Y,W-5,H, FL_ALIGN_CENTER); break;
			case 3: 
				{
					int nstatus = pt->GetStatus();
					if ((nstatus == CTask::COMPLETED) || 
						(nstatus == CTask::RUNNING  ) ||
						(nstatus == CTask::CANCELLED))
					{
						char sz[32];
						CTask::STATS& s = pt->m_stats;
						sprintf(sz, "%d:%02d:%02d", s.nhour, s.nmin, s.nsec);
						fl_draw(sz, X, Y, W-5, H, FL_ALIGN_RIGHT);
					}
				}
				break;
			}
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
    col_width(0,W-3*WIDTH);      // default width of columns
    col_width(1,WIDTH);			// default width of columns
    col_width(2,WIDTH);			// default width of columns
    col_width(3,WIDTH-4);      // default width of columns
	if (m_pg->visible()) show_progress(m_nrow);
}

//-----------------------------------------------------------------------------
void CTaskTable::show_progress(int nrow)
{
	int X, Y, W, H;
	m_nrow = nrow;
	find_cell(CONTEXT_CELL, nrow, 2, X, Y, W, H);
	m_pg->resize(X, Y, W, H);
	if (m_pg->visible() == 0)
	{
		m_pg->label("0%");
		m_pg->show();
	}
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
	m_pg->rows(pdoc->GetSession().Tasks());
	if (m_pg->rows() > 0)
	{
		m_pg->row_height_all(20); // default height of rows
		m_pg->select_row(0);
	}
}

//-----------------------------------------------------------------------------
void CTaskBrowser::AddTask(CTask *pt)
{
	// add a row
	int N = m_pg->rows();
	m_pg->rows(N+1);
    m_pg->row_height_all(20); // default height of rows
	m_pg->select_all_rows(0);
	m_pg->select_row(N);
}

//-----------------------------------------------------------------------------
void CTaskBrowser::RemoveTask(int n)
{
	int N = m_pg->rows();
	if (N > 0)
	{
		m_pg->rows(N - 1);
	    m_pg->row_height_all(20); // default height of rows
	}
}

//-----------------------------------------------------------------------------
void CTaskBrowser::SelectAll(int n)
{
	m_pg->select_all_rows(n);
}

//-----------------------------------------------------------------------------
void CTaskBrowser::SelectTask(int n)
{
	m_pg->select_row(n);
	m_pg->row_position(n);
}

//-----------------------------------------------------------------------------
int CTaskBrowser::SelectedTasks()
{
	int n = 0;
	for (int i=0; i<m_pg->rows(); ++i)
		if (m_pg->row_selected(i) == 1) n++;
	return n;
}

//-----------------------------------------------------------------------------
bool CTaskBrowser::IsTaskSelected(int n)
{
	return (m_pg->row_selected(n) == 1);
}

//-----------------------------------------------------------------------------
int CTaskBrowser::SelectedTask()
{
	for (int i=0; i<m_pg->rows(); ++i)
		if (m_pg->row_selected(i) == 1) return i;
	return -1;
}

//-----------------------------------------------------------------------------
Fl_Progress* CTaskBrowser::TrackTask(int nrow)
{
	m_pg->show_progress(nrow);
	return m_pg->GetProgressBar();
}

//-----------------------------------------------------------------------------
void CTaskBrowser::DoneTracking()
{
	m_pg->hide_progress();
}
