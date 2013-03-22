#include "stdafx.h"
#include "Task.h"
#include <assert.h>
#include "MainApp.h"
#include "Document.h"

//-----------------------------------------------------------------------------
// This function will be called when the text buffer is modified
void update_style(int pos, int nInserted, int nDeleted, int nRestyled, const char* deletedText, void* cbArg)
{
	CTask* pt = (CTask*) cbArg;
	pt->UpdateStyle(pos, nInserted, nDeleted, nRestyled, deletedText);
}

//-----------------------------------------------------------------------------
// initialize static CTask variables
CTask* CTask::m_prun = 0;

//-----------------------------------------------------------------------------
CTask::CTask()
{
	m_szfile[0] = 0; m_pFile = 0; m_pStyle = 0; m_nstatus = READY;

	// default FEBio options
	m_bdebug = false;
}

//-----------------------------------------------------------------------------
CTask::~CTask()
{ 
	m_pFile->remove_modify_callback(update_style, (void*)this);
	delete m_pFile; 
	delete m_pStyle;
}

//-----------------------------------------------------------------------------
void CTask::SetFileName(const char* szfile)
{
	m_szfile[0] = 0;
	int l = strlen(szfile)+1;
	assert((l>1) && (l<MAX_FILE));
	if ((l > 1) && (l<MAX_FILE)) strncpy(m_szfile, szfile, l);
}

//-----------------------------------------------------------------------------
const char* CTask::GetFileTitle()
{
	char* c1 = strrchr(m_szfile, '\\');
	char* c2 = strrchr(m_szfile, '/');
	if ((c1 == 0) && (c2 == 0)) return m_szfile;
	if (c1 == 0) return c2+1;
	if (c2 == 0) return c1+1;
	if (c1 < c2) return c2+1; else return c1+1;
}

//-----------------------------------------------------------------------------
void CTask::GetFilePath(char* szpath)
{
	char* c1 = strrchr(m_szfile, '\\');
	char* c2 = strrchr(m_szfile, '/');
	if ((c1 == 0) && (c2 == 0)) strcpy(szpath, m_szfile);
	if ((c1 == 0) || (c2 > c1)) strncpy(szpath, m_szfile, c2 - m_szfile);
	if ((c2 == 0) || (c1 > c2)) strncpy(szpath, m_szfile, c1 - m_szfile);
}

//-----------------------------------------------------------------------------
void CTask::Revert()
{
	// clear the file buffer
	m_pFile->select(0, m_pFile->length());
	m_pFile->remove_selection();	
	m_pFile->appendfile(m_szfile);
	SetStatus(READY);
}

//-----------------------------------------------------------------------------
bool is_comment_start(char* sz)
{
	if (sz[0] != '<') return false;
	if (sz[1] != '!') return false;
	if (sz[2] != '-') return false;
	if (sz[3] != '-') return false;
	return true;
}

//-----------------------------------------------------------------------------
bool is_comment_end(char* sz)
{
	if (sz[ 0] != '>') return false;
	if (sz[-1] != '-') return false;
	if (sz[-2] != '-') return false;
	return true;
}

//-----------------------------------------------------------------------------
void format_style(char* cs, char* cd, int l)
{
	char style = 'A';
	int nkey = 0;
	int ncomm = 0;
	for (int i=0; i<l; ++i, ++cd, ++cs)
	{
		if (ncomm == 0)
		{
			if ((*cs == '<') && is_comment_start(cs))
			{
				*cd = 'E';
				ncomm = 1;
			}
			else
			{
				switch (*cs)
				{
				case '<' : { *cd =  'B'; style = 'B'; nkey = 1; } break;
				case '>' : { *cd =  'B'; style = 'A'; nkey = 0; } break;
				case '"' : { *cd =  'C'; style = (style=='C'?'A':'C'); } break;
				case ' ' : { *cd =  'A'; style = (style=='C'?'C':(nkey==1?'D':'A')); } break;
				case '=' : { *cd =  'A'; } break;
				case '\n': { *cd = '\n'; style = 'A'; } break;
				case '/' : { *cd = 'B'; } break;
				default:
					*cd = style;
				}
			}
		}
		else
		{
			*cd = 'E';
			if ((*cs == '>') && is_comment_end(cs))
			{
				style = 'A';
				ncomm = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void CTask::SetTextBuffer(Fl_Text_Buffer* pb)
{
	// store the text buffer
	m_pFile = pb;
	m_pFile->tab_distance(4);

	// create a style buffer
	assert(m_pStyle == 0);
	m_pStyle = new Fl_Text_Buffer(m_pFile->length());

	// format the text
	char *style = new char[m_pFile->length() + 1];
	char *text = m_pFile->text();

	memset(style, 'A', m_pFile->length());
	style[m_pFile->length()] = '\0';

	// process the buffer
	format_style(text, style, m_pFile->length());

	// copy the style buffer
	m_pStyle->text(style);
	delete[] style;
	free(text);

	m_pFile->add_modify_callback(update_style, (void*)this);
}

//-----------------------------------------------------------------------------
void CTask::UpdateStyle(int npos, int nInserted, int nDeleted, int nRestyled, const char* deletedText)
{
  // If this is just a selection change, just unselect the style buffer...
  if (nInserted == 0 && nDeleted == 0) {
    m_pStyle->unselect();
    return;
  }

  // Track changes in the text buffer...
  if (nInserted > 0) {
    // Insert characters into the style buffer...
    char* style = new char[nInserted + 1];
    memset(style, 'A', nInserted);
    style[nInserted] = '\0';

    m_pStyle->replace(npos, npos + nDeleted, style);
    delete [] style;
  } else {
    // Just delete characters in the style buffer...
    m_pStyle->remove(npos, npos + nDeleted);
  }

  // Select the area that was just updated to avoid unnecessary callbacks...
  m_pStyle->select(npos, npos + nInserted - nDeleted);

  int start = m_pFile->line_start(npos);
  int end   = m_pFile->line_end(npos + nInserted);
  char* text  = m_pFile->text_range(start, end);
  char* style = m_pStyle->text_range(start, end);
  char last = (start==end?0:style[end - start - 1]);

  format_style(text, style, end - start);
  m_pStyle->replace(start, end, style);

  if (start==end || last != style[end - start - 1]) {
    // Either the user deleted some text, or the last character
    // on the line changed styles, so reparse the
    // remainder of the buffer...
    free(text);
    free(style);

    end   = m_pFile->length();
    text  = m_pFile->text_range(start, end);
    style = m_pStyle->text_range(start, end);

    format_style(text, style, end - start);

    m_pStyle->replace(start, end, style);
  }

  free(text);
  free(style);
}

//-----------------------------------------------------------------------------
// This callback will update the UI
void update_ui_cb(FEModel* pfem, void* pd)
{
	Progress* pp = (Progress*) pd;

	// get the number of steps
	int nsteps = pfem->Steps();
	int nstepi = pfem->m_nStep;
	double pct = 100.0/nsteps;

	// calculate progress
	double starttime = pfem->m_ftime0;
	double endtime = pfem->GetCurrentStep()->m_tend;
	double f = pct*nstepi + pct*(pfem->m_ftime - starttime) / (endtime - starttime);

	// set the progress (will also update UI)
	pp->SetProgress(f);
}

//-----------------------------------------------------------------------------
void CTask::Run(Progress& prg)
{
	// set the status to running
	SetStatus(RUNNING);

	// mark this task as the running task
	assert(m_prun == 0);
	m_prun = this;

	// setup the FE problem
	FEM fem(this);

	// set the callback function
	fem.AddCallback(update_ui_cb, &prg);

	// set the default output file names
	char szbase[1024] = {0}, szfile[1024] = {0};
	strcpy(szbase, GetFileName());
	char* ch = strrchr(szbase, '.'); assert(ch);
	if (ch) *ch = 0;
	sprintf(szfile, "%s.log", szbase); fem.SetLogFilename(szfile);
	sprintf(szfile, "%s.plt", szbase); fem.SetPlotFilename(szfile);
	sprintf(szfile, "%s.dmp", szbase); fem.SetDumpFilename(szfile);
	fem.SetInputFilename(GetFileName());

	// set command line options
	if (m_bdebug) fem.SetDebugFlag(true);

	// load the data from file
	if (fem.Input(GetFileName()) == false)
	{
		SetStatus(CTask::FAILED);
		m_prun = 0;
		return;
	}

	// initialize FE data
	if (fem.Init() == false) 
	{
		SetStatus(CTask::FAILED);
		m_prun = 0;
		return;
	}

	// run the problem
	bool bret = fem.Solve();

	// set the final status
	// Note that the user could have cancelled this task, so
	// we need to check whether the status is still RUNNING
	if (GetStatus() == CTask::RUNNING) SetStatus(bret?CTask::COMPLETED:CTask::FAILED);

	// reset running task
	m_prun = 0;
}
