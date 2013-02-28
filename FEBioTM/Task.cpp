#include "stdafx.h"
#include "Task.h"
#include <assert.h>
#include "MainApp.h"
#include "Document.h"

//-----------------------------------------------------------------------------
// initialize static CTask variables
CTask* CTask::m_prun = 0;

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
	m_pfile->select(0, m_pfile->length());
	m_pfile->remove_selection();	
	m_pfile->appendfile(m_szfile);
	SetStatus(READY);
}

//-----------------------------------------------------------------------------
// This callback will update the UI
void update_ui_cb(FEModel* pfem, void* pd)
{
	Progress* pp = (Progress*) pd;

	// get the number of steps
	int nsteps = pfem->Steps();

	// calculate progress
	double starttime = pfem->m_ftime0;
	double endtime = pfem->GetCurrentStep()->m_tend;
	double f = 100.f*(pfem->m_ftime - starttime) / (endtime - starttime);

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
