#pragma once
#include <FL/Fl_Text_Display.H>

class Progress;

//-----------------------------------------------------------------------------
// The task defines a problem to be run with FEBio. It has several states:
// - READY = File is read in and awaits to be scheduled.
// - MODIFIED = File has been modified. Will be saved automatically when run.
// - QUEUED = File is scheduled to be run.
// - RUNNING = File is being run.
// - COMPLETED = File has run and terminated successfully.
// - FAILED = File has run but did not terminate successfully.
// - CANCELLED = Run was cancelled by the user before terminating.
//
// Each task also manages its own input file through the Fl_Text_Buffer
// reference. 
//
// Only one task can be run at the same time. This limitation must be
// imposed because there is only one output window where the output from
// the run will be directed to. The task that is currently running is
// store in the CTask::m_prun member.
//
class CTask
{
	// max file name
	enum {MAX_FILE = 512};

public:
	// status values
	// note that at any given time, only one task may have its status set to RUNNING
	enum { READY, MODIFIED, QUEUED, RUNNING, COMPLETED, FAILED, CANCELLED };

public:
	// con-/de-structor
	CTask();
	~CTask();

	// file name handling
	void SetFileName(const char* szfile);
	const char* GetFileName() { return m_szfile; }
	const char* GetFileTitle();
	void GetFilePath(char* sz);

	// get/set the input file text buffer
	void SetTextBuffer(Fl_Text_Buffer* pb) { m_pfile = pb; }
	Fl_Text_Buffer* GetTextBuffer() { return m_pfile; }

	// get/set the task status
	void SetStatus(int n) { m_nstatus = n; }
	int GetStatus() { return m_nstatus; }

	// save the task input file
	void Save() { if (m_nstatus == MODIFIED) { m_pfile->savefile(m_szfile); SetStatus(READY);} }
	void Save(const char* szfile) { SetFileName(szfile); Save(); }

	// reload the file
	void Revert();

	// get/set progress indicator
	void SetProgress(float f) { m_prg = f; }
	float GetProgress() { return m_prg; }

	// run this task
	void Run(Progress& prg);
	static CTask* GetRunningTask() { return m_prun; }

protected:
	char			m_szfile[MAX_FILE];		//!< file name
	Fl_Text_Buffer*	m_pfile;				//!< text buffer for editing
	int				m_nstatus;				//!< status
	float			m_prg;					//!< progress indicator

public: // FEBio command line options
	bool	m_bdebug;	//!< debug mode

private:
	static CTask*	m_prun;	// this is the task that is running
};
