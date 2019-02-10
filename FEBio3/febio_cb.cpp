#include "stdafx.h"
#include "febio_cb.h"
#include <FEBioLib/FEBioModel.h>
#include <FEBioLib/version.h>
#include "console.h"
#include "Interrupt.h"
#include "breakpoint.h"
#include "FEBioApp.h"
#include <iostream>

//-----------------------------------------------------------------------------
// callback to update window title
bool update_console_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	FEBioModel& fem = static_cast<FEBioModel&>(*pfem);

	// get the number of steps
	int nsteps = fem.Steps();

	// calculate progress
	double starttime = fem.GetStartTime();
	double endtime = fem.GetEndTime();
	double f = 0.0;
	double ftime = fem.GetCurrentTime();
	if (endtime > 0.0) f = 100.f*(ftime - starttime) / (endtime - starttime);

	// check debug flag
	bool bdebug = fem.GetDebugFlag();

	// obtain a pointer to the console object. We'll use this to
	// set the title of the console window.
	Console* pShell = Console::GetHandle();

	char szvers[32] = {0};
#ifdef _DEBUG
	sprintf(szvers, "FEBio (DEBUG) %d.%d.%d", VERSION, SUBVERSION, SUBSUBVERSION);
#else
	sprintf(szvers, "FEBio %d.%d.%d", VERSION, SUBVERSION, SUBSUBVERSION);
#endif

	// print progress in title bar
	const char* szfile = fem.GetFileTitle();
	if (szfile == 0) szfile = "";

	if (nsteps > 1)
		pShell->SetTitle("(step %d/%d: %.f%%) %s - %s %s", fem.GetCurrentStepIndex() + 1, nsteps, f, szfile, szvers, (bdebug?"(debug mode)": ""));
	else
		pShell->SetTitle("(%.f%%) %s - %s %s", f, szfile, szvers, (bdebug?"(debug mode)": ""));

	// set progress (will print progress on task bar)
	if (nwhen == CB_SOLVED) pShell->SetProgress(100.0);
	else pShell->SetProgress(f);

	return true;
}

//-----------------------------------------------------------------------------
// break points cb
bool break_point_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	// get the current simulation time
	double t = pfem->GetTime().currentTime;

	// see if a break point was reached
	int bp = check_break(nwhen, t);
	if (bp >= 0)
	{
		std::cout << "breakpoint " << bp + 1 << " reached\n";

		// get a pointer to the app
		FEBioApp* app = FEBioApp::GetInstance();

		// process the commands
		if (app) app->ProcessCommands();
	}
	return true;
}

//-----------------------------------------------------------------------------
// callback for ctrl+c interruptions
bool interrupt_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	Interruption itr;
	if (itr.m_bsig)
	{
		itr.m_bsig = false;

		std::cout << "User interruption\n";

		// get a pointer to the app
		FEBioApp* app = FEBioApp::GetInstance();

		// process the commands
		if (app) app->ProcessCommands();
	}
	return true;
}
