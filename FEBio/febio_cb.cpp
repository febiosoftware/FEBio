/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



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
	if (nwhen != CB_INIT)
	{
		double ftime = fem.GetCurrentTime();
		if (endtime != starttime) f = (ftime - starttime) / (endtime - starttime);
		else
		{
			// this only happens (I think) when the model is solved
			f = 1.0;
		}
	}

	double pct = 100.0*f;

	// check debug flag
	int ndebug = fem.GetDebugLevel();

	// obtain a pointer to the console object. We'll use this to
	// set the title of the console window.
	Console* pShell = Console::GetHandle();

	
	char* szver = febio::getVersionString();

	char szvers[64] = {0};
#ifdef _DEBUG
	sprintf(szvers, "FEBio (DEBUG) %s", szver);
#else
	sprintf(szvers, "FEBio %s", szver);
#endif

	// print progress in title bar
	const std::string& sfile = fem.GetFileTitle();
	const char* szfile = sfile.c_str();

	if (nsteps > 1)
		pShell->SetTitle("(step %d/%d: %.f%%) %s - %s %s", fem.GetCurrentStepIndex() + 1, nsteps, pct, szfile, szvers, (ndebug?"(debug mode)": ""));
	else
		pShell->SetTitle("(%.f%%) %s - %s %s", pct, szfile, szvers, (ndebug?"(debug mode)": ""));

	if (nsteps > 1)
	{
		int step = fem.GetCurrentStepIndex();
		double w = (step + f) / nsteps;
		pct = 100.0*w;
	}

	// set progress (will print progress on task bar)
	if (nwhen == CB_SOLVED) pShell->SetProgress(100.0);
	else pShell->SetProgress(pct);

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
	// Do not process this when we are writing or reading the dump file since
	// this may cause problems
	if ((nwhen == CB_SERIALIZE_LOAD) || (nwhen == CB_SERIALIZE_SAVE)) return true;

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
