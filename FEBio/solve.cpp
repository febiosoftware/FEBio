#include "stdafx.h"
#include <math.h>
#include "Timer.h"
#include "PlotFile.h"
#include "stack.h"
#include "fem.h"
#include "console.h"

//-----------------------------------------------------------------------------
//! This is the main solve method. This function loops over all analysis steps
//! and solves each one in turn. 
//! \sa FEAnalysis

bool FEM::Solve()
{
	// echo fem data to the logfile
	// we do this here (and not e.g. directly after input)
	// since the data can be changed after input, which is the case,
	// for instance, in the parameter optimization module
	EchoInput();

	// obtain a pointer to the console window. We'll use this to 
	// set the title of the window
	Console* pShell = Console::GetHandle();

	// Create timer to track total running time
	Timer	m_TotalTime;

	// start the total time tracker
	m_TotalTime.start();

	// convergence flag
	bool bconv = true;

	// intialize time
	m_ftime = 0;

	// loop over all analysis steps
	for (int nstep=0; nstep < m_Step.size(); ++nstep)
	{
		// set the current analysis step
		m_nStep = nstep;
		m_pStep = &m_Step[nstep];

		// intitialize step data
		if (m_pStep->Init() == false)
		{
			bconv = false;
			break;
		}

		// solve the analaysis step
		bconv = m_pStep->Solve();

		// break if the step has failed
		if (bconv == false) break;

		// wrap it up
		m_pStep->Finish();
	}

	// close the plot file
	m_plot.Close();

	// stop total time tracker
	m_TotalTime.stop();

	// get and print elapsed time
	char sztime[64];

	m_TotalTime.time_str(sztime);
	m_log.printf("\n Elapsed time : %s\n\n", sztime);

	if (bconv)
	{
		m_log.printf("\n N O R M A L   T E R M I N A T I O N\n\n");
	}
	else
	{
		m_log.printf("\n E R R O R   T E R M I N A T I O N\n\n");
	}

	pShell->SetTitle("(%s) %s - FEBio", (bconv?"NT":"ET"), m_szfile_title);

	// We're done !
	return bconv;
}
