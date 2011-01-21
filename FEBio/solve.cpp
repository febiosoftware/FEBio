#include "stdafx.h"
#include <math.h>
#include "Timer.h"
#include "PlotFile.h"
#include "fem.h"
#include "console.h"
#include "log.h"

//-----------------------------------------------------------------------------
// echo the input data to the log file
extern void echo_input(FEM& fem);

//-----------------------------------------------------------------------------
bool solve(FEM& fem, const char* szfile)
{
	// read input data
	if (fem.Input(szfile) == false) return false;

	// initialize and check data
	if (fem.Init() == false) return false;

	// run the analysis
	return fem.Solve();
}

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
	echo_input(*this);

	// obtain a pointer to the console window. We'll use this to 
	// set the title of the window
	Console* pShell = Console::GetHandle();

	// start the total time tracker
	m_TotalTime.start();

	// convergence flag
	bool bconv = true;

	// loop over all analysis steps
	// Note that we don't necessarily from step 0.
	// This is because the user could have restarted
	// the analysis. 
	for (size_t nstep=m_nStep; nstep < m_Step.size(); ++nstep)
	{
		// set the current analysis step
		m_nStep = nstep;
		m_pStep = m_Step[nstep];

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
	if (m_plot) m_plot->Close();

	// stop total time tracker
	m_TotalTime.stop();

	// get and print elapsed time
	char sztime[64];

	// get the logfile
	Logfile& log = GetLogfile();

	m_TotalTime.time_str(sztime);
	log.printf("\n Elapsed time : %s\n\n", sztime);

	if (bconv)
	{
		log.printf("\n N O R M A L   T E R M I N A T I O N\n\n");
	}
	else
	{
		log.printf("\n E R R O R   T E R M I N A T I O N\n\n");
	}

	pShell->SetTitle("(%s) %s - FEBio", (bconv?"NT":"ET"), m_szfile_title);

	// We're done !
	return bconv;
}
