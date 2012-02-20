#include "stdafx.h"
#include "FEBioProgress.h"
#include "console.h"
#include "fem.h"

//-----------------------------------------------------------------------------
void FEBioProgress::SetProgress(double f)
{
	FEM& fem = dynamic_cast<FEM&>(m_fem);

	// get the number of steps
	int nsteps = m_fem.Steps();

	// check debug flag
	bool bdebug = m_fem.GetDebugFlag();

	// obtain a pointer to the console object. We'll use this to
	// set the title of the console window.
	Console* pShell = Console::GetHandle();

	// print progress in title bar
	if (nsteps > 1)
		pShell->SetTitle("(step %d/%d: %.f%%) %s - %s", m_fem.m_nStep+1, nsteps, f, fem.GetFileTitle(), (bdebug?"FEBio (debug mode)": "FEBio"));
	else
		pShell->SetTitle("(%.f%%) %s - %s", f, fem.GetFileTitle(), (bdebug?"FEBio (debug mode)": "FEBio"));
}
