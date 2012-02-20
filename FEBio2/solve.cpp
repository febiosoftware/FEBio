#include "stdafx.h"
#include <math.h>
#include "FEBioLib/Timer.h"
#include "fem.h"
#include "FEBioProgress.h"
#include "console.h"
#include "FEBioLib/log.h"

//-----------------------------------------------------------------------------
bool solve(FEM& fem, const char* szfile)
{
	// read input data
	if (fem.Input(szfile) == false) return false;

	// initialize and check data
	if (fem.Init() == false) return false;

	// create a progress tracker
	FEBioProgress prg(fem);

	// run the analysis
	bool bconv = fem.Solve(prg);

	// obtain a pointer to the console window. 
	// We'll use this to set the title of the window
	Console* pShell = Console::GetHandle();
	pShell->SetTitle("(%s) %s - FEBio", (bconv?"NT":"ET"), fem.GetFileTitle());

	return bconv;
}
