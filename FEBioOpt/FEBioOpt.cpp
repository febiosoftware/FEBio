#include "FEBioOpt.h"
#include "FEOptimizer.h"
#include <FECore/febio.h>
#include <FECore/Logfile.h>

//-----------------------------------------------------------------------------
// declared in dllmain.cpp
extern FEBioKernel* pFEBio;

static Logfile& GetLogfile() { return pFEBio->GetLogfile(); }

//-----------------------------------------------------------------------------
FEBioOpt::FEBioOpt(FEModel* pfem) : FEBioTask(pfem)
{

}

//-----------------------------------------------------------------------------
bool FEBioOpt::Run(const char* szfile)
{
	Logfile& log = GetLogfile();

	// TODO: the logfile has not been opened yet, so this will only get printed to the screen
	log.printbox("P A R A M E T E R   O P T I M I Z A T I O N   M O D U L E", "version 0.1");

	// create an optimizer object
	FEOptimizeData opt(*m_pfem);

	// read the data from the xml input file
	if (opt.Input(szfile) == false) return false;

	// do initialization
	if (opt.Init() == false) return false;

	// solve the problem
	bool bret = opt.Solve();

	if (bret)
		log.printf(" N O R M A L   T E R M I N A T I O N\n\n");
	else 
		log.printf(" E R R O R   T E R M I N A T I O N\n\n");

	return bret;
}
